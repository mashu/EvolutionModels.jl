module FASTXExt

using FASTX
using BioSequences
using Optim
using LinearAlgebra

import EvolutionModels: Model, sequence_likelihood, read_alignment, compute_distances, is_transition
import EvolutionModels: JC69Model, HKY85Model, GTRModel
import EvolutionModels: num_parameters, STANDARD_DNA
export LabeledAlignment
import EvolutionModels: WAGModel, LGModel, STANDARD_AA

# Custom error types for better error handling
struct AlignmentError <: Exception
    msg::String
end

struct SequenceConversionError <: Exception
    msg::String
end

"""
    LabeledAlignment{T<:LongSequence}

Structure to hold aligned sequences with their labels.

# Fields
- `sequences::Vector{T}`: Vector of aligned sequences
- `labels::Vector{String}`: Vector of sequence identifiers
- `length::Int`: Length of the alignment
"""
struct LabeledAlignment{T<:LongSequence}
    sequences::Vector{T}
    labels::Vector{String}
    length::Int

    function LabeledAlignment{T}(sequences::Vector{T}, labels::Vector{String}) where T<:LongSequence
        if isempty(sequences)
            throw(AlignmentError("Empty sequence alignment"))
        end

        if length(sequences) != length(labels)
            throw(AlignmentError("Number of sequences ($(length(sequences))) does not match number of labels ($(length(labels)))"))
        end

        aln_length = length(first(sequences))
        if !all(length(seq) == aln_length for seq in sequences)
            throw(AlignmentError("Sequences must be of equal length"))
        end

        new(sequences, labels, aln_length)
    end
end

"""
    read_alignment(fasta_file::String) -> LabeledAlignment

Read multiple sequence alignment from FASTA file.

# Arguments
- `fasta_file::String`: Path to FASTA format alignment file

# Returns
- `LabeledAlignment`: Structure containing sequences and their labels

# Throws
- `SystemError`: If file cannot be opened or read
- `AlignmentError`: If alignment is empty or sequences are not of equal length
- `SequenceConversionError`: If sequences cannot be converted to appropriate type
"""
function read_alignment(fasta_file::String)
    !isfile(fasta_file) && throw(SystemError("File not found: $fasta_file"))

    try
        reader = FASTA.Reader(open(fasta_file))

        # Pre-allocate vectors
        sequences = LongDNA{4}[]
        labels = String[]

        for record in reader
            try
                push!(sequences, LongDNA{4}(sequence(record)))
                push!(labels, identifier(record))
            catch e
                throw(SequenceConversionError("Failed to convert sequence '$(identifier(record))': $(e.msg)"))
            end
        end

        close(reader)
        return LabeledAlignment{LongDNA{4}}(sequences, labels)

    catch e
        if e isa SystemError || e isa AlignmentError || e isa SequenceConversionError
            rethrow(e)
        else
            throw(AlignmentError("Failed to read alignment: $(e.msg)"))
        end
    end
end

"""
    compute_distances(model::Model, aln::LabeledAlignment;
                     max_scale::Float64=100.0) -> NamedTuple

Compute pairwise distances between all sequences in the alignment using maximum likelihood.

# Arguments
- `model::Model`: Evolutionary model (JC69, HKY85, GTR, WAG, LG)
- `aln::LabeledAlignment`: Aligned sequences with labels
- `max_scale::Float64`: Maximum distance to consider (default=100.0)

# Returns
NamedTuple containing:
- `distances`: Distance matrix
- `labels`: Sequence labels
- `model`: Model name
"""
function compute_distances(
    model::Model,
    aln::LabeledAlignment;
    max_scale::Float64=100.0
)
    n = length(aln.sequences)
    D = zeros(n, n)
    model_name = get(model.params, :model, nothing)

    Threads.@threads for i in 1:n
        for j in i+1:n
            # Define negative log-likelihood function for optimization
            function neg_log_likelihood(t::Float64)
                if t < 0 || t > max_scale
                    return Inf
                end
                try
                    return -sequence_likelihood(model, aln.sequences[i], aln.sequences[j], t)
                catch
                    return Inf
                end
            end

            # Count valid positions for basic distance initialization
            total = 0
            diff = 0
            for (x, y) in zip(aln.sequences[i], aln.sequences[j])
                if (model_name in ("WAG", "LG") && x in STANDARD_AA && y in STANDARD_AA) ||
                   (model_name in ("JC69", "HKY85", "GTR") && x in STANDARD_DNA && y in STANDARD_DNA)
                    total += 1
                    if x != y
                        diff += 1
                    end
                end
            end

            if total == 0
                D[i,j] = D[j,i] = NaN
                continue
            end

            # Compute proportion of differences (p-distance)
            p_dist = diff / total
            
            # Handle identical sequences - distance is 0
            if diff == 0
                D[i,j] = D[j,i] = 0.0
                continue
            end
            
            # Initial estimate based on Jukes-Cantor correction as heuristic
            # For JC69: d = -3/4 * ln(1 - 4p/3) when p < 0.75
            # This provides a reasonable starting point for any model
            if p_dist < 0.75
                t_init = -0.75 * log(1.0 - 4.0 * p_dist / 3.0)
            else
                # For highly divergent sequences, use a larger estimate
                t_init = 2.0  # Roughly 2 substitutions per site
            end
            
            # Set search bounds using t_init as guidance
            # Lower bound: small positive to avoid log(0) issues
            # Upper bound: expand beyond t_init but cap at max_scale
            t_lower = 1e-10
            t_upper = min(max(t_init * 5.0, 1.0), max_scale)

            # Run optimization using Brent's method on informed interval
            res = try
                optimize(neg_log_likelihood, t_lower, t_upper, Brent())
            catch e
                # If first attempt fails, try full interval
                try
                    optimize(neg_log_likelihood, t_lower, max_scale, Brent())
                catch e2
                    @warn "Optimization failed for sequences $(aln.labels[i]) and $(aln.labels[j]): $(sprint(showerror, e2))"
                    nothing
                end
            end

            # Extract result, fall back to t_init if optimization didn't converge
            if res !== nothing && Optim.converged(res)
                D[i,j] = D[j,i] = Optim.minimizer(res)
            else
                # Use JC-corrected estimate as fallback
                D[i,j] = D[j,i] = t_init
                @warn "Using initial estimate for $(aln.labels[i]) vs $(aln.labels[j]): t=$t_init"
            end
        end
    end

    return (
        distances=D,
        labels=aln.labels,
        model=model_name
    )
end

"""
    num_parameters(model::Model) -> Int

Helper function to determine number of parameters in the model.
"""
function num_parameters(model::Model)
    model_name = get(model.params, :model, nothing)

    if model_name == "JC69"
        return 1  # Only evolutionary distance
    elseif model_name == "HKY85"
        return 2  # Distance and κ
    elseif model_name == "GTR"
        n = length(model.seq_type)
        return 1 + (n * (n-1))÷2  # Distance and symmetric rate matrix parameters
    elseif model_name == "WAG" || model_name == "LG"
        return 1  # Only evolutionary distance for empirical models
    else
        error("Unknown model type: $(model_name)")
    end
end

"""
    print_distance_matrix(result::NamedTuple; digits::Int=4)

Pretty print a labeled distance matrix.

# Arguments
- `result::NamedTuple`: Output from compute_distances
- `digits::Int`: Number of decimal places to display (default=4)
"""
function print_distance_matrix(result::NamedTuple; digits::Int=4)
    D, labels = result.distances, result.labels
    n = length(labels)

    # Calculate padding based on label lengths
    label_pad = maximum(length.(labels)) + 2
    num_pad = digits + 6

    # Print header row
    print(" " ^ label_pad)
    for label in labels
        print(rpad(label[1:min(end,20)], num_pad))
    end
    println()

    # Print matrix with row labels
    for i in 1:n
        print(rpad(labels[i][1:min(end,20)], label_pad))
        for j in 1:n
            val = D[i,j]
            str = isnan(val) ? "NA" : "$(round(val; digits=digits))"
            print(rpad(str, num_pad))
        end
        println()
    end
end

end # module