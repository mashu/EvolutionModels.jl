module FASTXExt

using FASTX
using BioSequences
using Optim

import EvolutionModels: Model, sequence_likelihood, read_alignment, compute_distances, is_transition
import EvolutionModels: JC69Model, HKY85Model, GTRModel
import EvolutionModels: num_parameters, STANDARD_DNA
export LabeledAlignment

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

# Example
```julia
aln = read_alignment("alignment.fasta")
println("Loaded \$(length(aln.sequences)) sequences of length \$(aln.length)")
```
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
    compute_analytical_distance(model::Model, seq1::LongSequence, seq2::LongSequence) -> Float64

Compute analytical distance between two sequences based on the model type.
Returns NaN if the distance cannot be computed.
"""
function compute_analytical_distance(model::Model, seq1::LongSequence, seq2::LongSequence)
    # Count sites and patterns, ignoring non-standard nucleotides
    total = 0
    S = 0  # transitions
    V = 0  # transversions

    for (x, y) in zip(seq1, seq2)
        if x in STANDARD_DNA && y in STANDARD_DNA
            total += 1
            if x != y
                if is_transition(x, y)
                    S += 1
                else
                    V += 1
                end
            end
        end
    end

    # Convert to proportions
    P = S / total  # proportion of transitions
    Q = V / total  # proportion of transversions

    # Get model type from params
    model_name = get(model.params, :model, nothing)

    if model_name == "JC69"
        # JC69: d = -3/4 * ln(1 - 4/3 * p), where p = (P + Q)
        p = (P + Q)
        return p >= 0.75 ? Inf : -3/4 * log(1 - 4/3 * p)

    elseif model_name == "K2P"  # Kimura 2-Parameter
        # K2P: d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
        if 2P + Q >= 1 || 2Q >= 1
            return Inf
        end
        return -1/2 * log(1 - 2P - Q) - 1/4 * log(1 - 2Q)

    elseif model_name == "HKY85"
        # For HKY85, use the modified formula taking into account base frequencies
        π = model.π  # stationary frequencies
        πY = π[2] + π[4]  # πC + πT
        πR = π[1] + π[3]  # πA + πG

        # Modified formula based on Felsenstein (2004)
        if 2P + Q >= 1 || 2Q >= 1
            return Inf
        end

        k = get(model.params, :kappa, 1.0)  # transition/transversion rate ratio

        return -2 * (πR * πY * k + πR * πY) * log(1 - P/(2πR * πY) - Q/(2πR * πY)) -
               2 * (πR * πY) * log(1 - Q/(2πR * πY))

    elseif model_name == "GTR"
        # For GTR, fall back to ML estimation as there's no simple analytical formula
        return NaN

    else
        error("Unsupported model for analytical distance: $model_name")
    end
end

"""
    compute_distances(model::Model, aln::LabeledAlignment;
                     method::Symbol=:analytical,
                     initial_scale::Float64=0.1,
                     max_scale::Float64=100.0) -> NamedTuple

Compute pairwise distances between sequences.

# Arguments
- `model::Model`: Evolution model (DNA or protein)
- `aln::LabeledAlignment`: Aligned sequences with labels
- `method::Symbol=:analytical`: Method to use (:analytical, :ml, or :auto)
- `initial_scale::Float64=0.1`: Initial guess for ML optimization
- `max_scale::Float64=100.0`: Maximum allowed distance

# Returns
NamedTuple containing:
- `distances::Matrix{Float64}`: Distance matrix
- `labels::Vector{String}`: Sequence labels

Note: :auto will use analytical formulas when available, falling back to ML otherwise
"""
function compute_distances(
    model::Model,
    aln::LabeledAlignment;
    method::Symbol=:analytical,
    initial_scale::Float64=0.1,
    max_scale::Float64=100.0
)
    n = length(aln.sequences)
    D = zeros(n, n)

    # Determine if analytical method is available
    model_name = get(model.params, :model, nothing)
    has_analytical = model_name in ("JC69", "K2P", "HKY85")

    # Validate method choice
    if method == :analytical && !has_analytical
        error("Analytical method not available for model: $model_name")
    end

    # Choose method
    use_ml = if method == :auto
        !has_analytical
    elseif method == :analytical
        false
    elseif method == :ml
        true
    else
        error("Unknown method: $method")
    end

    if use_ml
        # ML implementation
        n_params = num_parameters(model)

        if n_params == 1
            opt_options = Optim.Options(
                x_tol = 1e-8,
                f_tol = 1e-8,
                iterations = 1000
            )

            Threads.@threads for i in 1:n
                for j in i+1:n
                    # Define negative log-likelihood function for optimization
                    function f(t)
                        if t < 0 || t > max_scale
                            return Inf
                        end
                        try
                            return -sequence_likelihood(model, aln.sequences[i], aln.sequences[j], first(t))
                        catch
                            return Inf
                        end
                    end

                    # Run optimization
                    res = try
                        optimize(f, [0.0], [max_scale], [initial_scale])
                    catch e
                        @warn "Optimization failed for sequences $(aln.labels[i]) and $(aln.labels[j]): $(sprint(showerror, e))"
                        nothing
                    end

                    # Extract result
                    if res !== nothing && Optim.converged(res)
                        D[i,j] = D[j,i] = Optim.minimizer(res)[1]
                    else
                        D[i,j] = D[j,i] = NaN
                    end
                end
            end
        else
            @warn "Multi-parameter ML optimization not yet implemented for model type: $model_name"
            D .= NaN
        end
    else
        # Analytical formula implementation
        Threads.@threads for i in 1:n
            for j in i+1:n
                D[i,j] = D[j,i] = compute_analytical_distance(model, aln.sequences[i], aln.sequences[j])
            end
        end
    end

    # Add computation method to return value for reference
    return (
        distances=D,
        labels=aln.labels,
        method=use_ml ? :ml : :analytical,
        model=model_name
    )
end

# Helper function to determine number of parameters
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
            str = isnan(val) ? "NA" : "$(val)"
            print(rpad(str, num_pad))
        end
        println()
    end
end

end # module
