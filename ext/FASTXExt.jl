module FASTXExt

using FASTX
using BioSequences
using Optim

import EvolutionModels: Model, sequence_likelihood, read_alignment, compute_distances
import EvolutionModels: JC69Model, HKY85Model, GTRModel

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
    compute_distances(model::Model, aln::LabeledAlignment;
                     initial_scale::Float64=0.1,
                     max_scale::Float64=100.0) -> NamedTuple

Compute pairwise distances between sequences under the given model.
"""
function compute_distances(
    model::Model,
    aln::LabeledAlignment;
    initial_scale::Float64=0.1,
    max_scale::Float64=100.0
)
    n = length(aln.sequences)
    D = zeros(n, n)

    # Get number of parameters for this model
    n_params = num_parameters(model)

    # For single-parameter models
    if n_params == 1
        Threads.@threads for i in 1:n
            for j in i+1:n
                # Define negative log-likelihood function for optimization
                function f(t)
                    if t < 0 || t > max_scale
                        return Inf
                    end
                    try
                        return -sequence_likelihood(model, aln.sequences[i], aln.sequences[j], t)
                    catch
                        return Inf
                    end
                end

                # Run optimization using Brent's method
                res = try
                    optimize(t -> f(first(t)), [0.0], [max_scale], [initial_scale])
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
        # For multi-parameter models
        @warn "Multi-parameter optimization not yet implemented for model type: $(model.params[:model])"
        D .= NaN
    end

    return (distances=D, labels=aln.labels)
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
