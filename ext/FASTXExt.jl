module FASTXExt

using FASTX
using BioSequences
using EvolutionModels
using Optim

# This line is needed for extensions
using EvolutionModels: Model, sequence_likelihood

# Register the extension functions with the main module
import EvolutionModels: read_alignment, compute_distances

"""
    read_alignment(fasta_file::String) -> Vector{LongSequence}

Read multiple sequence alignment from FASTA file.
Validates that all sequences are of equal length.
"""
function read_alignment(fasta_file::String)
    reader = FASTA.Reader(open(fasta_file))
    sequences = [sequence(record) for record in reader]
    close(reader)
    
    # Validate alignment
    length(sequences) > 0 || error("Empty alignment file")
    ref_length = length(sequences[1])
    all(length(seq) == ref_length for seq in sequences) ||
        error("Sequences must be of equal length")
    
    return sequences
end

"""
    compute_distances(model::Model, alignment::Vector{LongSequence}; 
                     optimizer=LBFGS()) -> Matrix{Float64}

Compute pairwise distances between sequences under the given model.
Distances are estimated by maximum likelihood.

# Arguments
- `model::Model`: Evolution model (DNA or protein)
- `alignment::Vector{LongSequence}`: Vector of aligned sequences
- `optimizer=LBFGS()`: Optimization algorithm from Optim.jl

# Returns
- `Matrix{Float64}`: Distance matrix where d[i,j] is the ML estimate
                     of evolution time between sequences i and j

# Example
```julia
# Read alignment and compute distances
seqs = read_alignment("alignment.fasta")
model = create_model(JC69Model, 0.1)
D = compute_distances(model, seqs)
```
"""
function compute_distances(
    model::Model,
    alignment::Vector{LongSequence};
    optimizer=LBFGS()
)
    n = length(alignment)
    D = zeros(n, n)
    
    for i in 1:n, j in i+1:n
        # Define negative log-likelihood function
        f(t) = -sequence_likelihood(model, alignment[i], alignment[j], first(t))
        
        # Optimize time parameter
        res = optimize(f, [0.0], [100.0], [0.1], optimizer)
        D[i,j] = D[j,i] = Optim.minimizer(res)[1]
    end
    
    return D
end

end # module
