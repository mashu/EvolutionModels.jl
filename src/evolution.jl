using Random
using StatsBase

"""
    evolve_sequence(model::Model{DNAType}, seq::LongDNA{4}, t::Float64; 
                   rng::AbstractRNG=Random.GLOBAL_RNG) -> LongDNA{4}

Simulate DNA sequence evolution under a continuous-time Markov model for time t.

The function implements the following process:
1. Computes transition probability matrix P(t) = exp(Qt)
2. For each site, samples new nucleotide from the probability distribution
   given by the corresponding row of P(t)
3. Preserves non-standard nucleotides (gaps, ambiguous bases) unchanged

# Arguments
- `model::Model{DNAType}`: DNA evolution model (JC69, HKY85, or GTR)
- `seq::LongDNA{4}`: Input DNA sequence
- `t::Float64`: Evolution time (branch length)
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator

# Returns
- `LongDNA{4}`: Evolved DNA sequence of the same length as input

# Example
```julia
# Create model and evolve sequence
model = create_model(JC69Model, 0.1)
seq = dna"ATCG"
evolved = evolve_sequence(model, seq, 1.0)

# Use specific RNG for reproducibility
rng = MersenneTwister(42)
evolved = evolve_sequence(model, seq, 1.0, rng=rng)
```

Note: For very small t, few changes are expected. As t→∞, nucleotide frequencies
approach the model's stationary distribution π.
"""
function evolve_sequence(
    model::Model{DNAType},
    seq::LongDNA{4},
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    P = transition_probability_matrix(model, t)
    # Initialize with first sequence character
    result = LongDNA{4}([seq[1] for _ in 1:length(seq)])
    
    for (i, nt) in enumerate(seq)
        if nt in STANDARD_DNA
            current_idx = symbol_index(nt, model.seq_type)
            new_idx = sample(rng, 1:4, Weights(P[current_idx, :]))
            result[i] = model.seq_type.symbols[new_idx]
        else
            result[i] = nt
        end
    end
    
    return result
end

"""
    evolve_sequence(model::Model{ProteinType}, seq::LongAA, t::Float64; 
                   rng::AbstractRNG=Random.GLOBAL_RNG) -> LongAA

Simulate protein sequence evolution under empirical amino acid substitution models.

The function implements the following process:
1. Computes transition probability matrix P(t) = exp(Qt)
2. For each site, samples new amino acid from the probability distribution
   given by the corresponding row of P(t)
3. Preserves non-standard amino acids (gaps, ambiguous residues) unchanged

# Arguments
- `model::Model{ProteinType}`: Protein evolution model (WAG or LG)
- `seq::LongAA`: Input amino acid sequence
- `t::Float64`: Evolution time (branch length)
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator

# Returns
- `LongAA`: Evolved amino acid sequence of the same length as input

# Example
```julia
# Create model and evolve sequence
model = create_model(WAGModel, 0.1)
seq = aa"ARND"
evolved = evolve_sequence(model, seq, 1.0)
```

Note: The empirical models (WAG, LG) use pre-computed rate matrices derived
from large protein alignments, representing "average" evolutionary patterns.
"""
function evolve_sequence(
    model::Model{ProteinType},
    seq::LongAA,
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    P = transition_probability_matrix(model, t)
    # Initialize with first sequence character
    result = LongAA([seq[1] for _ in 1:length(seq)])
    
    for (i, aa) in enumerate(seq)
        if aa in STANDARD_AA
            current_idx = symbol_index(aa, model.seq_type)
            new_idx = sample(rng, 1:20, Weights(P[current_idx, :]))
            result[i] = model.seq_type.symbols[new_idx]
        else
            result[i] = aa
        end
    end
    
    return result
end

"""
    sequence_likelihood(model::Model, seq1::Union{LongDNA{4},LongAA}, 
                       seq2::Union{LongDNA{4},LongAA}, t::Float64) -> Float64

Compute the log-likelihood of evolving from sequence1 to sequence2 in time t
under the given evolutionary model.

The likelihood is computed as:
```math
L = ∏ᵢ P(t)[x₁ᵢ, x₂ᵢ]
```
where P(t) is the transition probability matrix and x₁ᵢ, x₂ᵢ are the states
at position i in sequences 1 and 2 respectively.

# Arguments
- `model::Model`: Evolution model (DNA or protein)
- `seq1::Union{LongDNA{4},LongAA}`: First sequence
- `seq2::Union{LongDNA{4},LongAA}`: Second sequence
- `t::Float64`: Evolution time between sequences

# Returns
- `Float64`: Log-likelihood of the evolution from seq1 to seq2

# Example
```julia
# Compare likelihoods for different evolution times
model = create_model(JC69Model, 0.1)
seq1 = dna"ATCG"
seq2 = dna"ATTG"

L1 = sequence_likelihood(model, seq1, seq2, 0.1)  # Short time
L2 = sequence_likelihood(model, seq1, seq2, 10.0)  # Long time
```

Note: Only standard nucleotides/amino acids contribute to the likelihood.
Non-standard characters (gaps, ambiguous bases) are ignored in the calculation.
"""
function sequence_likelihood(
    model::Model,
    seq1::Union{LongDNA{4},LongAA},
    seq2::Union{LongDNA{4},LongAA},
    t::Float64
)
    length(seq1) == length(seq2) || 
        throw(ArgumentError("Sequences must be same length"))
        
    P = transition_probability_matrix(model, t)
    logL = 0.0
    
    for (x, y) in zip(seq1, seq2)
        if (x in symbols(model.seq_type)) && (y in symbols(model.seq_type))
            i = symbol_index(x, model.seq_type)
            j = symbol_index(y, model.seq_type)
            logL += log(P[i,j])
        end
    end
    
    return logL
end
