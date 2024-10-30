using Random
using StatsBase

"""
Simulate sequence evolution under the model
"""
function evolve_sequence(
    model::Model{DNAType},
    seq::LongDNA{4},
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    P = transition_probability_matrix(model, t)
    result = LongDNA{4}(length(seq))
    
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
Simulate sequence evolution under the model
"""
function evolve_sequence(
    model::Model{ProteinType},
    seq::LongAA,
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    P = transition_probability_matrix(model, t)
    # Create an empty amino acid sequence with same length
    result = LongAA([AA_A for _ in 1:length(seq)])  # Initialize with temporary values
    
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
Compute the likelihood of observing sequence2 after time t given sequence1
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
