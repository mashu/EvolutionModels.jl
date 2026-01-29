using LinearAlgebra

"""
Get index of a symbol in the canonical ordering
"""
function symbol_index(symbol::DNA, seq_type::DNAType)
    idx = findfirst(==(symbol), seq_type.symbols)
    isnothing(idx) && throw(ArgumentError("Invalid DNA symbol: $symbol"))
    idx
end

function symbol_index(symbol::AminoAcid, seq_type::ProteinType)
    idx = findfirst(==(symbol), seq_type.symbols)
    isnothing(idx) && throw(ArgumentError("Invalid amino acid symbol: $symbol"))
    idx
end

"""
Check if two symbols are transitions (A↔G or C↔T)
"""
function is_transition(x::DNA, y::DNA)
    (x == DNA_A && y == DNA_G) || (x == DNA_G && y == DNA_A) ||
    (x == DNA_C && y == DNA_T) || (x == DNA_T && y == DNA_C)
end

"""
Validate frequencies vector
"""
function validate_frequencies(π::Vector{Float64}, seq_type::S) where S<:SequenceType
    length(π) == length(seq_type) || 
        throw(ArgumentError("π must have length $(length(seq_type))"))
    isapprox(sum(π), 1.0, rtol=1e-10) || 
        throw(ArgumentError("π must sum to 1"))
    all(x -> x ≥ 0, π) || 
        throw(ArgumentError("π must be non-negative"))
end

"""
Compute Q matrix (generator matrix) from R and P.
"""
function compute_q_matrix(R::Matrix{Float64}, P::Diagonal{Float64})
    Q = R * P
    for i in 1:size(Q,1)
        Q[i,i] = -sum(Q[i,:])
    end
    Q
end

"""
    normalize_rate_matrix!(Q::Matrix{Float64}, π::Vector{Float64})

Normalize the rate matrix Q in-place so that the expected number of 
substitutions per unit time equals 1.

The normalization factor is β = -∑ᵢ πᵢ Qᵢᵢ, and Q is divided by β.

# Arguments
- `Q::Matrix{Float64}`: Rate matrix to normalize (modified in-place)
- `π::Vector{Float64}`: Stationary frequencies

# Returns
- `Float64`: The normalization factor β used
"""
function normalize_rate_matrix!(Q::Matrix{Float64}, π::Vector{Float64})
    β = -sum(π[i] * Q[i,i] for i in 1:length(π))
    if β > 0
        Q ./= β
    end
    return β
end

"""
    expected_substitution_rate(Q::Matrix{Float64}, π::Vector{Float64}) -> Float64

Compute the expected number of substitutions per unit time for a rate matrix Q
with stationary distribution π.

# Formula
```math
\\beta = -\\sum_i \\pi_i Q_{ii}
```

# Arguments
- `Q::Matrix{Float64}`: Rate matrix
- `π::Vector{Float64}`: Stationary frequencies

# Returns
- `Float64`: Expected substitutions per unit time
"""
function expected_substitution_rate(Q::Matrix{Float64}, π::Vector{Float64})
    return -sum(π[i] * Q[i,i] for i in 1:length(π))
end
