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
