using BioSequences

# Canonical ordering of standard nucleotides and amino acids
const STANDARD_DNA = (DNA_A, DNA_C, DNA_G, DNA_T)
const STANDARD_AA = (
    AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I,
    AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V
)

export STANDARD_DNA, STANDARD_AA

"""
    SequenceType

Abstract type for sequence alphabets (DNA or Protein).
"""
abstract type SequenceType end

struct DNAType <: SequenceType 
    symbols::NTuple{4,DNA}
    DNAType() = new(STANDARD_DNA)
end

struct ProteinType <: SequenceType
    symbols::NTuple{20,AminoAcid}
    ProteinType() = new(STANDARD_AA)
end

Base.length(::DNAType) = 4
Base.length(::ProteinType) = 20
symbols(t::SequenceType) = t.symbols

"""
    EvolutionModel

Abstract type for all molecular evolution models.
"""
abstract type EvolutionModel end

"""
    DNAModel <: EvolutionModel

Abstract type for DNA evolution models.
"""
abstract type DNAModel <: EvolutionModel end

"""
    ProteinModel <: EvolutionModel

Abstract type for protein evolution models.
"""
abstract type ProteinModel <: EvolutionModel end

# Concrete model types
struct JC69Model <: DNAModel end
struct HKY85Model <: DNAModel end
struct GTRModel <: DNAModel end
struct WAGModel <: ProteinModel end
struct LGModel <: ProteinModel end

"""
    Model{S<:SequenceType}

Concrete implementation of an evolution model.
"""
struct Model{S<:SequenceType}
    seq_type::S
    μ::Float64                    # Overall rate scaling
    π::Vector{Float64}            # Stationary frequencies
    R::Matrix{Float64}            # Rate matrix (scaled by μ)
    P::Matrix{Float64}            # Diagonal frequency matrix
    Q::Matrix{Float64}            # Generator matrix
    params::Dict{Symbol,Any}      # Model-specific parameters
end
