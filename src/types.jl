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

"""
    GammaRateModel{S<:SequenceType}

Evolution model with discrete Gamma rate variation across sites.

Wraps a base Model and adds rate categories from a discretized Gamma distribution.
This accounts for site-to-site rate heterogeneity common in biological sequences.

# Fields
- `base_model::Model{S}`: The underlying substitution model
- `α::Float64`: Gamma shape parameter (smaller = more variation)
- `rates::Vector{Float64}`: Discrete rate categories
- `weights::Vector{Float64}`: Probability of each rate category (typically uniform)
"""
struct GammaRateModel{S<:SequenceType}
    base_model::Model{S}
    α::Float64
    rates::Vector{Float64}
    weights::Vector{Float64}
end

"""
    PartitionModel

Model that applies different evolutionary models to different sequence regions.

Useful for analyzing sequences with distinct evolutionary patterns, such as
antibody CDR vs framework regions.

# Fields
- `partitions::Vector{Tuple{UnitRange{Int}, Model}}`: List of (site_range, model) pairs
- `total_length::Int`: Total sequence length covered
"""
struct PartitionModel
    partitions::Vector{Tuple{UnitRange{Int}, Model}}
    total_length::Int
    
    function PartitionModel(partitions::Vector{<:Tuple{UnitRange{Int}, Model}})
        # Validate partitions don't overlap and cover sequence
        sorted = sort(partitions, by=p->first(p[1]))
        total_length = 0
        for (i, (range, _)) in enumerate(sorted)
            if i > 1
                prev_end = last(sorted[i-1][1])
                if first(range) <= prev_end
                    throw(ArgumentError("Partition ranges overlap"))
                end
            end
            total_length = max(total_length, last(range))
        end
        new(collect(partitions), total_length)
    end
end
