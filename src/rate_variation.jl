"""
Rate variation models: Gamma rate heterogeneity and partition models.
"""

using Random
using StatsBase

"""
    create_gamma_model(base_model::Model, α::Float64; n_categories::Int=4)

Create a model with discrete Gamma rate variation across sites.

# Arguments
- `base_model::Model`: The underlying substitution model (e.g., LG, WAG, HKY85)
- `α::Float64`: Gamma shape parameter. Smaller values = more rate variation:
  - α < 1: High variation (some sites very slow, some very fast)
  - α ≈ 1: Moderate variation
  - α > 2: Low variation
  - α → ∞: No variation (equivalent to base model)
- `n_categories::Int=4`: Number of discrete rate categories (typically 4-8)

# Returns
- `GammaRateModel`: Model with rate heterogeneity

# Example
```julia
# Create LG+G4 model (LG with 4 Gamma rate categories)
base = create_model(LGModel, 1.0, normalize=true)
model = create_gamma_model(base, 0.5)  # α=0.5 for high rate variation

# For antibodies, moderate variation is common
model = create_gamma_model(base, 1.0)  # α=1.0
```

# References
Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences
with variable rates over sites. J Mol Evol, 39:306-314.
"""
function create_gamma_model(base_model::Model{S}, α::Float64; n_categories::Int=4) where S
    α > 0 || throw(ArgumentError("α must be positive"))
    n_categories > 0 || throw(ArgumentError("n_categories must be positive"))
    
    rates, weights = discrete_gamma_rates(α, n_categories)
    
    return GammaRateModel{S}(base_model, α, rates, weights)
end

"""
    create_partition_model(partitions::Pair{UnitRange{Int}, Model}...)

Create a partition model with different models for different sequence regions.

# Arguments
- `partitions...`: Variable number of `range => model` pairs

# Returns
- `PartitionModel`: Model with region-specific parameters

# Example
```julia
# Antibody with different models for CDR and framework
framework_model = create_model(LGModel, 1.0, normalize=true)
cdr_model = create_model(LGModel, 2.0, normalize=true)  # Higher rate for CDRs

# Define regions (IMGT numbering example)
model = create_partition_model(
    1:26 => framework_model,     # FR1
    27:38 => cdr_model,          # CDR1
    39:55 => framework_model,    # FR2
    56:65 => cdr_model,          # CDR2
    66:104 => framework_model,   # FR3
    105:117 => cdr_model,        # CDR3
    118:128 => framework_model   # FR4
)
```
"""
function create_partition_model(partitions::Pair{UnitRange{Int}, <:Model}...)
    partition_vec = [(range, model) for (range, model) in partitions]
    return PartitionModel(partition_vec)
end

# Accessor functions for GammaRateModel
stationary_frequencies(model::GammaRateModel) = model.base_model.π
rate_matrix(model::GammaRateModel) = model.base_model.R

"""
    transition_probability_matrix(model::GammaRateModel, t::Float64, rate_category::Int)

Compute P(t) for a specific rate category.
"""
function transition_probability_matrix(model::GammaRateModel, t::Float64, rate_category::Int)
    r = model.rates[rate_category]
    return exp(model.base_model.Q * t * r)
end

"""
    sequence_likelihood(model::GammaRateModel, seq1, seq2, t::Float64; joint::Bool=false)

Compute likelihood integrating over Gamma rate categories.

For each site, the likelihood is averaged over all rate categories:
```math
P(y_k | x_k, t) = \\sum_r w_r P(t \\cdot r)_{x_k, y_k}
```
"""
function sequence_likelihood(
    model::GammaRateModel,
    seq1::Union{LongDNA{4},LongAA},
    seq2::Union{LongDNA{4},LongAA},
    t::Float64;
    joint::Bool=false
)
    length(seq1) == length(seq2) || 
        throw(ArgumentError("Sequences must be same length"))
    
    base = model.base_model
    π = base.π
    
    # Precompute P(t) for each rate category
    P_matrices = [transition_probability_matrix(model, t, r) for r in 1:length(model.rates)]
    
    logL = 0.0
    
    for (x, y) in zip(seq1, seq2)
        if (x in symbols(base.seq_type)) && (y in symbols(base.seq_type))
            i = symbol_index(x, base.seq_type)
            j = symbol_index(y, base.seq_type)
            
            # Average over rate categories
            site_prob = sum(model.weights[r] * P_matrices[r][i,j] for r in 1:length(model.rates))
            
            if joint
                logL += log(π[i]) + log(site_prob)
            else
                logL += log(site_prob)
            end
        end
    end
    
    return logL
end

"""
    sequence_likelihood(model::PartitionModel, seq1, seq2, t::Float64; joint::Bool=false)

Compute likelihood using region-specific models.
"""
function sequence_likelihood(
    model::PartitionModel,
    seq1::Union{LongDNA{4},LongAA},
    seq2::Union{LongDNA{4},LongAA},
    t::Float64;
    joint::Bool=false
)
    length(seq1) == length(seq2) || 
        throw(ArgumentError("Sequences must be same length"))
    
    logL = 0.0
    
    for (range, region_model) in model.partitions
        # Extract subsequences for this region
        if last(range) <= length(seq1)
            sub1 = seq1[range]
            sub2 = seq2[range]
            logL += sequence_likelihood(region_model, sub1, sub2, t, joint=joint)
        end
    end
    
    return logL
end

"""
    evolve_sequence(model::GammaRateModel{DNAType}, seq::LongDNA{4}, t::Float64; rng=Random.GLOBAL_RNG)

Simulate evolution with site-specific rates drawn from Gamma distribution.
"""
function evolve_sequence(
    model::GammaRateModel{DNAType},
    seq::LongDNA{4},
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    base = model.base_model
    n_sites = length(seq)
    n_cats = length(model.rates)
    
    # Assign each site to a rate category
    site_categories = sample(rng, 1:n_cats, Weights(model.weights), n_sites)
    
    # Initialize result
    result = LongDNA{4}([seq[1] for _ in 1:n_sites])
    
    for (i, nt) in enumerate(seq)
        if nt in STANDARD_DNA
            r = model.rates[site_categories[i]]
            P = exp(base.Q * t * r)
            
            current_idx = symbol_index(nt, base.seq_type)
            new_idx = sample(rng, 1:4, Weights(P[current_idx, :]))
            result[i] = base.seq_type.symbols[new_idx]
        else
            result[i] = nt
        end
    end
    
    return result
end

"""
    evolve_sequence(model::GammaRateModel{ProteinType}, seq::LongAA, t::Float64; rng=Random.GLOBAL_RNG)

Simulate protein evolution with site-specific rates.
"""
function evolve_sequence(
    model::GammaRateModel{ProteinType},
    seq::LongAA,
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    base = model.base_model
    n_sites = length(seq)
    n_cats = length(model.rates)
    
    # Assign each site to a rate category
    site_categories = sample(rng, 1:n_cats, Weights(model.weights), n_sites)
    
    # Initialize result
    result = LongAA([seq[1] for _ in 1:n_sites])
    
    for (i, aa) in enumerate(seq)
        if aa in STANDARD_AA
            r = model.rates[site_categories[i]]
            P = exp(base.Q * t * r)
            
            current_idx = symbol_index(aa, base.seq_type)
            new_idx = sample(rng, 1:20, Weights(P[current_idx, :]))
            result[i] = base.seq_type.symbols[new_idx]
        else
            result[i] = aa
        end
    end
    
    return result
end

"""
    evolve_sequence(model::PartitionModel, seq::Union{LongDNA{4},LongAA}, t::Float64; rng=Random.GLOBAL_RNG)

Simulate evolution with region-specific models.
"""
function evolve_sequence(
    model::PartitionModel,
    seq::LongDNA{4},
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    result = LongDNA{4}([seq[i] for i in 1:length(seq)])
    
    for (range, region_model) in model.partitions
        if last(range) <= length(seq)
            sub_seq = seq[range]
            evolved_sub = evolve_sequence(region_model, sub_seq, t, rng=rng)
            for (i, pos) in enumerate(range)
                result[pos] = evolved_sub[i]
            end
        end
    end
    
    return result
end

function evolve_sequence(
    model::PartitionModel,
    seq::LongAA,
    t::Float64;
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    result = LongAA([seq[i] for i in 1:length(seq)])
    
    for (range, region_model) in model.partitions
        if last(range) <= length(seq)
            sub_seq = seq[range]
            evolved_sub = evolve_sequence(region_model, sub_seq, t, rng=rng)
            for (i, pos) in enumerate(range)
                result[pos] = evolved_sub[i]
            end
        end
    end
    
    return result
end
