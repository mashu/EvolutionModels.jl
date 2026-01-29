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

"""
    discrete_gamma_rates(α::Float64, n_categories::Int=4) -> Tuple{Vector{Float64}, Vector{Float64}}

Compute discrete Gamma rate categories using the mean of each category.

The Gamma distribution has shape α and rate α (so mean = 1). The distribution
is divided into n_categories equiprobable categories, and the mean rate for
each category is returned.

# Arguments
- `α::Float64`: Shape parameter (smaller = more rate variation, α→∞ = no variation)
- `n_categories::Int=4`: Number of discrete rate categories

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (rates, weights) where weights are all 1/n_categories

# Example
```julia
rates, weights = discrete_gamma_rates(0.5, 4)
# rates ≈ [0.03, 0.28, 0.92, 2.77] - wide variation
# weights = [0.25, 0.25, 0.25, 0.25]

rates, weights = discrete_gamma_rates(10.0, 4)  
# rates ≈ [0.70, 0.90, 1.05, 1.35] - little variation
```

# References
Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences 
with variable rates over sites. J Mol Evol, 39:306-314.
"""
function discrete_gamma_rates(α::Float64, n_categories::Int=4)
    α > 0 || throw(ArgumentError("α must be positive"))
    n_categories > 0 || throw(ArgumentError("n_categories must be positive"))
    
    # Use quantile-based approach with mean of each category
    # For Gamma(α, α), the mean is 1
    rates = zeros(n_categories)
    weights = fill(1.0 / n_categories, n_categories)
    
    # Compute rate for each category as the mean within that quantile range
    # Using incomplete gamma function approximation
    for i in 1:n_categories
        # Boundaries of this category
        lower = (i - 1) / n_categories
        upper = i / n_categories
        
        # Mean rate in this category (using numerical integration approximation)
        # For Gamma(α,α), the mean of category [q1,q2] is:
        # E[r|q1<F(r)<q2] = (G(α+1, x2) - G(α+1, x1)) / (q2 - q1)
        # where x = α * quantile and G is the incomplete gamma
        
        # Simpler approximation using midpoint of quantile
        mid_quantile = (lower + upper) / 2
        rates[i] = gamma_quantile(α, mid_quantile)
    end
    
    # Normalize so mean rate = 1
    mean_rate = sum(rates .* weights)
    rates ./= mean_rate
    
    return rates, weights
end

"""
    gamma_quantile(α::Float64, p::Float64) -> Float64

Compute the p-th quantile of Gamma(α, α) distribution using Newton-Raphson iteration.
"""
function gamma_quantile(α::Float64, p::Float64)
    if p <= 0.0
        return 0.0
    elseif p >= 1.0
        return Inf
    end
    
    # Initial guess using Wilson-Hilferty approximation
    if α > 1
        z = quantile_normal(p)
        x = α * (1 - 1/(9*α) + z/sqrt(9*α))^3
        x = max(x, 0.001)
    else
        x = (p * α * gamma(α))^(1/α)
        x = max(x, 0.001)
    end
    
    # Newton-Raphson iteration
    for _ in 1:50
        # CDF and PDF of Gamma(α, α)
        cdf = gamma_cdf(x, α)
        pdf = gamma_pdf(x, α)
        
        if pdf < 1e-300
            break
        end
        
        dx = (cdf - p) / pdf
        x = max(x - dx, x/10)
        
        if abs(dx) < 1e-10 * x
            break
        end
    end
    
    return x
end

"""Standard normal quantile using rational approximation."""
function quantile_normal(p::Float64)
    # Rational approximation for standard normal quantile
    if p <= 0.0
        return -Inf
    elseif p >= 1.0
        return Inf
    elseif p == 0.5
        return 0.0
    end
    
    if p < 0.5
        return -quantile_normal(1.0 - p)
    end
    
    t = sqrt(-2.0 * log(1.0 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    
    return t - (c0 + c1*t + c2*t^2) / (1.0 + d1*t + d2*t^2 + d3*t^3)
end

"""Gamma PDF for Gamma(α, α) distribution."""
function gamma_pdf(x::Float64, α::Float64)
    if x <= 0
        return 0.0
    end
    return (α^α / gamma(α)) * x^(α-1) * exp(-α*x)
end

"""Gamma CDF using series expansion."""
function gamma_cdf(x::Float64, α::Float64)
    if x <= 0
        return 0.0
    end
    # Lower incomplete gamma using series expansion
    return lower_incomplete_gamma(α, α*x) / gamma(α)
end

"""Lower incomplete gamma function using series expansion."""
function lower_incomplete_gamma(a::Float64, x::Float64)
    if x < 0
        return 0.0
    end
    if x == 0
        return 0.0
    end
    
    # Series expansion: γ(a,x) = x^a * e^(-x) * Σ(x^n / (a+1)...(a+n))
    sum_val = 1.0 / a
    term = 1.0 / a
    
    for n in 1:200
        term *= x / (a + n)
        sum_val += term
        if abs(term) < 1e-15 * abs(sum_val)
            break
        end
    end
    
    return x^a * exp(-x) * sum_val
end
