"""
    create_model(::Type{JC69Model}, μ::Float64; normalize::Bool=false) -> Model{DNAType}

Create a Jukes-Cantor (1969) model of DNA evolution.

The JC69 model assumes:
- Equal base frequencies (π = 0.25 for all bases)
- Equal substitution rates between all nucleotides
- Single parameter μ controlling overall mutation rate

# Arguments
- `μ::Float64`: Overall mutation rate scaling factor (must be positive)
- `normalize::Bool=false`: If true, normalize Q so expected substitutions per unit time = 1

# Returns
- `Model{DNAType}`: A JC69 model with rate matrix Q where:
  * qᵢⱼ = μ/4 for i ≠ j (off-diagonal elements)  
  * qᵢᵢ = -3μ/4 (diagonal elements)
  * If normalized: qᵢⱼ = 1/3 for i ≠ j, qᵢᵢ = -1

# Example
```julia
# Create JC69 model with rate 0.1
model = create_model(JC69Model, 0.1)

# Create normalized model (1 expected substitution per unit time)
model_norm = create_model(JC69Model, 1.0, normalize=true)

# Access model properties
π = stationary_frequencies(model)  # All 0.25
R = rate_matrix(model)  # Symmetric rate matrix
```

# References
Jukes, T.H. and Cantor, C.R. (1969). Evolution of Protein Molecules.
"""
function create_model(::Type{JC69Model}, μ::Float64; normalize::Bool=false)
    μ > 0 || throw(ArgumentError("μ must be positive"))
    seq_type = DNAType()
    n = length(seq_type)
    
    π = fill(1.0/n, n)
    R = ones(n, n) - I
    R .*= μ
    
    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    
    if normalize
        normalize_rate_matrix!(Q, π)
    end
    
    params = Dict{Symbol,Any}(:model => "JC69", :normalized => normalize)
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end

"""
    create_model(::Type{HKY85Model}, μ::Float64, π::Vector{Float64}, κ::Float64; normalize::Bool=false) -> Model{DNAType}

Create a Hasegawa-Kishino-Yano (1985) model of DNA evolution.

The HKY85 model features:
- Arbitrary base frequencies π
- Different rates for transitions (A↔G, C↔T) vs transversions
- Rate ratio κ between transitions and transversions

# Arguments
- `μ::Float64`: Overall mutation rate scaling factor (must be positive)
- `π::Vector{Float64}`: Vector of base frequencies [πA, πC, πG, πT], must sum to 1
- `κ::Float64`: Transition/transversion rate ratio (must be positive)
- `normalize::Bool=false`: If true, normalize Q so expected substitutions per unit time = 1

# Returns
- `Model{DNAType}`: An HKY85 model with rate matrix Q where:
  * qᵢⱼ = μκπⱼ for transitions
  * qᵢⱼ = μπⱼ for transversions
  * qᵢᵢ = -∑ⱼ≠ᵢ qᵢⱼ

# Example
```julia
# Create HKY85 model with unequal base frequencies
π = [0.3, 0.2, 0.2, 0.3]  # A,C,G,T frequencies
model = create_model(HKY85Model, 0.1, π, 2.0)  # κ=2 means transitions occur 2× faster

# Create normalized model
model_norm = create_model(HKY85Model, 1.0, π, 2.0, normalize=true)
```

# References
Hasegawa, M., Kishino, H., and Yano, T. (1985). Dating of the human-ape splitting
by a molecular clock of mitochondrial DNA. J. Mol. Evol., 22(2):160-174.
"""
function create_model(::Type{HKY85Model}, μ::Float64, π::Vector{Float64}, κ::Float64; normalize::Bool=false)
    μ > 0 || throw(ArgumentError("μ must be positive"))
    κ > 0 || throw(ArgumentError("κ must be positive"))
    
    seq_type = DNAType()
    validate_frequencies(π, seq_type)
    
    n = length(seq_type)
    R = zeros(n, n)
    
    for i in 1:n, j in i+1:n
        if is_transition(seq_type.symbols[i], seq_type.symbols[j])
            R[i,j] = R[j,i] = κ
        else
            R[i,j] = R[j,i] = 1.0
        end
    end
    
    R .*= μ
    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    
    if normalize
        normalize_rate_matrix!(Q, π)
    end
    
    params = Dict{Symbol,Any}(:model => "HKY85", :kappa => κ, :normalized => normalize)
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end

"""
    create_model(::Type{GTRModel}, μ::Float64, π::Vector{Float64}, rates::Matrix{Float64}; normalize::Bool=false) -> Model{DNAType}

Create a General Time-Reversible (GTR) model of DNA evolution.

The GTR model is the most general neutral, independent sites, reversible model:
- Arbitrary base frequencies π
- Arbitrary symmetric rate matrix
- Satisfies detailed balance: πᵢqᵢⱼ = πⱼqⱼᵢ

# Arguments
- `μ::Float64`: Overall mutation rate scaling factor (must be positive)
- `π::Vector{Float64}`: Vector of base frequencies [πA, πC, πG, πT], must sum to 1
- `rates::Matrix{Float64}`: 4×4 symmetric matrix of relative substitution rates
- `normalize::Bool=false`: If true, normalize Q so expected substitutions per unit time = 1

# Returns
- `Model{DNAType}`: A GTR model with rate matrix Q = μR⋅diag(π) where:
  * R is the symmetric rate matrix
  * qᵢⱼ = μrᵢⱼπⱼ for i ≠ j
  * qᵢᵢ = -∑ⱼ≠ᵢ qᵢⱼ

# Example
```julia
# Create GTR model with custom rates
π = [0.3, 0.2, 0.2, 0.3]  # Base frequencies
rates = [0.0 1.0 2.0 1.0;  # Symmetric rate matrix
         1.0 0.0 1.0 2.0;
         2.0 1.0 0.0 1.0;
         1.0 2.0 1.0 0.0]
model = create_model(GTRModel, 0.1, π, rates)

# Create normalized model
model_norm = create_model(GTRModel, 1.0, π, rates, normalize=true)
```

# References
Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis
of DNA sequences. Lectures on Mathematics in the Life Sciences, 17:57-86.
"""
function create_model(::Type{GTRModel}, μ::Float64, π::Vector{Float64}, rates::Matrix{Float64}; normalize::Bool=false)
    μ > 0 || throw(ArgumentError("μ must be positive"))
    seq_type = DNAType()
    validate_frequencies(π, seq_type)
    
    n = length(seq_type)
    size(rates) == (n,n) || throw(DimensionMismatch("rates matrix must be $(n)×$(n)"))
    issymmetric(rates) || throw(ArgumentError("rates matrix must be symmetric"))
    
    R = copy(rates)
    R .*= μ
    
    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    
    if normalize
        normalize_rate_matrix!(Q, π)
    end
    
    params = Dict{Symbol,Any}(:model => "GTR", :normalized => normalize)
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end
