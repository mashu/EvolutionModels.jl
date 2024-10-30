"""
Create a JC69 model
"""
function create_model(::Type{JC69Model}, μ::Float64)
    μ > 0 || throw(ArgumentError("μ must be positive"))
    seq_type = DNAType()
    n = length(seq_type)
    
    π = fill(1.0/n, n)
    R = ones(n, n) - I
    R .*= μ
    
    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    params = Dict{Symbol,Any}(:model => "JC69")
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end

"""
Create an HKY85 model
"""
function create_model(::Type{HKY85Model}, μ::Float64, π::Vector{Float64}, κ::Float64)
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
    params = Dict{Symbol,Any}(:model => "HKY85", :kappa => κ)
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end

"""
Create a GTR model
"""
function create_model(::Type{GTRModel}, μ::Float64, π::Vector{Float64}, rates::Matrix{Float64})
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
    params = Dict{Symbol,Any}(:model => "GTR")
    
    Model{DNAType}(seq_type, μ, π, R, P, Q, params)
end
