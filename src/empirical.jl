"""
    create_model(::Type{WAGModel}, μ::Float64)

Create a WAG (Whelan And Goldman) protein evolution model with specified rate scaling factor μ.
Returns a `Model{ProteinType}` with empirical substitution rates and frequencies.

The WAG model is derived from 182 protein families using maximum likelihood methods,
providing both substitution rates and equilibrium frequencies estimated from the data.
The amino acid order follows STANDARD_AA: A R N D C Q E G H I L K M F P S T W Y V

# Citation
Whelan, S., & Goldman, N. (2001). A general empirical model of protein evolution
derived from multiple protein families using a maximum-likelihood approach.
Molecular Biology and Evolution, 18(5), 691-699.
DOI: 10.1093/oxfordjournals.molbev.a003851

# Arguments
- `μ::Float64`: Rate scaling factor (must be positive)

# Returns
- `Model{ProteinType}`: Configured WAG model
"""
function create_model(::Type{WAGModel}, μ::Float64)
    μ > 0 || throw(ArgumentError("μ must be positive"))

    seq_type = ProteinType()
    π = load_wag_frequencies()
    R = load_wag_rates() .* μ

    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    params = Dict{Symbol,Any}(:model => "WAG")

    Model{ProteinType}(seq_type, μ, π, R, P, Q, params)
end

"""
    create_model(::Type{LGModel}, μ::Float64)

Create an LG (Le-Gascuel) protein evolution model with specified rate scaling factor μ.
Returns a `Model{ProteinType}` with empirical substitution rates and frequencies.

The LG model is derived from 3,912 protein families using maximum likelihood methods,
improving upon previous models by incorporating rate variation across sites.
The amino acid order follows STANDARD_AA: A R N D C Q E G H I L K M F P S T W Y V

# Citation
Le, S. Q., & Gascuel, O. (2008). An improved general amino acid replacement matrix.
Molecular Biology and Evolution, 25(7), 1307-1320.
DOI: 10.1093/molbev/msn067

# Arguments
- `μ::Float64`: Rate scaling factor (must be positive)

# Returns
- `Model{ProteinType}`: Configured LG model
"""
function create_model(::Type{LGModel}, μ::Float64)
    μ > 0 || throw(ArgumentError("μ must be positive"))

    seq_type = ProteinType()
    π = load_lg_frequencies()
    R = load_lg_rates() .* μ

    P = Diagonal(π)
    Q = compute_q_matrix(R, P)
    params = Dict{Symbol,Any}(:model => "LG")

    Model{ProteinType}(seq_type, μ, π, R, P, Q, params)
end

"""
    load_empirical_data(filepath::String)

Load empirical substitution matrix and frequencies from a data file.
Internal helper function for specific model loaders.

# Format
File should contain:
- Lines 1-19: Lower triangular elements of the 20×20 substitution rate matrix
- Line 22: Vector of 20 equilibrium frequencies
The amino acid order follows STANDARD_AA: A R N D C Q E G H I L K M F P S T W Y V

# Arguments
- `filepath::String`: Path to the data file

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}}`: Rate matrix and frequency vector
"""
function load_empirical_data(filepath::String)
    S_data = zeros(Float64, 20, 20)
    π = zeros(Float64, 20)

    open(filepath) do f
        row = 2  # Start from second row since first row is empty
        for line in eachline(f)
            if row < 21
                rate = parse.(Float64, split(line))

                # Fill the lower triangular part for the current row
                for col in eachindex(rate)
                    S_data[row, col] = rate[col]
                    # Fill the symmetric upper triangular part
                    S_data[col, row] = rate[col]
                end
            elseif row == 22
                π = parse.(Float64, split(line))
            end
            row += 1
        end
    end

    return S_data, π
end

"""
    load_wag_rates()

Load WAG substitution rates matrix from WAG.dat.
Returns a 20×20 symmetric matrix of amino acid substitution rates.

# Citation
Whelan & Goldman (2001). DOI: 10.1093/oxfordjournals.molbev.a003851

# Returns
- `Matrix{Float64}`: 20×20 symmetric rate matrix
"""
function load_wag_rates()
    S, _ = load_empirical_data("../data/WAG.dat")
    return S
end

"""
    load_wag_frequencies()

Load WAG amino acid frequencies from WAG.dat.
Returns equilibrium frequencies in STANDARD_AA order.

# Citation
Whelan & Goldman (2001). DOI: 10.1093/oxfordjournals.molbev.a003851

# Returns
- `Vector{Float64}`: Vector of 20 equilibrium frequencies
"""
function load_wag_frequencies()
    _, π = load_empirical_data("../data/WAG.dat")
    return π
end

"""
    load_lg_rates()

Load LG substitution rates matrix from LG.dat.
Returns a 20×20 symmetric matrix of amino acid substitution rates.

# Citation
Le & Gascuel (2008). DOI: 10.1093/molbev/msn067

# Returns
- `Matrix{Float64}`: 20×20 symmetric rate matrix
"""
function load_lg_rates()
    S, _ = load_empirical_data("../data/LG.dat")
    return S
end

"""
    load_lg_frequencies()

Load LG amino acid frequencies from LG.dat.
Returns equilibrium frequencies in STANDARD_AA order.

# Citation
Le & Gascuel (2008). DOI: 10.1093/molbev/msn067

# Returns
- `Vector{Float64}`: Vector of 20 equilibrium frequencies
"""
function load_lg_frequencies()
    _, π = load_empirical_data("../data/LG.dat")
    return π
end