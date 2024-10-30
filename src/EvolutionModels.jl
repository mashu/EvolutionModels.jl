module EvolutionModels

using LinearAlgebra
using BioSequences
using Random
using StatsBase

# Extension interface functions that will be implemented by extensions
function read_alignment end
function compute_distances end

# Export the interface functions
export read_alignment, compute_distances

include("types.jl")
include("utils.jl")
include("parametric.jl")
include("empirical.jl")
include("evolution.jl")

# Export types
export SequenceType, DNAType, ProteinType
export EvolutionModel, DNAModel, ProteinModel, Model

# Export model types
export JC69Model, HKY85Model, GTRModel, WAGModel, LGModel

# Export constants
export STANDARD_DNA, STANDARD_AA

# Export functions
export create_model
export transition_probability_matrix
export stationary_frequencies, rate_matrix
export evolve_sequence, sequence_likelihood
export symbols  # Export symbols function

# Interface functions
"""
Get the stationary distribution π.
"""
stationary_frequencies(model::Model) = model.π

"""
Get the rate matrix R (includes μ scaling).
"""
rate_matrix(model::Model) = model.R

"""
Compute transition probability matrix P(t) = exp(Qt).
"""
transition_probability_matrix(model::Model, t::Float64) = exp(model.Q * t)

end # module
