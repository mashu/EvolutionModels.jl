# EvolutionModels.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/EvolutionModels.jl/dev/)
[![Build Status](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/EvolutionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/EvolutionModels.jl)
[![Experimental](https://img.shields.io/badge/status-experimental-orange.svg)](https://github.com/mashu/EvolutionModels.jl)

> **Warning**
> 
> This package is a **personal experiment** and is **not ready for production use**. The implementation may contain errors, the API is unstable, and results have not been independently validated. Do not use for research, clinical, or any application where correctness matters. If you need reliable phylogenetic tools, consider established software.

Continuous-time Markov models for molecular sequence evolution in Julia.

## Features

- **DNA models**: JC69, HKY85, GTR
- **Protein models**: WAG, LG
- Sequence evolution simulation
- Pairwise likelihood computation
- ML distance estimation
- Integration with BioSequences.jl and FASTX.jl

## Installation

```julia
using Pkg
Pkg.add("EvolutionModels")
```

## Quick Start

```julia
using EvolutionModels
using BioSequences

# Create a model
model = create_model(LGModel, 1.0, normalize=true)

# Simulate evolution
seq = aa"EVQLVESGGGLVQPGGSLRL"
evolved = evolve_sequence(model, seq, 0.1)

# Compute likelihood
logL = sequence_likelihood(model, seq, evolved, 0.1)
```

## Documentation

See the [documentation](https://mashu.github.io/EvolutionModels.jl/dev/) for model selection guidance, distance interpretation, and antibody sequence analysis examples.

## License

MIT
