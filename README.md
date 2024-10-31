# EvolutionModels.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/EvolutionModels.jl/dev/)
[![Build Status](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/EvolutionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/EvolutionModels.jl)

A Julia package implementing continuous-time Markov models for molecular sequence evolution. Provides both parametric nucleotide models (JC69, HKY85, GTR) and empirical amino acid models (WAG, LG) with simulation and likelihood computation capabilities.

> ⚠️ **Note: This package is currently under development and extension (loaded with FASTX and Optim) is experimental. Use with caution in production environments.**

## Features

- Nucleotide substitution models: JC69, HKY85, GTR
- Amino acid substitution models: WAG, LG
- Sequence evolution simulation
- Likelihood computation
- Distance matrix computation (requires FASTX.jl and Optim.jl)
- Integration with BioSequences.jl

## Installation

```julia
using Pkg
Pkg.add("EvolutionModels")
```

## Basic Usage

```julia
using EvolutionModels
using BioSequences

# Create JC69 model
model = create_model(JC69Model, 0.1)

# Simulate evolution
seq = dna"ATCG"
evolved = evolve_sequence(model, seq, 1.0)

# Compute likelihood
logL = sequence_likelihood(model, seq, evolved, 1.0)
```

## Computing Distances from FASTA Alignment

The package can compute pairwise evolutionary distances from multiple sequence alignments in FASTA format. This functionality requires both FASTX.jl and Optim.jl:

```julia
using EvolutionModels
using FASTX
using Optim

# Read alignment
seqs = read_alignment("alignment.fasta")

# Create model and compute distances
model = create_model(JC69Model, 0.1)
result = compute_distances(model, seqs)

# Print the distance matrix
print_distance_matrix(result)
```

The resulting distance matrix contains maximum likelihood estimates of evolutionary times between sequences under the specified model. This can be used for:
- Phylogenetic tree reconstruction
- Sequence clustering
- Evolutionary rate estimation

## License

MIT

