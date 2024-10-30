# EvolutionModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/EvolutionModels.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/EvolutionModels.jl/dev/)
[![Build Status](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mashu/EvolutionModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mashu/EvolutionModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/EvolutionModels.jl)

A Julia package implementing continuous-time Markov models for molecular sequence evolution. Provides both parametric nucleotide models (JC69, HKY85, GTR) and empirical amino acid models (WAG, LG) with simulation and likelihood computation capabilities.

## Features

- Nucleotide substitution models: JC69, HKY85, GTR
- Amino acid substitution models: WAG, LG
- Sequence evolution simulation
- Likelihood computation
- Distance matrix computation (requires FASTX.jl)
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

The package can compute pairwise evolutionary distances from multiple sequence alignments in FASTA format. This functionality requires FASTX.jl:

```julia
using EvolutionModels
using FASTX

# Read alignment
seqs = read_alignment("alignment.fasta")

# Create model and compute distances
model = create_model(JC69Model, 0.1)
D = compute_distances(model, seqs)
```

The resulting distance matrix `D` contains maximum likelihood estimates of evolutionary times between sequences under the specified model. This can be used for:
- Phylogenetic tree reconstruction
- Sequence clustering
- Evolutionary rate estimation

## License

MIT
