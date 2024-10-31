```@meta
CurrentModule = EvolutionModels
```

# EvolutionModels

Documentation for [EvolutionModels](https://github.com/mashu/EvolutionModels.jl).

> 🚧 **Work in Progress**: This package is under active development.

## Introduction to DNA Substitution Models

### Jukes-Cantor (JC69)
The simplest model of DNA evolution, assuming:
- Equal base frequencies (πA = πC = πG = πT = 1/4)
- Equal substitution rates between all nucleotides
- Single parameter μ controlling overall mutation rate

Rate matrix Q has form:
```math
Q = \begin{pmatrix}
-3\alpha & \alpha & \alpha & \alpha \\
\alpha & -3\alpha & \alpha & \alpha \\
\alpha & \alpha & -3\alpha & \alpha \\
\alpha & \alpha & \alpha & -3\alpha
\end{pmatrix}
```
where α = μ/4

### HKY85 Model
Extends JC69 by introducing:
- Arbitrary base frequencies (πA, πC, πG, πT)
- Different rates for transitions (A↔G, C↔T) vs transversions
- Rate ratio κ between transitions and transversions

Rate matrix elements:
```math
q_{ij} = \begin{cases}
\mu\kappa\pi_j & \text{for transitions} \\
\mu\pi_j & \text{for transversions} \\
-\sum_{k\neq i} q_{ik} & \text{for } i = j
\end{cases}
```

### General Time-Reversible (GTR)
Most general neutral model with properties:
- Arbitrary base frequencies π
- Six substitution rate parameters (rAC, rAG, rAT, rCG, rCT, rGT)
- Satisfies detailed balance: πᵢqᵢⱼ = πⱼqⱼᵢ

Rate matrix Q = {qᵢⱼ} where:
```math
q_{ij} = \begin{cases}
\mu r_{ij}\pi_j & \text{for } i \neq j \\
-\sum_{k\neq i} q_{ik} & \text{for } i = j
\end{cases}
```

Note: Both JC69 and HKY85 are special cases of GTR with appropriate constraints on rates and frequencies.

## Usage

### Basic Example
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

### Distance Computation
```julia
using EvolutionModels
using FASTX
using Optim

# Read alignment
seqs = read_alignment("alignment.fasta")

# Create model and compute distances
model = create_model(JC69Model, 0.1)
result = compute_distances(model, seqs)

# Print distance matrix
print_distance_matrix(result)
```

### Using GTR Model
```julia
# Define base frequencies
π = [0.3, 0.2, 0.2, 0.3]  # A,C,G,T

# Define GTR rates
rates = GTRRates(
    1.0,  # A↔C
    2.0,  # A↔G
    1.0,  # A↔T
    1.0,  # C↔G
    2.0,  # C↔T
    1.0   # G↔T
)

# Create model
model = create_model(GTRModel, 0.1, π, rates)
```

## API Reference

```@autodocs
Modules = [EvolutionModels]
```
