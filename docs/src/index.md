```@meta
CurrentModule = EvolutionModels
```

# EvolutionModels

Documentation for [EvolutionModels](https://github.com/mashu/EvolutionModels.jl).

> 🚧 **Work in Progress**: This package is under active development.

## Overview

EvolutionModels.jl implements continuous-time Markov chain (CTMC) models for molecular sequence evolution. The package supports:

- **DNA Models**: JC69, HKY85, GTR
- **Protein Models**: WAG, LG (empirical amino acid substitution matrices)
- **Simulation**: Evolve sequences along branches
- **Likelihood**: Compute probabilities of sequence pairs
- **Distance Estimation**: ML-based pairwise evolutionary distances

## Mathematical Foundation

### Continuous-Time Markov Chains

Sequence evolution is modeled as a continuous-time Markov process on the state space of nucleotides (A, C, G, T) or amino acids. The key mathematical objects are:

- **Q (Rate Matrix)**: Instantaneous rates of change, with Qᵢⱼ ≥ 0 for i ≠ j and row sums = 0
- **P(t) = exp(Qt)**: Transition probability matrix for evolution over time t
- **π (Stationary Distribution)**: Equilibrium frequencies satisfying πQ = 0

For time-reversible models, detailed balance holds:
```math
\pi_i Q_{ij} = \pi_j Q_{ji}
```

This is achieved by constructing Q from symmetric exchangeabilities S:
```math
Q_{ij} = S_{ij} \cdot \pi_j \quad \text{for } i \neq j
```

### Rate Normalization

By default, the overall evolutionary rate depends on μ. Use `normalize=true` to scale Q so that:
```math
-\sum_i \pi_i Q_{ii} = 1
```
This means branch length t directly represents the expected number of substitutions per site.

## DNA Substitution Models

### Jukes-Cantor (JC69)

The simplest model assuming:
- Equal base frequencies: πA = πC = πG = πT = 1/4
- Equal substitution rates between all nucleotides

Rate matrix:
```math
Q = \begin{pmatrix}
-3\alpha & \alpha & \alpha & \alpha \\
\alpha & -3\alpha & \alpha & \alpha \\
\alpha & \alpha & -3\alpha & \alpha \\
\alpha & \alpha & \alpha & -3\alpha
\end{pmatrix}
```
where α = μ/4 (or α = 1/3 when normalized)

### HKY85 Model

Extends JC69 with:
- Arbitrary base frequencies (πA, πC, πG, πT)
- Different rates for transitions (A↔G, C↔T) vs transversions
- Rate ratio κ between transitions and transversions

```math
q_{ij} = \begin{cases}
\mu\kappa\pi_j & \text{for transitions} \\
\mu\pi_j & \text{for transversions} \\
-\sum_{k\neq i} q_{ik} & \text{for } i = j
\end{cases}
```

### General Time-Reversible (GTR)

Most general neutral, reversible model:
- Arbitrary base frequencies π
- Six exchangeability parameters (rAC, rAG, rAT, rCG, rCT, rGT)
- Satisfies detailed balance

## Protein Substitution Models

### WAG Model

Empirical amino acid substitution matrix derived from 182 protein families (Whelan & Goldman, 2001). Provides both exchangeabilities and equilibrium frequencies estimated via maximum likelihood.

### LG Model

Improved empirical matrix from 3,912 protein families (Le & Gascuel, 2008). Generally recommended for modern protein sequence analysis.

## Basic Usage

### Creating Models

```julia
using EvolutionModels
using BioSequences

# JC69 - simplest DNA model
model_jc = create_model(JC69Model, 1.0)

# JC69 normalized (1 expected substitution per unit time)
model_jc_norm = create_model(JC69Model, 1.0, normalize=true)

# HKY85 with custom frequencies and transition bias
π = [0.3, 0.2, 0.2, 0.3]  # A,C,G,T
κ = 2.0  # transitions 2× faster than transversions
model_hky = create_model(HKY85Model, 1.0, π, κ, normalize=true)

# GTR with full parameterization
rates = [0.0 1.0 2.0 1.0;
         1.0 0.0 1.0 2.0;
         2.0 1.0 0.0 1.0;
         1.0 2.0 1.0 0.0]
model_gtr = create_model(GTRModel, 1.0, π, rates, normalize=true)

# Protein models
model_wag = create_model(WAGModel, 1.0, normalize=true)
model_lg = create_model(LGModel, 1.0, normalize=true)
```

### Sequence Evolution Simulation

```julia
# DNA evolution
model = create_model(JC69Model, 1.0, normalize=true)
seq = dna"ATCGATCGATCGATCG"
evolved = evolve_sequence(model, seq, 0.1)  # t=0.1 expected subs/site

# Protein evolution
model_aa = create_model(LGModel, 1.0, normalize=true)
protein = aa"EVQLVESGGGLVQPGGSLRL"
evolved_aa = evolve_sequence(model_aa, protein, 0.1)
```

### Likelihood Computation

```julia
model = create_model(JC69Model, 1.0, normalize=true)
seq1 = dna"ATCGATCG"
seq2 = dna"ATTGATCG"

# Compute log P(seq2 | seq1, t)
logL = sequence_likelihood(model, seq1, seq2, 0.1)
```

### Distance Matrix Computation

Requires FASTX.jl and Optim.jl:

```julia
using EvolutionModels
using FASTX
using Optim

# Read alignment
aln = read_alignment("sequences.fasta")

# Compute pairwise distances
model = create_model(JC69Model, 1.0, normalize=true)
result = compute_distances(model, aln)

# Access results
D = result.distances  # Distance matrix
labels = result.labels  # Sequence names
```

## Antibody Sequence Analysis

EvolutionModels.jl is well-suited for analyzing antibody evolution, including somatic hypermutation and B cell lineage analysis.

### Analyzing VH Region Evolution

```julia
using EvolutionModels
using BioSequences

# Use LG model for protein-level analysis
model = create_model(LGModel, 1.0, normalize=true)

# Germline VH sequence (IGHV3-23*01)
germline = aa"EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"

# Somatically mutated antibody
mutated = aa"EVQLVESGGGLVQPGRSLRLSCAASGFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAK"

# Compute likelihood at different evolutionary distances
for t in [0.01, 0.05, 0.10, 0.20]
    ll = sequence_likelihood(model, germline, mutated, t)
    println("t=$t: log-likelihood = $(round(ll, digits=2))")
end
```

### Simulating Somatic Hypermutation

```julia
using EvolutionModels
using BioSequences
using Random

Random.seed!(42)

model = create_model(LGModel, 1.0, normalize=true)

# Germline antibody
germline = aa"EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"

# Simulate evolution at increasing distances
println("Simulating somatic hypermutation:")
for t in [0.02, 0.05, 0.10, 0.15]
    evolved = evolve_sequence(model, germline, t)
    n_mut = count(g != e for (g, e) in zip(germline, evolved))
    pct = round(100 * n_mut / length(germline), digits=1)
    println("  t=$t: $n_mut mutations ($pct%)")
end
```

### DNA-Level Analysis with Transition Bias

Somatic hypermutation shows bias toward transitions. Model this with HKY85:

```julia
using EvolutionModels
using BioSequences

# HKY85 with transition bias (κ > 1)
π = [0.25, 0.25, 0.25, 0.25]
κ = 2.5  # Transitions favored
model = create_model(HKY85Model, 1.0, π, κ, normalize=true)

# VH nucleotide sequence
germline_nt = dna"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT"

# Simulate mutations
mutated_nt = evolve_sequence(model, germline_nt, 0.05)

# Analyze mutation spectrum
transitions = 0
transversions = 0
for (g, m) in zip(germline_nt, mutated_nt)
    if g != m
        if (g == DNA_A && m == DNA_G) || (g == DNA_G && m == DNA_A) ||
           (g == DNA_C && m == DNA_T) || (g == DNA_T && m == DNA_C)
            transitions += 1
        else
            transversions += 1
        end
    end
end
println("Transitions: $transitions, Transversions: $transversions")
println("Ti/Tv ratio: $(round(transitions/max(transversions,1), digits=2))")
```

### B Cell Lineage Distance Matrix

```julia
using EvolutionModels
using FASTX
using Optim

# Read aligned antibody sequences from a B cell lineage
# (germline + clonally-related sequences)
aln = read_alignment("bcell_lineage.fasta")

# Compute evolutionary distances
model = create_model(WAGModel, 1.0, normalize=true)
result = compute_distances(model, aln)

# Print formatted distance matrix
print_distance_matrix(result)

# Distances can be used for:
# - Phylogenetic tree reconstruction (e.g., with Phylogenetics.jl)
# - Hierarchical clustering
# - Identifying clonal relationships
```

### Interpreting Distances for Antibodies

When using normalized models, the evolutionary distance t represents expected amino acid substitutions per site:

| Distance (t) | Interpretation |
|-------------|----------------|
| 0.01-0.03 | Low SHM, early in affinity maturation |
| 0.05-0.10 | Moderate SHM, typical memory B cell |
| 0.10-0.20 | High SHM, extensively mutated |
| > 0.20 | Very high mutation load |

## Model Comparison

| Model | States | Parameters | Use Case |
|-------|--------|------------|----------|
| JC69 | DNA (4) | 1 (μ) | Quick estimates, equal rates |
| HKY85 | DNA (4) | 4+κ | Transition/transversion bias |
| GTR | DNA (4) | 4+6 | Full flexibility |
| WAG | Protein (20) | 1 (μ) | General protein evolution |
| LG | Protein (20) | 1 (μ) | Modern protein analysis |

## API Reference

```@autodocs
Modules = [EvolutionModels]
```
