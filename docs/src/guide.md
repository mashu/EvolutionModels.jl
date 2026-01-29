# Practical Guide

This guide covers model selection, interpretation of results, and practical applications for antibody sequence analysis.

## Model Selection

### DNA Models

#### JC69 (Jukes-Cantor 1969)

The simplest nucleotide substitution model with equal rates between all bases and uniform base frequencies (π = 0.25 for all).

**Use when:**
- Quick distance estimates are needed
- Sequences are closely related (low divergence)
- Base composition is approximately uniform

```julia
model = create_model(JC69Model, 1.0, normalize=true)
```

#### HKY85 (Hasegawa-Kishino-Yano 1985)

Extends JC69 with a transition/transversion rate ratio (κ) and arbitrary base frequencies.

**Use when:**
- Transition bias is expected (common in most biological data)
- Base frequencies differ from uniform
- Analyzing somatic hypermutation in antibodies (κ ≈ 2-3 captures AID bias)

```julia
π = [0.25, 0.25, 0.25, 0.25]  # A, C, G, T frequencies
κ = 2.5  # Transition/transversion ratio
model = create_model(HKY85Model, 1.0, π, κ, normalize=true)
```

#### GTR (General Time-Reversible)

The most general neutral, reversible model with 6 exchangeability parameters and arbitrary base frequencies.

**Use when:**
- Full flexibility is required
- Model selection indicates simpler models are inadequate
- Sufficient data exists to estimate additional parameters

```julia
π = [0.3, 0.2, 0.2, 0.3]
rates = [0.0 1.0 2.0 1.0;
         1.0 0.0 1.0 2.0;
         2.0 1.0 0.0 1.0;
         1.0 2.0 1.0 0.0]
model = create_model(GTRModel, 1.0, π, rates, normalize=true)
```

### Protein Models

#### LG (Le-Gascuel 2008)

Empirical amino acid substitution matrix derived from 3,912 protein families. Generally the recommended default for protein sequence analysis.

**Use when:**
- Analyzing protein sequences (default choice)
- Antibody VH/VL region analysis
- General protein evolution studies

```julia
model = create_model(LGModel, 1.0, normalize=true)
```

#### WAG (Whelan-Goldman 2001)

Earlier empirical matrix from 182 protein families. Provides an alternative to LG with slightly different exchangeabilities.

**Use when:**
- Alternative to LG for comparison
- Consistency with older analyses that used WAG

```julia
model = create_model(WAGModel, 1.0, normalize=true)
```

### Recommended Defaults

| Data Type | Model | Configuration |
|-----------|-------|---------------|
| Protein sequences | LG | `create_model(LGModel, 1.0, normalize=true)` |
| DNA (general) | HKY85 | `create_model(HKY85Model, 1.0, π, 2.0, normalize=true)` |
| DNA (quick estimate) | JC69 | `create_model(JC69Model, 1.0, normalize=true)` |

---

## Understanding Distances

### Normalized vs Unnormalized Models

The `normalize=true` option scales the rate matrix so that the expected number of substitutions per unit time equals 1. This makes distances directly interpretable.

```julia
# Normalized: distance = expected substitutions per site
model = create_model(LGModel, 1.0, normalize=true)

# Verify normalization
rate = expected_substitution_rate(model.Q, model.π)  # ≈ 1.0
```

**With normalization**: A distance of 0.1 means approximately 10% of sites have experienced at least one substitution.

**Without normalization**: Distances depend on the model's inherent scaling and cannot be directly compared across models.

### Distance Interpretation

| Distance | Sequence Divergence | Typical Scenario |
|----------|--------------------|--------------------|
| 0.00-0.01 | Nearly identical | Recent duplication, sequencing of same clone |
| 0.01-0.05 | Low divergence | Closely related sequences, early divergence |
| 0.05-0.15 | Moderate divergence | Typical within-species variation |
| 0.15-0.50 | High divergence | Between-species comparisons |
| >0.50 | Saturated | Multiple substitutions obscure signal |

### Handling Saturated Sequences

For highly divergent sequences (p-distance > 0.7), distance estimates become unreliable due to multiple substitutions at the same site. Consider:

1. Using protein sequences instead of DNA (larger alphabet reduces saturation)
2. Focusing on conserved regions
3. Interpreting distances with caution

---

## Antibody Sequence Analysis

### Measuring Somatic Hypermutation

Evolutionary distance between germline and mutated antibody sequences quantifies the extent of somatic hypermutation (SHM).

```julia
using EvolutionModels
using BioSequences
using Optim

model = create_model(LGModel, 1.0, normalize=true)

germline = aa"EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
mutated  = aa"EVQLVESGGGLVQPGRSLRLSCAASGFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAK"

# Maximum likelihood distance estimation
neg_ll(t) = t < 0 ? Inf : -sequence_likelihood(model, germline, mutated, t)
result = optimize(neg_ll, 0.0, 1.0, Brent())
distance = Optim.minimizer(result)

println("SHM distance: $(round(distance, digits=4))")
```

### SHM Distance Interpretation

| Distance | SHM Level | Biological Context |
|----------|-----------|-------------------|
| 0.00-0.02 | Minimal | Near-germline, naive B cell |
| 0.02-0.05 | Low | Early memory, limited affinity maturation |
| 0.05-0.10 | Moderate | Typical memory B cell |
| 0.10-0.15 | High | Extensively matured, long-lived plasma cell |
| >0.15 | Very high | Chronic infection, autoimmune conditions |

### B Cell Lineage Analysis

Computing pairwise distances for clonally-related sequences enables lineage reconstruction and subclone identification.

```julia
using EvolutionModels
using FASTX
using Optim

# Load aligned sequences from a B cell lineage
aln = read_alignment("bcell_lineage.fasta")

model = create_model(LGModel, 1.0, normalize=true)
result = compute_distances(model, aln)

# Distance matrix for downstream analysis
D = result.distances
labels = result.labels

print_distance_matrix(result, digits=3)
```

The resulting distance matrix can be used for:
- Hierarchical clustering to identify subclones
- Input to tree reconstruction algorithms
- Quantifying within-lineage diversity

### DNA-Level Analysis

For nucleotide-level analysis capturing SHM transition bias:

```julia
using EvolutionModels
using BioSequences

# HKY85 with transition bias typical of AID-mediated SHM
π = [0.25, 0.25, 0.25, 0.25]
κ = 2.5  # Transitions favored over transversions
model = create_model(HKY85Model, 1.0, π, κ, normalize=true)

germline_nt = dna"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAG"
mutated_nt  = dna"GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGCAG"

# Compute likelihood at estimated distance
distance = 0.05
ll = sequence_likelihood(model, germline_nt, mutated_nt, distance)
```

### Simulating Sequence Evolution

Generate hypothetical evolved sequences for method validation or hypothesis testing:

```julia
using EvolutionModels
using BioSequences
using Random

Random.seed!(42)

model = create_model(LGModel, 1.0, normalize=true)
germline = aa"EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"

# Simulate evolution at different distances
for t in [0.02, 0.05, 0.10, 0.15]
    evolved = evolve_sequence(model, germline, t)
    n_mut = count(g != e for (g, e) in zip(germline, evolved))
    println("t=$t: $n_mut mutations ($(round(100*n_mut/length(germline), digits=1))%)")
end
```

---

## Likelihood Computation

### Conditional Probability

The `sequence_likelihood` function computes the conditional log-probability:

```math
\log P(\text{seq2} | \text{seq1}, t) = \sum_{k=1}^{L} \log P(t)_{x_{1k}, x_{2k}}
```

where P(t) = exp(Qt) is the transition probability matrix.

```julia
model = create_model(LGModel, 1.0, normalize=true)
seq1 = aa"EVQLVESGGGLVQPGGSLRL"
seq2 = aa"EVQLVESGGGLIQPGGSLRL"

# Log-probability of seq2 given seq1 evolved for time t
logL = sequence_likelihood(model, seq1, seq2, 0.05)
```

### Maximum Likelihood Distance Estimation

Find the evolutionary distance that maximizes the likelihood:

```julia
using Optim

model = create_model(LGModel, 1.0, normalize=true)
seq1 = aa"EVQLVESGGGLVQPGGSLRL"
seq2 = aa"EVQLVESGGGLIQPGGSLRL"

neg_ll(t) = t < 0 ? Inf : -sequence_likelihood(model, seq1, seq2, t)
result = optimize(neg_ll, 0.0, 2.0, Brent())

ml_distance = Optim.minimizer(result)
max_logL = -Optim.minimum(result)
```

### Comparing Evolutionary Hypotheses

Use likelihood ratios to compare different evolutionary scenarios:

```julia
model = create_model(LGModel, 1.0, normalize=true)

# Two sequences
seqA = aa"EVQLVESGGGLVQPGGSLRL"
seqB = aa"EVQLVESGGGLIQPGGSLRL"

# Compare hypotheses: recent vs ancient divergence
L_recent = sequence_likelihood(model, seqA, seqB, 0.02)
L_ancient = sequence_likelihood(model, seqA, seqB, 0.20)

log_ratio = L_recent - L_ancient
println("Log likelihood ratio: $log_ratio")
# Positive = recent divergence more likely
```

---

## Limitations and Considerations

### Model Assumptions

All models assume:
- **Site independence**: Each position evolves independently
- **Time reversibility**: Process is stationary and reversible
- **Homogeneous rates**: Same process at all sites (no rate variation)

These assumptions may be violated in biological data, particularly for antibodies where CDR regions experience selection.

### Empirical Models and SHM

The WAG and LG matrices were derived from germline protein evolution across diverse protein families. Somatic hypermutation has distinct characteristics:
- AID-mediated with specific hotspot motifs
- Strong selection during affinity maturation
- Rapid timescale (days to weeks vs millions of years)

The empirical models provide useful approximations but may not perfectly capture SHM patterns. For specialized applications, consider SHM-specific models from the literature.

### Gaps and Missing Data

Non-standard characters (gaps, ambiguous bases) are excluded from likelihood calculations. Heavily gapped alignments may produce unreliable results. Consider:
- Quality filtering before analysis
- Analyzing gap patterns separately
- Using alignment-free methods for highly divergent sequences

---

## References

1. Jukes, T.H. & Cantor, C.R. (1969). Evolution of Protein Molecules. *Mammalian Protein Metabolism*, 3:21-132.

2. Hasegawa, M., Kishino, H. & Yano, T. (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. *J Mol Evol*, 22:160-174.

3. Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17:57-86.

4. Whelan, S. & Goldman, N. (2001). A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Mol Biol Evol*, 18:691-699.

5. Le, S.Q. & Gascuel, O. (2008). An improved general amino acid replacement matrix. *Mol Biol Evol*, 25:1307-1320.

6. Yang, Z. (2014). *Molecular Evolution: A Statistical Approach*. Oxford University Press.

7. Yaari, G. & Kleinstein, S.H. (2015). Practical guidelines for B-cell receptor repertoire sequencing analysis. *Genome Med*, 7:121.
