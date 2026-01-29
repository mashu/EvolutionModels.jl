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
| Protein (simple) | LG | `create_model(LGModel, 1.0, normalize=true)` |
| Protein (rate variation) | LG+G | `create_gamma_model(base_lg, 1.0)` |
| Antibodies | LG+G or Partition | See Antibody section |
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

### Recommended Models for Antibodies

For antibody sequences, we recommend using rate variation models since CDR and framework regions evolve at very different rates:

```julia
# Option 1: LG+G (simplest, good default)
base = create_model(LGModel, 1.0, normalize=true)
model = create_gamma_model(base, 0.8)  # α=0.8 captures CDR/FR heterogeneity

# Option 2: Partition model (if you know CDR boundaries)
framework = create_model(LGModel, 1.0, normalize=true)
cdr = create_model(LGModel, 2.5, normalize=true)  # CDRs evolve faster
model = create_partition_model(
    1:26 => framework, 27:38 => cdr,      # FR1, CDR1
    39:55 => framework, 56:65 => cdr,     # FR2, CDR2
    66:104 => framework, 105:117 => cdr,  # FR3, CDR3
    118:128 => framework                   # FR4
)
```

### Measuring Somatic Hypermutation

Evolutionary distance between germline and mutated antibody sequences quantifies the extent of somatic hypermutation (SHM).

```julia
using EvolutionModels
using BioSequences
using Optim

# Use Gamma model for better handling of rate variation
base = create_model(LGModel, 1.0, normalize=true)
model = create_gamma_model(base, 1.0)

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

### Conditional vs Joint Likelihood

The `sequence_likelihood` function supports two modes:

**Conditional likelihood** (default):
```math
\log P(\text{seq2} | \text{seq1}, t) = \sum_{k=1}^{L} \log P(t)_{x_{1k}, x_{2k}}
```

**Joint likelihood** (`joint=true`):
```math
\log P(\text{seq1}, \text{seq2} | t) = \sum_{k=1}^{L} \left[ \log \pi_{x_{1k}} + \log P(t)_{x_{1k}, x_{2k}} \right]
```

```julia
model = create_model(LGModel, 1.0, normalize=true)
seq1 = aa"EVQLVESGGGLVQPGGSLRL"
seq2 = aa"EVQLVESGGGLIQPGGSLRL"

# Conditional: P(seq2 | seq1, t)
L_cond = sequence_likelihood(model, seq1, seq2, 0.05)

# Joint: P(seq1, seq2 | t)  
L_joint = sequence_likelihood(model, seq1, seq2, 0.05, joint=true)
```

### When Joint Likelihood Matters

#### Use Conditional (default) When:

**Estimating evolutionary distance between two sequences**

The stationary frequency term ∑ log(πᵢ) is constant for fixed sequences regardless of t. Therefore, finding the distance that maximizes likelihood gives the same answer whether you use conditional or joint:

```julia
# Both give the same optimal distance
argmax_t P(seq2|seq1,t) = argmax_t P(seq1,seq2|t)
```

For distance estimation, conditional is computationally simpler and equally correct.

**Comparing likelihoods at different times for the same sequence pair**

If you're asking "is t=0.05 or t=0.10 more likely for these sequences?", the π term cancels out in the comparison.

#### Use Joint Likelihood When:

**Comparing different models**

When comparing models with different stationary frequencies (e.g., LG vs WAG, or models with estimated π), the joint likelihood is required for proper comparison:

```julia
model_lg = create_model(LGModel, 1.0, normalize=true)
model_wag = create_model(WAGModel, 1.0, normalize=true)

# WRONG: Conditional doesn't account for different π
L_lg_cond = sequence_likelihood(model_lg, seq1, seq2, 0.1)
L_wag_cond = sequence_likelihood(model_wag, seq1, seq2, 0.1)

# CORRECT: Joint includes π contribution
L_lg_joint = sequence_likelihood(model_lg, seq1, seq2, 0.1, joint=true)
L_wag_joint = sequence_likelihood(model_wag, seq1, seq2, 0.1, joint=true)

# Model with higher joint likelihood better explains the data
better_model = L_lg_joint > L_wag_joint ? "LG" : "WAG"
```

**Computing AIC/BIC for model selection**

Information criteria require the actual likelihood, not conditional probability:

```julia
# AIC = 2k - 2ln(L) where L is the likelihood
k = 1  # number of parameters (just distance for empirical models)
AIC_lg = 2*k - 2*L_lg_joint
AIC_wag = 2*k - 2*L_wag_joint
```

**Bayesian inference**

Posterior probabilities require the joint likelihood:

```julia
# P(model|data) ∝ P(data|model) × P(model)
# P(data|model) is the joint likelihood
```

**Comparing sequences with different compositions**

If comparing likelihood across different sequence pairs, joint likelihood accounts for how "probable" each sequence is under the model's equilibrium:

```julia
# Sequence pair A (common amino acids)
seq1a = aa"AAAAAAAAA"
seq2a = aa"AAAAAAAAA"

# Sequence pair B (rare amino acids)  
seq1b = aa"WWWWWWWWW"
seq2b = aa"WWWWWWWWW"

# Conditional treats both equally (identical sequences)
L_a_cond = sequence_likelihood(model, seq1a, seq2a, 0.01)  # High
L_b_cond = sequence_likelihood(model, seq1b, seq2b, 0.01)  # High

# Joint reflects that W is rare in natural proteins
L_a_joint = sequence_likelihood(model, seq1a, seq2a, 0.01, joint=true)
L_b_joint = sequence_likelihood(model, seq1b, seq2b, 0.01, joint=true)
# L_a_joint > L_b_joint because A is more common than W
```

#### Lineage Reconstruction and Clustering

For B cell lineage reconstruction or antibody clustering, use **conditional likelihood** (the default). This is what `compute_distances()` does internally.

When computing pairwise distances for a distance matrix:
1. Each pair gets its own ML distance estimate
2. The π term doesn't affect which distance is optimal for each pair
3. All distances are in the same units and directly comparable

```julia
# Correct approach for lineage/clustering
model = create_model(LGModel, 1.0, normalize=true)
result = compute_distances(model, alignment)
D = result.distances  # Ready for tree building or clustering
```

#### Comparing Models for a Lineage

If you want to ask "does LG or WAG better explain my antibody lineage?", use `alignment_likelihood` to compute total likelihood over all pairs:

```julia
using EvolutionModels
using FASTX
using Optim

aln = read_alignment("lineage.fasta")

model_lg = create_model(LGModel, 1.0, normalize=true)
model_wag = create_model(WAGModel, 1.0, normalize=true)

# Compute total joint log-likelihood for each model
result_lg = alignment_likelihood(model_lg, aln)
result_wag = alignment_likelihood(model_wag, aln)

println("LG:  total=$(round(result_lg.total_logL, digits=2)), mean=$(round(result_lg.mean_logL, digits=2))")
println("WAG: total=$(round(result_wag.total_logL, digits=2)), mean=$(round(result_wag.mean_logL, digits=2))")
println("Better model: ", result_lg.total_logL > result_wag.total_logL ? "LG" : "WAG")
```

This answers "which substitution model better captures the evolutionary patterns in my antibody data?" - useful when you have lineage-specific amino acid preferences or want to validate model choice.

### Summary

| Question | Use | Reason |
|----------|-----|--------|
| "What is the evolutionary distance?" | Conditional | π doesn't affect optimal t |
| "Distance matrix for clustering/trees?" | Conditional | Each pair optimized independently |
| "Which model fits better?" | Joint | Must account for different π |
| "Is this sequence pair unusual?" | Joint | Accounts for amino acid frequencies |
| "What's the AIC/BIC?" | Joint | Information criteria need true likelihood |
| "Posterior probability of model?" | Joint | Bayesian inference requirement |

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

## Advanced Models

### Gamma Rate Variation (+G Models)

Real sequences have sites that evolve at different rates - some positions are highly conserved while others are hypervariable. The Gamma model accounts for this by drawing site rates from a discretized Gamma distribution.

```julia
# Create base model
base = create_model(LGModel, 1.0, normalize=true)

# Add Gamma rate variation with 4 categories
# α controls rate variation: smaller = more variation
model = create_gamma_model(base, 0.5)   # High variation
model = create_gamma_model(base, 1.0)   # Moderate variation  
model = create_gamma_model(base, 2.0)   # Low variation
```

**When to use:**
- Sequences with mix of conserved and variable sites
- Antibody analysis (CDRs are hypervariable, frameworks conserved)
- Any analysis where assuming uniform rates seems unrealistic

**Parameter α guidance:**
| α value | Rate variation | Typical use |
|---------|---------------|-------------|
| 0.3-0.5 | Very high | Highly heterogeneous data |
| 0.5-1.0 | High | Antibodies, proteins with functional constraints |
| 1.0-2.0 | Moderate | General protein analysis |
| >2.0 | Low | Relatively homogeneous rates |

### Partition Models (Region-Specific Evolution)

For antibodies, CDR and framework regions evolve differently. Partition models apply different evolutionary parameters to different sequence regions.

```julia
# Different models for different regions
framework = create_model(LGModel, 1.0, normalize=true)
cdr = create_model(LGModel, 2.0, normalize=true)  # 2× rate for CDRs

# IMGT numbering example for VH
model = create_partition_model(
    1:26 => framework,     # FR1
    27:38 => cdr,          # CDR1
    39:55 => framework,    # FR2
    56:65 => cdr,          # CDR2
    66:104 => framework,   # FR3
    105:117 => cdr,        # CDR3
    118:128 => framework   # FR4
)

# Use like any other model
seq1 = aa"..."  # Your antibody sequence
seq2 = aa"..."
logL = sequence_likelihood(model, seq1, seq2, 0.1)
```

**When to use:**
- Antibody sequences with defined CDR/framework boundaries
- Multi-domain proteins with different evolutionary pressures
- Any sequence with distinct functional regions

### Combining Approaches

You can use Gamma models as the base for partition models:

```julia
# Framework: conserved, low rate variation
framework_base = create_model(LGModel, 1.0, normalize=true)
framework = create_gamma_model(framework_base, 2.0)  # α=2, low variation

# CDR: hypervariable, high rate variation
cdr_base = create_model(LGModel, 2.0, normalize=true)
cdr = create_gamma_model(cdr_base, 0.5)  # α=0.5, high variation

model = create_partition_model(
    1:26 => framework,
    27:38 => cdr,
    # ... etc
)
```

---

## Limitations and Considerations

### Model Assumptions

**Basic models** (JC69, HKY85, GTR, WAG, LG) assume:
- Site independence: Each position evolves independently
- Time reversibility: Process is stationary and reversible
- Homogeneous rates: Same process at all sites

**Gamma models** (+G) relax the homogeneous rates assumption by allowing site-specific rate variation.

**Partition models** allow region-specific parameters but still assume independence within each region.

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
