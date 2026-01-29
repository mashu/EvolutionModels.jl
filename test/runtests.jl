using Test
using EvolutionModels
using BioSequences
using LinearAlgebra
using Random
using Statistics: std
import EvolutionModels: symbols
import EvolutionModels: is_transition
using FASTX
using Optim

# Helper function to check if matrix is symmetric in off-diagonal elements
is_symmetric_off_diagonal(matrix) = all(matrix[i, j] ≈ matrix[j, i]
    for i in 1:size(matrix, 1), j in 1:size(matrix, 2) if i != j)

import EvolutionModels: expected_substitution_rate

@testset verbose=true "EvolutionModels.jl" begin
    @testset "Types" begin
        @test length(DNAType()) == 4
        @test length(ProteinType()) == 20

        dna_type = DNAType()
        protein_type = ProteinType()

        @test length(symbols(dna_type)) == 4
        @test length(symbols(protein_type)) == 20

        @test all(nt isa DNA for nt in symbols(dna_type))
        @test all(aa isa AminoAcid for aa in symbols(protein_type))
    end

    @testset "DNA Models" begin
        @testset "JC69" begin
            μ = 0.1
            model = create_model(JC69Model, μ)

            @test model.μ == μ
            @test all(model.π .≈ 0.25)

            # Test R matrix structure
            @test size(model.R) == (4, 4)
            @test all(diag(model.R) .== 0.0)
            @test all(x -> x ≈ μ, filter(!iszero, model.R))
            @test issymmetric(model.R)
            @test sum(triu(model.R, 1)) ≈ 6μ

            # Test Q matrix properties
            @test all(isapprox.(sum(model.Q, dims=2), 0.0, atol=1e-14))
            @test issymmetric(model.Q ./ model.π')
            @test all(isapprox.(model.Q[1,2:4], μ/4, atol=1e-14))
            @test all(isapprox.(diag(model.Q), fill(-3μ/4, 4), atol=1e-14))
        end

        @testset "HKY85" begin
            μ = 0.1
            π = [0.3, 0.2, 0.2, 0.3]
            κ = 2.0
            model = create_model(HKY85Model, μ, π, κ)

            @test model.μ == μ
            @test model.π == π

            # Test R matrix structure
            @test size(model.R) == (4, 4)
            @test all(diag(model.R) .== 0.0)
            @test issymmetric(model.R)

            # Test transition/transversion rates
            dna_type = DNAType()
            for i in 1:4, j in i+1:4
                x, y = symbols(dna_type)[i], symbols(dna_type)[j]
                if is_transition(x, y)
                    @test model.R[i,j] ≈ κ * μ
                else
                    @test model.R[i,j] ≈ μ
                end
            end

            # Test Q matrix properties
            @test all(isapprox.(sum(model.Q, dims=2), 0.0, atol=1e-14))
            @test is_symmetric_off_diagonal(model.Q ./ model.π')
        end

        @testset "GTR" begin
            μ = 0.1
            π = [0.3, 0.2, 0.2, 0.3]
            rates = [0.0 1.0 2.0 1.0;
                    1.0 0.0 1.0 2.0;
                    2.0 1.0 0.0 1.0;
                    1.0 2.0 1.0 0.0]

            model = create_model(GTRModel, μ, π, rates)

            @test model.μ == μ
            @test model.π == π
            @test model.R == rates .* μ

            # Test Q matrix properties
            @test all(isapprox.(sum(model.Q, dims=2), 0.0, atol=1e-14))
            @test is_symmetric_off_diagonal(model.Q ./ model.π')
        end
    end

    @testset "Protein Models" begin
        @testset "WAG" begin
            μ = 0.1
            model = create_model(WAGModel, μ)

            @test model.μ == μ
            @test length(model.π) == 20
            @test isapprox(sum(model.π), 1.0, atol=1e-3)

            # Test matrix dimensions
            @test size(model.R) == (20, 20)
            @test size(model.Q) == (20, 20)

            # Test basic matrix properties
            @test all(diag(model.R) .== 0.0)
            @test issymmetric(model.R)
            @test all(isapprox.(sum(model.Q, dims=2), 0.0, atol=1e-14))
            @test is_symmetric_off_diagonal(model.Q ./ model.π')
        end

        @testset "LG" begin
            μ = 0.1
            model = create_model(LGModel, μ)

            @test model.μ == μ
            @test length(model.π) == 20
            @test isapprox(sum(model.π), 1.0, atol=1e-3)

            # Test matrix dimensions
            @test size(model.R) == (20, 20)
            @test size(model.Q) == (20, 20)

            # Test basic matrix properties
            @test all(diag(model.R) .== 0.0)
            @test issymmetric(model.R)
            @test all(isapprox.(sum(model.Q, dims=2), 0.0, atol=1e-14))
            @test is_symmetric_off_diagonal(model.Q ./ model.π')
        end
    end

    @testset "Normalization" begin
        @testset "JC69 Normalized" begin
            model = create_model(JC69Model, 1.0, normalize=true)
            rate = expected_substitution_rate(model.Q, model.π)
            @test isapprox(rate, 1.0, atol=1e-10)
        end

        @testset "HKY85 Normalized" begin
            π = [0.3, 0.2, 0.2, 0.3]
            model = create_model(HKY85Model, 1.0, π, 2.0, normalize=true)
            rate = expected_substitution_rate(model.Q, model.π)
            @test isapprox(rate, 1.0, atol=1e-10)
        end

        @testset "GTR Normalized" begin
            π = [0.3, 0.2, 0.2, 0.3]
            rates = [0.0 1.0 2.0 1.0;
                    1.0 0.0 1.0 2.0;
                    2.0 1.0 0.0 1.0;
                    1.0 2.0 1.0 0.0]
            model = create_model(GTRModel, 1.0, π, rates, normalize=true)
            rate = expected_substitution_rate(model.Q, model.π)
            @test isapprox(rate, 1.0, atol=1e-10)
        end

        @testset "WAG Normalized" begin
            model = create_model(WAGModel, 1.0, normalize=true)
            rate = expected_substitution_rate(model.Q, model.π)
            @test isapprox(rate, 1.0, atol=1e-10)
        end

        @testset "LG Normalized" begin
            model = create_model(LGModel, 1.0, normalize=true)
            rate = expected_substitution_rate(model.Q, model.π)
            @test isapprox(rate, 1.0, atol=1e-10)
        end

        @testset "Non-normalized rates differ" begin
            model_norm = create_model(JC69Model, 1.0, normalize=true)
            model_raw = create_model(JC69Model, 1.0, normalize=false)
            
            rate_norm = expected_substitution_rate(model_norm.Q, model_norm.π)
            rate_raw = expected_substitution_rate(model_raw.Q, model_raw.π)
            
            @test isapprox(rate_norm, 1.0, atol=1e-10)
            @test rate_raw ≈ 0.75  # 3μ/4 for JC69 with μ=1
        end
    end

    @testset "Gamma Rate Models" begin
        @testset "discrete_gamma_rates" begin
            rates, weights = discrete_gamma_rates(1.0, 4)
            @test length(rates) == 4
            @test length(weights) == 4
            @test all(weights .≈ 0.25)
            @test isapprox(sum(rates .* weights), 1.0, atol=0.1)  # Mean ≈ 1
            
            # Higher α = less variation
            rates_low_var, _ = discrete_gamma_rates(10.0, 4)
            rates_high_var, _ = discrete_gamma_rates(0.5, 4)
            @test std(rates_low_var) < std(rates_high_var)
        end

        @testset "GammaRateModel creation" begin
            base = create_model(LGModel, 1.0, normalize=true)
            model = create_gamma_model(base, 1.0)
            
            @test model.α == 1.0
            @test length(model.rates) == 4
            @test length(model.weights) == 4
            @test model.base_model === base
        end

        @testset "GammaRateModel likelihood" begin
            base = create_model(LGModel, 1.0, normalize=true)
            model = create_gamma_model(base, 1.0)
            
            seq1 = aa"EVQLVESGGGLVQPGGSLRL"
            seq2 = aa"EVQLVESGGGLIQPGGSLRL"
            
            # Gamma model likelihood should be computable
            L_gamma = sequence_likelihood(model, seq1, seq2, 0.1)
            L_base = sequence_likelihood(base, seq1, seq2, 0.1)
            
            @test isfinite(L_gamma)
            @test L_gamma < 0  # Log-likelihood is negative
            
            # Joint vs conditional
            L_joint = sequence_likelihood(model, seq1, seq2, 0.1, joint=true)
            L_cond = sequence_likelihood(model, seq1, seq2, 0.1, joint=false)
            @test L_joint < L_cond
        end

        @testset "GammaRateModel evolution" begin
            Random.seed!(42)
            base = create_model(LGModel, 1.0, normalize=true)
            model = create_gamma_model(base, 1.0)
            
            seq = aa"EVQLVESGGGLVQPGGSLRL"
            evolved = evolve_sequence(model, seq, 0.1)
            
            @test length(evolved) == length(seq)
            @test all(aa in STANDARD_AA for aa in evolved)
        end
    end

    @testset "Partition Models" begin
        @testset "PartitionModel creation" begin
            model1 = create_model(LGModel, 1.0, normalize=true)
            model2 = create_model(LGModel, 2.0, normalize=true)
            
            partition = create_partition_model(
                1:10 => model1,
                11:20 => model2
            )
            
            @test length(partition.partitions) == 2
            @test partition.total_length == 20
        end

        @testset "PartitionModel likelihood" begin
            model_slow = create_model(LGModel, 1.0, normalize=true)
            model_fast = create_model(LGModel, 1.0, normalize=true)
            
            partition = create_partition_model(
                1:10 => model_slow,
                11:20 => model_fast
            )
            
            seq1 = aa"EVQLVESGGGLVQPGGSLRL"
            seq2 = aa"EVQLVESGGGLIQPGGSLRL"
            
            L = sequence_likelihood(partition, seq1, seq2, 0.1)
            @test isfinite(L)
            @test L < 0
        end

        @testset "PartitionModel evolution" begin
            Random.seed!(42)
            model1 = create_model(LGModel, 1.0, normalize=true)
            model2 = create_model(LGModel, 2.0, normalize=true)
            
            partition = create_partition_model(
                1:10 => model1,
                11:20 => model2
            )
            
            seq = aa"EVQLVESGGGLVQPGGSLRL"
            evolved = evolve_sequence(partition, seq, 0.1)
            
            @test length(evolved) == length(seq)
            @test all(aa in STANDARD_AA for aa in evolved)
        end
    end

    @testset "Input Validation" begin
        # Invalid mutation rate
        @test_throws ArgumentError create_model(JC69Model, -0.1)
        @test_throws ArgumentError create_model(JC69Model, 0.0)
        @test_throws ArgumentError create_model(WAGModel, -0.1)

        # Invalid frequencies
        @test_throws ArgumentError create_model(HKY85Model, 0.1, [0.3, 0.2, 0.2], 2.0)
        @test_throws ArgumentError create_model(HKY85Model, 0.1, [0.3, 0.2, 0.2, 0.4], 2.0)
        @test_throws ArgumentError create_model(HKY85Model, 0.1, [-0.1, 0.3, 0.4, 0.4], 2.0)

        # Invalid κ
        @test_throws ArgumentError create_model(HKY85Model, 0.1, [0.25, 0.25, 0.25, 0.25], -1.0)
        @test_throws ArgumentError create_model(HKY85Model, 0.1, [0.25, 0.25, 0.25, 0.25], 0.0)

        # Invalid GTR rates
        @test_throws DimensionMismatch create_model(GTRModel, 0.1, [0.25, 0.25, 0.25, 0.25], ones(3,3))
    end

    @testset "Evolution" begin
        Random.seed!(42)  # For reproducibility

        @testset "DNA Evolution" begin
            model = create_model(JC69Model, 0.1)
            seq = dna"ATCGATCGATCG"

            evolved = evolve_sequence(model, seq, 1.0)
            @test length(evolved) == length(seq)
            @test all(nt in STANDARD_DNA for nt in evolved)

            # Test that very short time results in few changes
            evolved_short = evolve_sequence(model, seq, 0.0001)
            differences = count(x != y for (x,y) in zip(seq, evolved_short))
            @test differences < length(seq) ÷ 2

            # Test that long time approaches stationary distribution
            long_seq = dna"A"^1000
            evolved_long = evolve_sequence(model, long_seq, 1000.0)
            counts = zeros(Int, 4)
            for nt in evolved_long
                counts[findfirst(==(nt), STANDARD_DNA)] += 1
            end
            frequencies = counts ./ length(long_seq)
            @test all(isapprox.(frequencies, 0.25, atol=0.2))
        end

        @testset "Protein Evolution" begin
            model = create_model(WAGModel, 0.1)
            seq = aa"ARNDCQEGHILKMFPSTWYV"

            evolved = evolve_sequence(model, seq, 1.0)
            @test length(evolved) == length(seq)
            @test all(aa in STANDARD_AA for aa in evolved)

            # Test that very short time results in few changes
            evolved_short = evolve_sequence(model, seq, 0.0001)
            differences = count(x != y for (x,y) in zip(seq, evolved_short))
            @test differences < length(seq) ÷ 2
        end
    end

    @testset "Likelihood" begin
        @testset "DNA Likelihood" begin
            model = create_model(JC69Model, 0.1)
            seq1 = dna"ATCG"
            seq2 = dna"ATCG"  # Identical sequences
            seq3 = dna"GGGG"  # Very different sequence

            # Identical sequences should have higher likelihood at short times
            L_identical = sequence_likelihood(model, seq1, seq2, 0.1)
            L_different = sequence_likelihood(model, seq1, seq3, 0.1)
            @test L_identical > L_different

            # Test that likelihood decreases with time for identical sequences
            L_short = sequence_likelihood(model, seq1, seq1, 0.1)
            L_long = sequence_likelihood(model, seq1, seq1, 10.0)
            @test L_short > L_long
        end

        @testset "Protein Likelihood" begin
            model = create_model(WAGModel, 0.1)
            seq1 = aa"ARND"
            seq2 = aa"ARND"  # Identical sequences
            seq3 = aa"YYYY"  # Different sequence

            # Identical sequences should have higher likelihood at short times
            L_identical = sequence_likelihood(model, seq1, seq2, 0.1)
            L_different = sequence_likelihood(model, seq1, seq3, 0.1)
            @test L_identical > L_different
        end

        @testset "Joint vs Conditional Likelihood" begin
            model = create_model(LGModel, 1.0, normalize=true)
            seq1 = aa"ARNDCQEG"
            seq2 = aa"ARNDCQEG"

            L_cond = sequence_likelihood(model, seq1, seq2, 0.1)
            L_joint = sequence_likelihood(model, seq1, seq2, 0.1, joint=true)

            # Joint includes log(π) terms, so joint = cond + Σ log(πᵢ)
            π_contribution = sum(log(model.π[findfirst(==(aa), STANDARD_AA)]) 
                                 for aa in seq1 if aa in STANDARD_AA)
            @test isapprox(L_joint, L_cond + π_contribution, atol=1e-10)

            # Joint likelihood should be less than conditional (π < 1)
            @test L_joint < L_cond

            # Both should decrease with time for identical sequences
            L_cond_long = sequence_likelihood(model, seq1, seq2, 1.0)
            L_joint_long = sequence_likelihood(model, seq1, seq2, 1.0, joint=true)
            @test L_cond > L_cond_long
            @test L_joint > L_joint_long
        end

        # Test input validation
        @test_throws ArgumentError sequence_likelihood(
            create_model(JC69Model, 0.1),
            dna"AT", dna"ATCG", 0.1
        )
    end

    @testset "Interface Functions" begin
        model_dna = create_model(JC69Model, 0.1)
        model_protein = create_model(WAGModel, 0.1)

        # Test stationary frequencies
        @test length(stationary_frequencies(model_dna)) == 4
        @test length(stationary_frequencies(model_protein)) == 20
        @test all(isapprox.(sum(stationary_frequencies(model_dna)), 1.0))
        @test all(isapprox.(sum(stationary_frequencies(model_protein)), 1.0, atol=1e-3))

        # Test rate matrix
        @test size(rate_matrix(model_dna)) == (4, 4)
        @test size(rate_matrix(model_protein)) == (20, 20)
        @test issymmetric(rate_matrix(model_dna))
        @test issymmetric(rate_matrix(model_protein))

        # Test transition probability matrix
        P_dna = transition_probability_matrix(model_dna, 1.0)
        P_protein = transition_probability_matrix(model_protein, 1.0)

        @test size(P_dna) == (4, 4)
        @test size(P_protein) == (20, 20)
        @test all(isapprox.(sum(P_dna, dims=2), 1.0))
        @test all(isapprox.(sum(P_protein, dims=2), 1.0))
        @test all(P_dna .>= 0.0)
        @test all(P_protein .>= 0.0)
    end

    @testset "FASTX Extension" begin
        # Create a temporary FASTA file
        test_fasta = tempname() * ".fasta"
        open(test_fasta, "w") do io
            write(io, """>seq1
ATCGATCG
>seq2
ATCGATCC
>seq3
ATCGATCT
""")
        end

        try
            # Test read_alignment
            aln = read_alignment(test_fasta)
            @test length(aln.sequences) == 3
            @test length(aln.labels) == 3
            @test aln.length == 8
            @test aln.labels == ["seq1", "seq2", "seq3"]
            @test all(seq isa LongDNA{4} for seq in aln.sequences)

            # Test compute_distances with JC69 model
            model = create_model(JC69Model, 0.1)
            result = compute_distances(model, aln)

            # Test result structure
            @test result isa NamedTuple
            @test haskey(result, :distances)
            @test haskey(result, :labels)
            @test size(result.distances) == (3,3)
            @test result.labels == aln.labels

            # Test distance matrix properties
            D = result.distances
            @test issymmetric(D)
            @test all(isapprox.(diag(D), 0.0, atol=1e-10))
            @test all(D[i,j] > 0 for i in 1:3 for j in i+1:3)  # Off-diagonal elements should be positive

            # Test error handling
            @test_throws SystemError read_alignment("nonexistent.fasta")

            # Test unequal sequence lengths
            unequal_fasta = tempname() * ".fasta"
            open(unequal_fasta, "w") do io
                write(io, """>seq1
ATCG
>seq2
ATCGATCG
""")
            end
        finally
            rm(test_fasta)
        end

        @testset "alignment_likelihood" begin
            # Create a temporary FASTA file
            test_fasta = tempname() * ".fasta"
            open(test_fasta, "w") do io
                write(io, """>seq1
ATCGATCG
>seq2
ATCGATCC
>seq3
ATCGATCT
""")
            end

            try
                aln = read_alignment(test_fasta)
                model = create_model(JC69Model, 1.0, normalize=true)

                # Test alignment_likelihood
                result = alignment_likelihood(model, aln)
                
                @test haskey(result, :total_logL)
                @test haskey(result, :n_pairs)
                @test haskey(result, :mean_logL)
                @test haskey(result, :model)
                
                @test result.n_pairs == 3  # 3 choose 2 = 3 pairs
                @test result.total_logL < 0  # Log-likelihoods are negative
                @test isapprox(result.mean_logL, result.total_logL / result.n_pairs)
                
                # Joint likelihood should be less than conditional
                result_joint = alignment_likelihood(model, aln, joint=true)
                result_cond = alignment_likelihood(model, aln, joint=false)
                @test result_joint.total_logL < result_cond.total_logL
            finally
                rm(test_fasta)
            end
        end
    end
end
