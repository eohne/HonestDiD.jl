# Delta matrix builder checks. The builders were validated byte-for-byte against
# R across 149 configurations during development; here we keep a representative
# set with inline expected matrices (plus the numPre==1 bug fix). End-to-end
# exactness of the RM/SDRM builders is covered by validate_variants.jl.
using HonestDiD, Test

const H = HonestDiD

@testset "delta builders" begin
    @testset "create_A_SD / create_d_SD" begin
        @test H.create_A_SD(2, 2) == [1.0 -2 0 0; 0 1 1 0; 0 0 -2 1; -1 2 0 0; 0 -1 -1 0; 0 0 2 -1]
        @test H.create_d_SD(2, 2, 0.1) == fill(0.1, 6)
        @test size(H.create_A_SD(4, 4)) == (14, 8)            # 2*(numPre+numPost-1)
    end

    @testset "create_A_B (sign restriction)" begin
        @test H.create_A_B(3, 2, "positive") == [0.0 0 0 -1 0; 0 0 0 0 -1]
        @test H.create_A_B(3, 2, "negative") == [0.0 0 0 1 0; 0 0 0 0 1]
    end

    @testset "create_A_M (monotonicity)" begin
        @test H.create_A_M(2, 2, "increasing", false) == [1.0 -1 0 0; 0 1 0 0; 0 0 -1 0; 0 0 1 -1]
        @test H.create_A_M(2, 2, "decreasing", false) == -[1.0 -1 0 0; 0 1 0 0; 0 0 -1 0; 0 0 1 -1]
    end

    @testset "create_A_M numPre==1 (R bug fix)" begin
        # δ_pre ≤ 0 (no stray coupling from R's 1:0 loop).
        @test H.create_A_M(1, 2, "increasing", false) == [1.0 0 0; 0 -1 0; 0 1 -1]
    end

    @testset "B/M variants append correctly" begin
        A_SDB = H.create_A_SDB(3, 2; biasDirection="positive")
        A_SD = H.create_A_SD(3, 2)
        @test A_SDB[1:size(A_SD, 1), :] == A_SD
        @test A_SDB[(size(A_SD, 1)+1):end, :] == H.create_A_B(3, 2, "positive")
        @test H.create_d_SDB(3, 2, 0.1) == vcat(fill(0.1, size(A_SD, 1)), zeros(2))
    end

    @testset "RM/SDRM builders: no all-zero rows, d all zero" begin
        for s in -2:0
            A = H.create_A_RM(3, 2; Mbar=1.5, s=s, max_positive=true)
            @test all(sum(abs2, A[r, :]) > 1e-10 for r in 1:size(A, 1))
            @test size(A, 2) == 5
        end
        @test all(==(0), H.create_d_RM(3, 2))
        @test all(==(0), H.create_d_SDRM(3, 2))
    end
end
