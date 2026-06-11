# End-to-end validation against R ground truth (embedded in reference.jl).
using HonestDiD, Test
import Tables
isdefined(Main, :BC_BETAHAT) || include(joinpath(@__DIR__, "reference.jl"))

bh, sg, np, nq = BC_BETAHAT, BC_SIGMA, BC_NPRE, BC_NPOST

@testset "HonestDiD vs R (BCdata)" begin
    @testset "test_delta_lp_fn (exact LP + duals)" begin
        r = HonestDiD.test_delta_lp_fn(LP_Y, LP_X, LP_SIGMA)
        @test isapprox(r.eta_star, LP_ETA; atol=1e-7)
        @test isapprox(r.delta_star, LP_DELTA; atol=1e-6)
        @test isapprox(r.lambda, LP_LAMBDA; atol=1e-6)
    end

    @testset "constructOriginalCS (exact)" begin
        oc = constructOriginalCS(bh, sg, np, nq)
        @test isapprox(oc.lb[1], ORIG_LB; atol=1e-9)
        @test isapprox(oc.ub[1], ORIG_UB; atol=1e-9)
    end

    @testset "findOptimalFLCI (analytic qfoldednormal, atol 1e-3)" begin
        for r in axes(FLCI_BC, 1)
            res = findOptimalFLCI(bh, sg, FLCI_BC[r, 1], np, nq)
            @test isapprox(res.FLCI[1], FLCI_BC[r, 2]; atol=1e-3)
            @test isapprox(res.FLCI[2], FLCI_BC[r, 3]; atol=1e-3)
            @test isapprox(res.optimalHalfLength, FLCI_BC[r, 4]; atol=1e-3)
        end
        @test isapprox(findOptimalFLCI(bh, sg, 0.1, np, nq).optimalVec, FLCI_OPTVEC_M01; atol=1e-3)
    end

    @testset "computeConditionalCS_DeltaSD ARP (exact)" begin
        cs = computeConditionalCS_DeltaSD(bh, sg, np, nq; M=0.1, hybrid_flag="ARP", gridPoints=200)
        lb, ub = confidence_interval(cs)
        @test isapprox(lb, CCS_SD_ARP[1]; atol=1e-6)
        @test isapprox(ub, CCS_SD_ARP[2]; atol=1e-6)
    end

    @testset "computeConditionalCS_DeltaSD FLCI hybrid (atol 2e-3)" begin
        cs = computeConditionalCS_DeltaSD(bh, sg, np, nq; M=0.1, hybrid_flag="FLCI", gridPoints=200)
        lb, ub = confidence_interval(cs)
        @test isapprox(lb, CCS_SD_FLCI[1]; atol=2e-3)
        @test isapprox(ub, CCS_SD_FLCI[2]; atol=2e-3)
    end

    @testset "createSensitivityResults Conditional (atol 2e-3)" begin
        res = createSensitivityResults(bh, sg, np, nq; method="Conditional", Mvec=[0.0, 0.1, 0.2])
        for r in axes(SENS_COND, 1)
            @test isapprox(res.lb[r], SENS_COND[r, 2]; atol=2e-3)
            @test isapprox(res.ub[r], SENS_COND[r, 3]; atol=2e-3)
        end
    end

    @testset "createSensitivityResults FLCI (atol 1e-3)" begin
        res = createSensitivityResults(bh, sg, np, nq; method="FLCI", Mvec=[0.0, 0.1, 0.2])
        for r in axes(SENS_FLCI, 1)
            @test isapprox(res.lb[r], SENS_FLCI[r, 2]; atol=1e-3)
            @test isapprox(res.ub[r], SENS_FLCI[r, 3]; atol=1e-3)
        end
    end

    @testset "DeltaSD M bounds" begin
        @test isapprox(DeltaSD_upperBound_Mpre(bh, sg, np), M_UPPERBOUND; atol=1e-7)
        @test isinf(M_LOWERBOUND) && isinf(DeltaSD_lowerBound_Mpre(bh, sg, np; gridPoints=200))
    end

    @testset "Tables.jl interop + pretty show" begin
        cs = computeConditionalCS_DeltaSD(bh, sg, np, nq; M=0.1, hybrid_flag="ARP", gridPoints=50)
        @test Tables.istable(typeof(cs))
        @test Tables.columnnames(Tables.columns(cs)) == (:grid, :accept)
        res = constructOriginalCS(bh, sg, np, nq)
        @test Tables.istable(typeof(res))
        @test :lb in propertynames(res)
        @test occursin("SensitivityResults", sprint(show, MIME("text/plain"), res))
        @test occursin("ConditionalCS", sprint(show, MIME("text/plain"), cs))
    end
end
