# Exercises code paths the main BC suite (4 post-periods) never touches: the
# numPostPeriods == 1 no-nuisance path, the sensitivity orchestrators, the RM/SDRM
# identified sets, and the backend-agnostic plotting helpers. Kept on small inputs.
using HonestDiD, Test
isdefined(Main, :BC_BETAHAT) || include(joinpath(@__DIR__, "reference.jl"))

const H = HonestDiD

@testset "numPostPeriods == 1 (no-nuisance path)" begin
    bh = BC_BETAHAT[1:3]; sg = BC_SIGMA[1:3, 1:3]   # 2 pre, 1 post

    arp = computeConditionalCS_DeltaSD(bh, sg, 2, 1; M=0.1, hybrid_flag="ARP", gridPoints=200)
    @test all(isapprox.(confidence_interval(arp), NP1_SD_ARP; atol=1e-6))

    flci = computeConditionalCS_DeltaSD(bh, sg, 2, 1; M=0.1, hybrid_flag="FLCI", gridPoints=200)
    @test all(isapprox.(confidence_interval(flci), NP1_SD_FLCI; atol=2e-3))

    lf = computeConditionalCS_DeltaSD(bh, sg, 2, 1; M=0.1, hybrid_flag="LF", gridPoints=200, seed=0)
    @test all(isapprox.(confidence_interval(lf), NP1_SD_LF; atol=0.02))

    # the SDB / SDM / RM single-post-period branches (sanity)
    for cs in (
        computeConditionalCS_DeltaSDB(bh, sg, 2, 1; M=0.1, hybrid_flag="LF", biasDirection="positive", gridPoints=100),
        computeConditionalCS_DeltaSDM(bh, sg, 2, 1; M=0.1, hybrid_flag="ARP", monotonicityDirection="increasing", gridPoints=100),
        computeConditionalCS_DeltaRM(bh, sg, 2, 1; Mbar=1, hybrid_flag="LF", gridPoints=100),
    )
        lb, ub = confidence_interval(cs)
        @test lb <= ub
    end
end

@testset "sensitivity orchestrators (smoke)" begin
    bh = BC_BETAHAT[1:5]; sg = BC_SIGMA[1:5, 1:5]; np, nq = 3, 2

    # smoothness with sign / monotonicity (FLCI keeps it fast; it warns, as in R)
    r1 = createSensitivityResults(bh, sg, np, nq; method="FLCI", Mvec=[0.0, 0.1], biasDirection="negative")
    @test all(==("DeltaSDNB"), r1.Delta) && all(r1.lb .<= r1.ub)
    r2 = createSensitivityResults(bh, sg, np, nq; method="FLCI", Mvec=[0.0, 0.1], monotonicityDirection="increasing")
    @test all(==("DeltaSDI"), r2.Delta)

    # relative magnitudes: each bound / restriction branch (Conditional avoids the
    # 1000-draw least-favorable simulation, so this stays fast)
    rmkw = (; method="Conditional", Mbarvec=[0.5], gridPoints=40)
    rm   = createSensitivityResults_relativeMagnitudes(bh, sg, np, nq; rmkw...)
    rmm  = createSensitivityResults_relativeMagnitudes(bh, sg, np, nq; rmkw..., monotonicityDirection="decreasing")
    rmb  = createSensitivityResults_relativeMagnitudes(bh, sg, np, nq; rmkw..., biasDirection="positive")
    sdrm = createSensitivityResults_relativeMagnitudes(bh, sg, np, nq; rmkw..., bound="deviation from linear trend")
    @test all(==("DeltaRM"), rm.Delta)
    @test all(==("DeltaRMD"), rmm.Delta)
    @test all(==("DeltaRMPB"), rmb.Delta)
    @test all(==("DeltaSDRM"), sdrm.Delta)
end

@testset "identified sets (RM / SDRM)" begin
    bh = BC_BETAHAT[1:5]; lv = [1.0, 0.0]
    @test (x -> x[1] <= x[2])(H.compute_IDset_DeltaRM(1.0, bh, lv, 3, 2))
    @test (x -> x[1] <= x[2])(H.compute_IDset_DeltaSDRM(1.0, bh, lv, 3, 2))
end

@testset "plotting helpers (no backend)" begin
    res = H.SensitivityResults((lb=[0.1, 0.05], ub=[0.3, 0.35], method=["FLCI", "FLCI"],
                                Delta=["DeltaSD", "DeltaSD"], M=[0.0, 0.1]))
    orig = H.SensitivityResults((lb=[0.15], ub=[0.25], method=["Original"], Delta=Union{String,Missing}[missing]))

    d = H._sensitivity_plot_data(res, orig, :M)
    @test length(d.x) == 3 && "Original" in d.method && length(d.lb) == length(d.ub)

    ev = H._eventstudy_plot_data([0.01, 0.03, 0.07]; sigma=[0.01 0.0 0.0; 0.0 0.01 0.0; 0.0 0.0 0.01],
                                 numPrePeriods=2, numPostPeriods=1, timeVec=[-2, -1, 1], referencePeriod=0)
    @test length(ev.t) == 4 && ev.beta[3] == 0.0   # reference pinned to 0

    # the stub functions error helpfully when no backend is loaded
    @test_throws ErrorException createSensitivityPlot(res, orig)
    @test_throws ErrorException createEventStudyPlot([0.1]; numPrePeriods=0, numPostPeriods=1, timeVec=[1], referencePeriod=0)
end
