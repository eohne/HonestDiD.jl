# StatsAPI / StagDiDModels integration + least-favorable (simulation) path checks.
using HonestDiD, Test
import StatsAPI
isdefined(Main, :BC_BETAHAT) || include(joinpath(@__DIR__, "reference.jl"))

# A minimal event-study RegressionModel mimicking a StagDiDModels.jl dynamic fit
# (coefficient names "τ::<event-time>", reference period dropped).
struct MockES <: StatsAPI.RegressionModel
    beta::Vector{Float64}
    V::Matrix{Float64}
    names::Vector{String}
end
StatsAPI.coef(m::MockES) = m.beta
StatsAPI.vcov(m::MockES) = m.V
StatsAPI.coefnames(m::MockES) = m.names

@testset "StatsAPI integration (eventstudy_inputs / honest_did)" begin
    # Event times -3,-2,(ref -1 dropped),0,1,2 ⇒ numPre=2, numPost=3. Scrambled
    # order + a non-event-study coef to exercise parsing/sorting/filtering.
    order = [3, 1, 5, 2, 4]
    taus = [-3, -2, 0, 1, 2]
    betas = [0.2, -0.1, 0.5, 0.4, 0.35]
    M = [2.0 0.3 0.1 0.0 0.05; 0.3 1.5 0.2 0.1 0.0; 0.1 0.2 1.8 0.25 0.1;
         0.0 0.1 0.25 1.6 0.2; 0.05 0.0 0.1 0.2 1.4]
    names = vcat(["(Intercept)"], ["τ::$(taus[i])" for i in order])
    beta = vcat([99.0], betas[order])
    V = zeros(6, 6); V[2:6, 2:6] = M[order, order]; V[1, 1] = 7.0

    m = MockES(beta, V, names)
    inp = eventstudy_inputs(m)
    @test inp.numPrePeriods == 2
    @test inp.numPostPeriods == 3
    @test inp.taus == taus
    @test inp.betahat ≈ betas
    @test inp.sigma ≈ M

    bv = HonestDiD.basisVector(1, 3)
    hd = honest_did(m; e=0, type="smoothness", method="FLCI", Mvec=[0.0, 0.05])
    direct = createSensitivityResults(betas, M, 2, 3; l_vec=bv, method="FLCI", Mvec=[0.0, 0.05])
    @test hd.orig_ci.lb ≈ constructOriginalCS(betas, M, 2, 3; l_vec=bv).lb
    @test hd.robust_ci.lb ≈ direct.lb
    @test hd.robust_ci.ub ≈ direct.ub
end

@testset "Least-favorable (simulation) path vs R (loose)" begin
    bh, sg, np, nq = BC_BETAHAT, BC_SIGMA, BC_NPRE, BC_NPOST
    # RNG differs from R, so LF critical values (hence CIs) match only approximately.
    for (cs, ref) in (
        (computeConditionalCS_DeltaSD(bh, sg, np, nq; M=0.1, hybrid_flag="LF", gridPoints=200, seed=0), LF_SD),
        (computeConditionalCS_DeltaRM(bh, sg, np, nq; Mbar=1, hybrid_flag="LF", gridPoints=200, seed=0), LF_RM),
    )
        lb, ub = confidence_interval(cs)
        @test isapprox(lb, ref[1]; atol=0.02)
        @test isapprox(ub, ref[2]; atol=0.02)
    end
end
