# All nine computeConditionalCS_Delta* variants vs R (deterministic ARP test,
# gridPoints=100 to keep runtime modest). Reference values in reference.jl.
using HonestDiD, Test
isdefined(Main, :BC_BETAHAT) || include(joinpath(@__DIR__, "reference.jl"))

bh, sg, np, nq = BC_BETAHAT, BC_SIGMA, BC_NPRE, BC_NPOST

function check(cs, ref; atol=1e-6)
    lb, ub = confidence_interval(cs)
    @test isapprox(lb, ref[1]; atol=atol)
    @test isapprox(ub, ref[2]; atol=atol)
end

@testset "computeConditionalCS variants vs R (ARP, gp100)" begin
    check(computeConditionalCS_DeltaSD(bh, sg, np, nq; l_vec=fill(1 / nq, nq), M=0.1, hybrid_flag="ARP", gridPoints=100), CCS_SD_AVG_ARP)
    check(computeConditionalCS_DeltaSDB(bh, sg, np, nq; M=0.1, hybrid_flag="ARP", biasDirection="positive", gridPoints=100), CCS_SDB)
    check(computeConditionalCS_DeltaSDM(bh, sg, np, nq; M=0.1, hybrid_flag="ARP", monotonicityDirection="increasing", gridPoints=100), CCS_SDM)
    check(computeConditionalCS_DeltaRMB(bh, sg, np, nq; Mbar=1, hybrid_flag="ARP", biasDirection="positive", gridPoints=100), CCS_RMB)
    check(computeConditionalCS_DeltaRMM(bh, sg, np, nq; Mbar=1, hybrid_flag="ARP", monotonicityDirection="increasing", gridPoints=100), CCS_RMM)
    check(computeConditionalCS_DeltaSDRM(bh, sg, np, nq; Mbar=1, hybrid_flag="ARP", gridPoints=100), CCS_SDRM)
    check(computeConditionalCS_DeltaSDRMB(bh, sg, np, nq; Mbar=1, hybrid_flag="ARP", biasDirection="negative", gridPoints=100), CCS_SDRMB)
    check(computeConditionalCS_DeltaSDRMM(bh, sg, np, nq; Mbar=1, hybrid_flag="ARP", monotonicityDirection="decreasing", gridPoints=100), CCS_SDRMM)
end
