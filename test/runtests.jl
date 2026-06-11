using HonestDiD
using Test

# Reference values (extracted from the R HonestDiD package on BCdata) are
# embedded in reference.jl — no CSV fixtures or R installation required.
include(joinpath(@__DIR__, "reference.jl"))

@testset "HonestDiD.jl" begin
    include("validate_builders.jl")     # delta matrix builders + numPre==1 bug fix
    include("validate_e2e.jl")          # LP/duals, FLCI, ARP grids, sensitivity, M bounds, Tables/show
    include("validate_variants.jl")     # all 9 computeConditionalCS variants (ARP, exact)
    include("validate_integration.jl")  # StatsAPI/honest_did + least-favorable (loose)
end
