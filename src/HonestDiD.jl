module HonestDiD

# Julia port of the R package HonestDiD (Rambachan & Roth 2023, "A More Credible
# Approach to Parallel Trends"). Robust inference and sensitivity analysis for
# difference-in-differences / event-study designs.

using LinearAlgebra
using Random
using Statistics
using Printf
using Distributions
using JuMP
import HiGHS
import ECOS
import StatsAPI
import Tables

const MOI = JuMP.MOI

include("utilities.jl")
include("results.jl")
include("deltas.jl")
include("solvers.jl")
include("identified_sets.jl")
include("flci.jl")
include("arp_nonuisance.jl")
include("arp_nuisance.jl")
include("bounds.jl")
include("conditional_cs.jl")
include("sensitivity.jl")
include("statsapi.jl")
include("plots.jl")

# Public API
export basisVector
export constructOriginalCS
export createSensitivityResults, createSensitivityResults_relativeMagnitudes
export findOptimalFLCI
export DeltaSD_upperBound_Mpre, DeltaSD_lowerBound_Mpre
export computeConditionalCS_DeltaSD, computeConditionalCS_DeltaSDB, computeConditionalCS_DeltaSDM
export computeConditionalCS_DeltaRM, computeConditionalCS_DeltaRMB, computeConditionalCS_DeltaRMM
export computeConditionalCS_DeltaSDRM, computeConditionalCS_DeltaSDRMB, computeConditionalCS_DeltaSDRMM
export eventstudy_inputs, honest_did
export createSensitivityPlot, createSensitivityPlot_relativeMagnitudes, createEventStudyPlot
export ConditionalCS, SensitivityResults, confidence_interval

end # module
