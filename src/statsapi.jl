# Integration with the StatsAPI / StatsBase ecosystem, in particular
# StagDiDModels.jl event-study models. Dispatch is generic on any
# StatsAPI.RegressionModel whose dynamic coefficients are named "τ::<int>".

"""
    eventstudy_inputs(model; ref_p=-1, pattern=r"τ::(-?\\d+)")

Extract `(betahat, sigma, taus, numPrePeriods, numPostPeriods)` from a fitted
event-study model (anything implementing `StatsAPI.coef`/`vcov`/`coefnames`
with `τ::<event-time>` coefficient names, e.g. a StagDiDModels.jl dynamic model).
Coefficients are ordered ascending by event time. The reference period `ref_p`
(default `-1`) is assumed already absent from the coefficients; pre/post are split
as event times `< ref_p` / `> ref_p`.
"""
function eventstudy_inputs(model; ref_p::Integer=-1, pattern::Regex=r"τ::(-?\d+)")
    names = StatsAPI.coefnames(model)
    idx = Int[]; taus = Int[]
    for (i, nm) in enumerate(names)
        m = match(pattern, string(nm))
        m === nothing && continue
        push!(idx, i); push!(taus, parse(Int, m.captures[1]))
    end
    isempty(idx) && error("no event-study coefficients (pattern $(pattern)) found in model")
    perm = sortperm(taus)
    idx = idx[perm]; taus = taus[perm]
    betahat = collect(float.(StatsAPI.coef(model)[idx]))
    sigma = Matrix{Float64}(StatsAPI.vcov(model)[idx, idx])
    numPrePeriods = count(<(ref_p), taus)
    numPostPeriods = count(>(ref_p), taus)
    return (betahat=betahat, sigma=sigma, taus=taus,
            numPrePeriods=numPrePeriods, numPostPeriods=numPostPeriods)
end

"""
    honest_did(model; e=0, type="smoothness", gridPoints=100, ref_p=-1, kwargs...)

Sensitivity analysis (Rambachan & Roth) for the event-time-`e` effect of a fitted
event-study `model`. `type` is `"smoothness"` (Delta^SD) or `"relative_magnitude"`
(Delta^RM). Extra `kwargs` are forwarded to `createSensitivityResults` /
`createSensitivityResults_relativeMagnitudes`. Returns
`(; robust_ci, orig_ci, type)`.
"""
function honest_did(model; e::Integer=0, type::AbstractString="smoothness",
                    gridPoints::Integer=100, ref_p::Integer=-1, kwargs...)
    inp = eventstudy_inputs(model; ref_p=ref_p)

    # Consecutive-time-period check (mirrors R honest_did.AGGTEobj).
    full = sort(vcat(inp.taus, ref_p))
    all(diff(full) .== 1) || error("honest_did expects consecutive event-time periods " *
                                   "(with $(ref_p) as the reference); please re-code your event study.")
    (inp.numPostPeriods > 0) || error("not enough post-periods (reference period is $(ref_p))")
    (inp.numPrePeriods > 0) || error("not enough pre-periods (reference period is $(ref_p))")

    baseVec = basisVector(e + 1, inp.numPostPeriods)
    orig_ci = constructOriginalCS(inp.betahat, inp.sigma, inp.numPrePeriods, inp.numPostPeriods; l_vec=baseVec)

    if type == "relative_magnitude"
        robust_ci = createSensitivityResults_relativeMagnitudes(inp.betahat, inp.sigma,
            inp.numPrePeriods, inp.numPostPeriods; l_vec=baseVec, gridPoints=gridPoints, kwargs...)
    elseif type == "smoothness"
        robust_ci = createSensitivityResults(inp.betahat, inp.sigma,
            inp.numPrePeriods, inp.numPostPeriods; l_vec=baseVec, kwargs...)
    else
        error("type must be 'smoothness' or 'relative_magnitude'")
    end
    return (robust_ci=robust_ci, orig_ci=orig_ci, type=type)
end
