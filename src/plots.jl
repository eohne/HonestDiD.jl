# Plotting entry points. The real methods live in `ext/HonestDiDPlotsExt.jl`
# and are loaded automatically when `Plots` is available (`using Plots`).

"""
    createSensitivityPlot(robustResults, originalResults; rescaleFactor=1, maxM=Inf, add_xAxis=true)

Plot the robust confidence intervals from [`createSensitivityResults`](@ref)
against the smoothness bound `M`, with the original (non-robust) interval from
[`constructOriginalCS`](@ref) drawn at the far left. The robust CIs are blue, the
original is red, and a dashed line is drawn at 0 (unless `add_xAxis = false`).

Requires a plotting backend - `using Plots` returns a `Plots.Plot`, `using Makie`
(e.g. CairoMakie) returns a Makie `Figure`. `rescaleFactor` multiplies all values;
`maxM` truncates the x-axis.
"""
function createSensitivityPlot end

"""
    createSensitivityPlot_relativeMagnitudes(robustResults, originalResults; rescaleFactor=1, maxMbar=Inf, add_xAxis=true)

As [`createSensitivityPlot`](@ref) but for relative-magnitude results
(x-axis `M̄`, from [`createSensitivityResults_relativeMagnitudes`](@ref)).
Requires a plotting backend (`using Plots` or `using Makie`).
"""
function createSensitivityPlot_relativeMagnitudes end

"""
    createEventStudyPlot(betahat; sigma=nothing, stdErrors=nothing, numPrePeriods, numPostPeriods, alpha=0.05, timeVec, referencePeriod, useRelativeEventTime=false)

Event-study coefficient plot: the point estimates with pointwise `1 - alpha`
confidence intervals, including the reference period pinned at 0.

Supply uncertainty via either `sigma` (covariance matrix) or `stdErrors`.
`timeVec` is the vector of event times (excluding the reference), `referencePeriod`
the omitted period; set `useRelativeEventTime = true` to re-centre time at the
reference. Requires a plotting backend (`using Plots` or `using Makie`).
"""
function createEventStudyPlot end

_needs_plots(name) = error("$(name) requires a plotting backend; run `using Plots` or `using Makie` (e.g. CairoMakie) first.")
createSensitivityPlot(args...; kwargs...) = _needs_plots("createSensitivityPlot")
createSensitivityPlot_relativeMagnitudes(args...; kwargs...) = _needs_plots("createSensitivityPlot_relativeMagnitudes")
createEventStudyPlot(args...; kwargs...) = _needs_plots("createEventStudyPlot")

const ROBUST_COLOR = "#01a2d9"

# Backend-agnostic data preparation shared by the Plots and Makie extensions.
function _sensitivity_plot_data(robustResults, originalResults, xcol::Symbol; rescaleFactor=1, maxX=Inf)
    robx = collect(float.(getproperty(robustResults, xcol)))
    xgap = minimum(diff(sort(robx)))
    xmin = minimum(robx)
    n_orig = length(originalResults.lb)
    x = vcat(fill(xmin - xgap, n_orig), robx) .* rescaleFactor
    lb = vcat(collect(float.(originalResults.lb)), collect(float.(robustResults.lb))) .* rescaleFactor
    ub = vcat(collect(float.(originalResults.ub)), collect(float.(robustResults.ub))) .* rescaleFactor
    method = vcat(collect(string.(originalResults.method)), collect(string.(robustResults.method)))
    keep = x .<= maxX
    return (x=x[keep], lb=lb[keep], ub=ub[keep], method=method[keep])
end

function _eventstudy_plot_data(betahat; stdErrors=nothing, sigma=nothing, numPrePeriods, numPostPeriods,
                               alpha=0.05, timeVec, referencePeriod, useRelativeEventTime=false)
    bh = collect(float.(betahat isa AbstractMatrix ? vec(betahat) : betahat))
    if stdErrors === nothing && sigma === nothing
        error("User must specify either vector of standard errors or vcv matrix!")
    elseif stdErrors === nothing
        stdErrors = sqrt.(diag(sigma))
    end
    tv = collect(float.(timeVec)); ref = float(referencePeriod)
    if useRelativeEventTime
        tv = tv .- ref; ref = 0.0
    end
    pre = 1:numPrePeriods
    post = (numPrePeriods+1):(numPrePeriods+numPostPeriods)
    t = vcat(tv[pre], ref, tv[post])
    beta = vcat(bh[pre], 0.0, bh[post])
    se = vcat(collect(float.(stdErrors[pre])), NaN, collect(float.(stdErrors[post])))
    z = quantile(Normal(), 1 - alpha / 2)
    return (t=t, beta=beta, se=se, z=z)
end
