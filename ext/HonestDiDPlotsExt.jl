module HonestDiDPlotsExt

# Plots.jl rendering for the HonestDiD plot functions. Data preparation lives in
# HonestDiD core (backend-agnostic); this module only draws.

using HonestDiD
using Plots

const HD = HonestDiD

function _sensitivity(robustResults, originalResults, xcol, xlabel; rescaleFactor, maxX, add_xAxis)
    d = HD._sensitivity_plot_data(robustResults, originalResults, xcol; rescaleFactor=rescaleFactor, maxX=maxX)
    plt = plot(; xlabel=xlabel, ylabel="", legend=:bottom)
    for m in unique(d.method)
        idx = d.method .== m
        color = m == "Original" ? :red : HD.ROBUST_COLOR
        mid = (d.lb[idx] .+ d.ub[idx]) ./ 2
        scatter!(plt, d.x[idx], mid; yerror=(mid .- d.lb[idx], d.ub[idx] .- mid),
                 label=m, color=color, markerstrokecolor=color, markersize=3)
    end
    add_xAxis && hline!(plt, [0]; color=:black, linestyle=:dash, label="")
    return plt
end

HD.createSensitivityPlot(robustResults, originalResults; rescaleFactor=1, maxM=Inf, add_xAxis=true) =
    _sensitivity(robustResults, originalResults, :M, "M"; rescaleFactor=rescaleFactor, maxX=maxM, add_xAxis=add_xAxis)

HD.createSensitivityPlot_relativeMagnitudes(robustResults, originalResults; rescaleFactor=1, maxMbar=Inf, add_xAxis=true) =
    _sensitivity(robustResults, originalResults, :Mbar, "Mbar"; rescaleFactor=rescaleFactor, maxX=maxMbar, add_xAxis=add_xAxis)

function HD.createEventStudyPlot(betahat; stdErrors=nothing, sigma=nothing,
        numPrePeriods, numPostPeriods, alpha=0.05, timeVec, referencePeriod, useRelativeEventTime=false)
    d = HD._eventstudy_plot_data(betahat; stdErrors=stdErrors, sigma=sigma, numPrePeriods=numPrePeriods,
        numPostPeriods=numPostPeriods, alpha=alpha, timeVec=timeVec, referencePeriod=referencePeriod,
        useRelativeEventTime=useRelativeEventTime)
    fin = .!isnan.(d.se)
    plt = plot(; xlabel="Event time", ylabel="", legend=false)
    scatter!(plt, d.t, d.beta; color=:red, markersize=4)
    scatter!(plt, d.t[fin], d.beta[fin]; yerror=d.z .* d.se[fin],
             color=HD.ROBUST_COLOR, markerstrokecolor=HD.ROBUST_COLOR)
    return plt
end

end # module
