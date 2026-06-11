# Upper/lower bounds on the smoothness parameter M from the pre-period,
# ported from R/ublbM_functions.R.

# .testInIdentifiedSet_Max : APR test of E[A Y - 1*M] <= 0 (no max(0,.) on CV).
function testInIdentifiedSet_Max(M, y, sigma, A, alpha, d)
    yv = _asvec(y); dv = _asvec(d)
    d_mod = dv .* M
    sigmaTilde = sqrt.(diag(A * sigma * A'))
    Atilde = Diagonal(1.0 ./ sigmaTilde) * A
    dtilde = (1.0 ./ sigmaTilde) .* d_mod

    normalizedMoments = Atilde * yv - dtilde
    maxLocation = argmax(normalizedMoments)
    maxMoment = normalizedMoments[maxLocation]

    m = size(Atilde, 1)
    T_B = selectionMat(maxLocation, m; select="rows")
    iota = ones(m, 1)
    gamma = vec(T_B * Atilde)
    Abar = Atilde - iota * (T_B * Atilde)
    dbar = (Matrix{Float64}(I, m, m) - iota * T_B) * dtilde

    sigmabar = sqrt(_scalar(gamma' * sigma * gamma))
    VLo, VUp = VLoVUpFN(gamma, sigma, Abar, dbar,
                        (Matrix{Float64}(I, length(yv), length(yv)) -
                         ((sigma * gamma) ./ _scalar(gamma' * sigma * gamma)) * gamma') * yv)
    mu = _scalar(T_B * dtilde)
    criticalVal = norminvp_generalized(1 - alpha, VLo, VUp; mu=mu, sd=sigmabar)
    return (maxMoment + mu) > criticalVal
end

# .create_A_and_D_SD_prePeriods
function create_A_and_D_SD_prePeriods(numPrePeriods)
    numPrePeriods < 2 && error("Can't estimate M in pre-period with < 2 pre-period coeffs")
    Atilde = zeros(numPrePeriods - 1, numPrePeriods)
    Atilde[numPrePeriods-1, (numPrePeriods-1):numPrePeriods] = [1.0, -2.0]
    for r in 1:(numPrePeriods-2)
        Atilde[r, r:(r+2)] = [1.0, -2.0, 1.0]
    end
    A_pre = vcat(Atilde, -Atilde)
    d = ones(size(A_pre, 1))
    return (A=A_pre, d=d)
end

# .estimate_lowerBound_M_conditionalTest
function estimate_lowerBound_M_conditionalTest(prePeriodCoef, prePeriodCovar, grid_ub; alpha=0.05, gridPoints)
    numPre = length(_asvec(prePeriodCoef))
    Ad = create_A_and_D_SD_prePeriods(numPre)
    mGrid = collect(range(0.0, grid_ub; length=gridPoints))
    accept = [(testInIdentifiedSet_Max(maxVal, prePeriodCoef, prePeriodCovar, Ad.A, alpha, Ad.d) ? 0.0 : 1.0)
              for maxVal in mGrid]
    if sum(accept) == 0
        @warn "ARP conditional test rejects all values of M provided. User should increase upper bound of grid."
        return Inf
    else
        return minimum(mGrid[accept .== 1])
    end
end

"""
    DeltaSD_upperBound_Mpre(betahat, sigma, numPrePeriods; alpha=0.05)

A data-driven `1 - alpha` upper bound on the smoothness parameter `M`, formed from
the largest second difference of the observed *pre-period* coefficients plus a
one-sided normal critical value. Useful for choosing a default range of `M` in a
sensitivity analysis (this is what [`createSensitivityResults`](@ref) uses when
`Mvec` is not supplied). Requires `numPrePeriods > 1`.
"""
function DeltaSD_upperBound_Mpre(betahat, sigma, numPrePeriods; alpha=0.05)
    numPrePeriods > 1 || error("numPrePeriods must be > 1")
    bh = _asvec(betahat)
    prePeriodCoef = bh[1:numPrePeriods]
    prePeriodSigma = sigma[1:numPrePeriods, 1:numPrePeriods]
    A_SD = create_A_SD(numPrePeriods, 0)
    prePeriodCoefDiffs = A_SD * prePeriodCoef
    prePeriodSigmaDiffs = A_SD * prePeriodSigma * A_SD'
    seDiffs = sqrt.(diag(prePeriodSigmaDiffs))
    upperBoundVec = prePeriodCoefDiffs .+ quantile(Normal(), 1 - alpha) .* seDiffs
    return maximum(upperBoundVec)
end

"""
    DeltaSD_lowerBound_Mpre(betahat, sigma, numPrePeriods; alpha=0.05, grid_ub=nothing, gridPoints=1000)

A one-sided `1 - alpha` *lower* bound on the smoothness parameter `M`, obtained by
inverting the Andrews–Roth–Pakes (2019) conditional test on the maximal second
difference of the observed pre-period coefficients. Returns `Inf` (with a warning)
if every `M` on the grid is rejected - increase `grid_ub` in that case. Requires
`numPrePeriods > 1`.
"""
function DeltaSD_lowerBound_Mpre(betahat, sigma, numPrePeriods; alpha=0.05, grid_ub=nothing, gridPoints=1000)
    numPrePeriods > 1 || error("numPrePeriods must be > 1")
    bh = _asvec(betahat)
    prePeriodCoef = bh[1:numPrePeriods]
    prePeriodSigma = sigma[1:numPrePeriods, 1:numPrePeriods]
    gub = grid_ub === nothing ? 3 * maximum(sqrt.(diag(prePeriodSigma))) : grid_ub
    return estimate_lowerBound_M_conditionalTest(prePeriodCoef, prePeriodSigma, gub;
                                                 alpha=alpha, gridPoints=gridPoints)
end
