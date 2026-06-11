# ARP test without nuisance parameters, ported from R/arp-nonuisance.R.
# hybrid_list is a Dict{Symbol,Any} mirroring the R list (keys :hybrid_kappa,
# :lf_cv, :flci_l, :flci_halflength, :vbar, :dbar).

# .testInIdentifiedSet
function testInIdentifiedSet(y, sigma, A, d; Abar_additional=nothing, dbar_additional=nothing, alpha)
    yv = _asvec(y); dv = _asvec(d)
    sigmaTilde = sqrt.(diag(A * sigma * A'))
    Atilde = Diagonal(1.0 ./ sigmaTilde) * A
    dtilde = (1.0 ./ sigmaTilde) .* dv

    normalizedMoments = Atilde * yv - dtilde
    maxLocation = argmax(normalizedMoments)
    maxMoment = normalizedMoments[maxLocation]

    m = size(Atilde, 1)
    T_B = selectionMat(maxLocation, m; select="rows")          # 1×m
    iota = ones(m, 1)
    gamma = vec(T_B * Atilde)                                   # length(y)
    Abar = Atilde - iota * (T_B * Atilde)
    dbar = (Matrix{Float64}(I, m, m) - iota * T_B) * dtilde

    if Abar_additional !== nothing
        Abar = vcat(Abar, Abar_additional)
        dbar = vcat(dbar, _asvec(dbar_additional))
    end

    sigmabar = sqrt(_scalar(gamma' * sigma * gamma))
    c = (sigma * gamma) ./ _scalar(gamma' * sigma * gamma)
    z = (Matrix{Float64}(I, length(yv), length(yv)) - c * gamma') * yv
    VLo, VUp = VLoVUpFN(gamma, sigma, Abar, dbar, z)

    mu = _scalar(T_B * dtilde)
    criticalVal = max(0.0, norminvp_generalized(1 - alpha, VLo, VUp; mu=mu, sd=sigmabar))
    return (maxMoment + mu) > criticalVal
end

# .testInIdentifiedSet_FLCI_Hybrid
function testInIdentifiedSet_FLCI_Hybrid(y, sigma, A, d, alpha, hybrid_list)
    yv = _asvec(y)
    flci_l = _asvec(hybrid_list[:flci_l])
    A_firststage = vcat(reshape(flci_l, 1, :), reshape(-flci_l, 1, :))
    d_firststage = [hybrid_list[:flci_halflength], hybrid_list[:flci_halflength]]
    if maximum(A_firststage * yv - d_firststage) > 0
        return true
    else
        alphatilde = (alpha - hybrid_list[:hybrid_kappa]) / (1 - hybrid_list[:hybrid_kappa])
        return testInIdentifiedSet(yv, sigma, A, d;
                                   Abar_additional=A_firststage, dbar_additional=d_firststage, alpha=alphatilde)
    end
end

# .testInIdentifiedSet_LF_Hybrid
function testInIdentifiedSet_LF_Hybrid(y, sigma, A, d, alpha, hybrid_list)
    yv = _asvec(y); dv = _asvec(d)
    sigmaTilde = sqrt.(diag(A * sigma * A'))
    Atilde = Diagonal(1.0 ./ sigmaTilde) * A
    dtilde = (1.0 ./ sigmaTilde) .* dv

    normalizedMoments = Atilde * yv - dtilde
    maxLocation = argmax(normalizedMoments)
    maxMoment = normalizedMoments[maxLocation]

    if maxMoment > hybrid_list[:lf_cv]
        return true
    end

    m = size(Atilde, 1)
    T_B = selectionMat(maxLocation, m; select="rows")
    iota = ones(m, 1)
    gamma = vec(T_B * Atilde)
    Abar = Atilde - iota * (T_B * Atilde)
    dbar = (Matrix{Float64}(I, m, m) - iota * T_B) * dtilde

    sigmabar = sqrt(_scalar(gamma' * sigma * gamma))
    c = (sigma * gamma) ./ _scalar(gamma' * sigma * gamma)
    z = (Matrix{Float64}(I, length(yv), length(yv)) - c * gamma') * yv
    VLo, VUp = VLoVUpFN(gamma, sigma, Abar, dbar, z)

    alphatilde = (alpha - hybrid_list[:hybrid_kappa]) / (1 - hybrid_list[:hybrid_kappa])
    mu = _scalar(T_B * dtilde)
    criticalVal = max(0.0, norminvp_generalized(1 - alphatilde, VLo, VUp; mu=mu, sd=sigmabar))
    return (maxMoment + mu) > criticalVal
end

# .testOverThetaGrid
function testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha;
                           testFn=testInIdentifiedSet_ARP, hybrid_list=nothing)
    bh = _asvec(betahat)
    e = basisVector(numPrePeriods + 1, length(bh))
    accept = Float64[]
    for theta in thetaGrid
        y = bh - e .* theta
        reject = testFn(y, sigma, A, d, alpha, hybrid_list)
        push!(accept, reject ? 0.0 : 1.0)
    end
    return accept
end

# Uniform-signature wrapper for the plain ARP test.
testInIdentifiedSet_ARP(y, sigma, A, d, alpha, hybrid_list) =
    testInIdentifiedSet(y, sigma, A, d; alpha=alpha)

# .APR_computeCI_NoNuis
function APR_computeCI_NoNuis(betahat, sigma, A, d, numPrePeriods, numPostPeriods, l_vec,
                              alpha, returnLength, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    thetaGrid = collect(range(grid_lb, grid_ub; length=gridPoints))
    if hybrid_flag == "ARP"
        accept = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha;
                                   testFn=testInIdentifiedSet_ARP)
    elseif hybrid_flag == "FLCI"
        accept = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha;
                                   testFn=testInIdentifiedSet_FLCI_Hybrid, hybrid_list=hybrid_list)
    elseif hybrid_flag == "LF"
        accept = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha;
                                   testFn=testInIdentifiedSet_LF_Hybrid, hybrid_list=hybrid_list)
    else
        error("hybrid_flag must equal 'ARP' or 'FLCI' or 'LF'")
    end

    if accept[1] == 1 || accept[end] == 1
        @warn "CI is open at one of the endpoints; CI length may not be accurate"
    end

    if returnLength
        gridLength = 0.5 .* (vcat(0.0, diff(thetaGrid)) .+ vcat(diff(thetaGrid), 0.0))
        return sum(accept .* gridLength)
    else
        return ConditionalCS(thetaGrid, accept)
    end
end
