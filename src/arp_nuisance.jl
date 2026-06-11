# ARP test with nuisance parameters, ported from R/arp-nuisance.R.
# HiGHS LPs (with duals) replace lpSolveAPI/Rglpk; ECOS handles the bias programs.
# hybrid_list is a Dict{Symbol,Any}.
#
# Speed note: across the theta grid (X_T and sigma fixed) only the RHS of the
# eta-LP and the objective of the dual max-program change with theta. So
# LPWorkspace builds both HiGHS models once and just updates the RHS/objective
# per solve, instead of building a fresh JuMP model at every grid point. That
# rebuild is what made the naive version slower than R.

mutable struct LPWorkspace
    sdVec::Vector{Float64}
    Wt::Matrix{Float64}              # [sdVec X_T]  (M × k+1)
    eta_model::JuMP.Model
    eta_x::Vector{VariableRef}
    eta_con::Vector{<:ConstraintRef}
    max_model::JuMP.Model
    max_lam::Vector{VariableRef}
end

function LPWorkspace(X_T::AbstractMatrix, sigma::AbstractMatrix)
    M = size(sigma, 1)
    k = size(X_T, 2)
    sdVec = sqrt.(diag(sigma))
    Wt = hcat(sdVec, Matrix{Float64}(X_T))
    C = -Wt

    em = _lp_model()
    ex = @variable(em, [1:(k+1)])
    econ = [@constraint(em, sum(C[i, j] * ex[j] for j in 1:(k+1)) <= 0.0) for i in 1:M]
    @objective(em, Min, 1.0 * ex[1])

    mm = _lp_model()
    ml = @variable(mm, [1:M], lower_bound = 0.0)
    beq = vcat(1.0, zeros(k))
    for i in 1:(k+1)
        @constraint(mm, sum(Wt[r, i] * ml[r] for r in 1:M) == beq[i])
    end

    return LPWorkspace(sdVec, Wt, em, collect(ex), econ, mm, collect(ml))
end

function _solve_eta!(ws::LPWorkspace, y_T)
    b = .-_asvec(y_T)
    for i in eachindex(ws.eta_con)
        set_normalized_rhs(ws.eta_con[i], b[i])
    end
    optimize!(ws.eta_model)
    ok = _is_optimal(ws.eta_model)
    eta = ok ? objective_value(ws.eta_model) : NaN
    xval = ok ? JuMP.value.(ws.eta_x) : fill(NaN, length(ws.eta_x))
    lambda = ok ? -JuMP.dual.(ws.eta_con) : fill(NaN, length(ws.eta_con))
    return (eta_star=eta, delta_star=xval[2:end], lambda=lambda, error_flag=(ok ? 0 : 1))
end

function _solve_max!(ws::LPWorkspace, f)
    @objective(ws.max_model, Max, sum(f[r] * ws.max_lam[r] for r in eachindex(f)))
    optimize!(ws.max_model)
    ok = _is_optimal(ws.max_model)
    return (objective=(ok ? objective_value(ws.max_model) : NaN),
            solution=(ok ? JuMP.value.(ws.max_lam) : fill(NaN, length(ws.max_lam))), ok=ok)
end

# .test_delta_lp_fn : min_{eta,delta} eta s.t. y_T - X_T delta <= eta*sdVec.
function test_delta_lp_fn(y_T, X_T, sigma; ws=nothing)
    ws !== nothing && return _solve_eta!(ws, y_T)
    yv = _asvec(y_T)
    Xt = X_T isa AbstractVector ? reshape(X_T, :, 1) : X_T
    dimDelta = size(Xt, 2)
    sdVec = sqrt.(diag(sigma))
    f = vcat(1.0, zeros(dimDelta))
    C = -hcat(sdVec, Xt)
    b = -yv

    model = _lp_model()
    @variable(model, x[1:(dimDelta+1)])
    @constraint(model, con, C * x .<= b)
    @objective(model, Min, f' * x)
    optimize!(model)

    ok = _is_optimal(model)
    eta = ok ? objective_value(model) : NaN
    xval = ok ? JuMP.value.(x) : fill(NaN, dimDelta + 1)
    delta = xval[2:end]
    lambda = ok ? -JuMP.dual.(con) : fill(NaN, length(b))
    return (eta_star=eta, delta_star=delta, lambda=lambda, error_flag=(ok ? 0 : 1))
end

_roundeps(x; eps_=eps(Float64)^(3 / 4)) = abs(x) < eps_ ? 0.0 : x

# max_{lambda>=0} f'lambda  s.t.  W_T' lambda = e_1, where
# f = s_T + (gamma' sigma gamma)^(-1) (sigma gamma) c. Returns the max and the lambda.
function max_program(s_T, gamma_tilde, sigma, W_T, c; ws=nothing)
    f = s_T .+ (1.0 / _scalar(gamma_tilde' * sigma * gamma_tilde)) .* (sigma * gamma_tilde) .* c
    # Degenerate dual (gamma_tilde ~ 0) yields non-finite objective coefficients.
    # R's GLPK tolerates this and the bisection falls through; we mirror that by
    # reporting "not a solution" so the caller early-returns gracefully.
    if !all(isfinite, f)
        return (objective=NaN, solution=fill(NaN, size(W_T, 1)), ok=false)
    end
    ws !== nothing && return _solve_max!(ws, f)
    nrowAeq = size(W_T, 2)
    beq = vcat(1.0, zeros(nrowAeq - 1))
    model = _lp_model()
    @variable(model, lam[1:size(W_T, 1)] >= 0)
    @constraint(model, W_T' * lam .== beq)
    @objective(model, Max, f' * lam)
    optimize!(model)
    ok = _is_optimal(model)
    return (objective=(ok ? objective_value(model) : NaN), solution=(ok ? JuMP.value.(lam) : fill(NaN, size(W_T, 1))), ok=ok)
end

function check_if_solution_helper(c, tol, s_T, gamma_tilde, sigma, W_T; ws=nothing)
    lp = max_program(s_T, gamma_tilde, sigma, W_T, c; ws=ws)
    honest = lp.ok && (abs(c - lp.objective) <= tol)
    return (objective=lp.objective, solution=lp.solution, honestsolution=honest)
end

# .vlo_vup_dual_fn : bisection (ARP 2021, Appendix D, Algorithm 1).
function vlo_vup_dual_fn(eta, s_T, gamma_tilde, sigma, W_T; ws=nothing)
    tol_c = 1e-6; tol_equality = 1e-6
    sigma_B = sqrt(_scalar(gamma_tilde' * sigma * gamma_tilde))
    low_initial = min(-100.0, eta - 20 * sigma_B)
    high_initial = max(100.0, eta + 20 * sigma_B)
    maxiters = 10000; switchiters = 10
    bvec = (1.0 / _scalar(gamma_tilde' * sigma * gamma_tilde)) .* (sigma * gamma_tilde)

    checksol = check_if_solution_helper(eta, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws).honestsolution
    if !checksol
        return (vlo=eta, vup=Inf)
    end

    # --- vup ---
    local vup
    if check_if_solution_helper(high_initial, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws).honestsolution
        vup = Inf
    else
        dif = 0.0; iters = 1
        lp = check_if_solution_helper(high_initial, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws)
        mid = _roundeps(_scalar(lp.solution' * s_T)) / (1 - _scalar(lp.solution' * bvec))
        while iters < maxiters
            lp = check_if_solution_helper(mid, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws)
            lp.honestsolution && break
            iters += 1
            if iters >= switchiters
                dif = tol_c + 1
                break
            end
            mid = _roundeps(_scalar(lp.solution' * s_T)) / (1 - _scalar(lp.solution' * bvec))
        end
        low = eta; high = mid
        while dif > tol_c && iters < maxiters
            iters += 1
            mid = (high + low) / 2
            if check_if_solution_helper(mid, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws).honestsolution
                low = mid
            else
                high = mid
            end
            dif = high - low
        end
        vup = mid
    end

    # --- vlo ---
    local vlo
    if check_if_solution_helper(low_initial, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws).honestsolution
        vlo = -Inf
    else
        dif = 0.0; iters = 1
        lp = check_if_solution_helper(low_initial, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws)
        mid = _roundeps(_scalar(lp.solution' * s_T)) / (1 - _scalar(lp.solution' * bvec))
        while iters < maxiters
            lp = check_if_solution_helper(mid, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws)
            lp.honestsolution && break
            iters += 1
            if iters >= switchiters
                dif = tol_c + 1
                break
            end
            mid = _roundeps(_scalar(lp.solution' * s_T)) / (1 - _scalar(lp.solution' * bvec))
        end
        low = mid; high = eta
        while dif > tol_c && iters < maxiters
            mid = (low + high) / 2
            iters += 1
            if check_if_solution_helper(mid, tol_equality, s_T, gamma_tilde, sigma, W_T; ws=ws).honestsolution
                high = mid
            else
                low = mid
            end
            dif = high - low
        end
        vlo = mid
    end

    return (vlo=vlo, vup=vup)
end

# .lp_dual_fn
function lp_dual_fn(y_T, X_T, eta, gamma_tilde, sigma; ws=nothing)
    yv = _asvec(y_T)
    Xt = X_T isa AbstractVector ? reshape(X_T, :, 1) : X_T
    sdVec = sqrt.(diag(sigma))
    W_T = hcat(sdVec, Xt)
    gt = _asvec(gamma_tilde)
    s_T = (Matrix{Float64}(I, length(yv), length(yv)) -
           (1.0 / _scalar(gt' * sigma * gt)) .* (sigma * (gt * gt'))) * yv
    v = vlo_vup_dual_fn(eta, s_T, gt, sigma, W_T; ws=ws)
    return (vlo=v.vlo, vup=v.vup, eta=eta, gamma_tilde=gt)
end

# .FLCI_computeVloVup
function FLCI_computeVloVup(vbar, dbar, S, c)
    vb = _asvec(vbar)
    VbarMat = vcat(reshape(vb, 1, :), reshape(-vb, 1, :))
    VbarS = VbarMat * _asvec(S)
    Vbarc = VbarMat * _asvec(c)
    max_or_min = (_asvec(dbar) .- VbarS) ./ Vbarc
    vlo = maximum(max_or_min[Vbarc .< 0])
    vup = minimum(max_or_min[Vbarc .> 0])
    return (vlo=vlo, vup=vup)
end

# .compute_least_favorable_cv
function compute_least_favorable_cv(X_T, sigma, hybrid_kappa; sims=1000, rowsForARP=nothing, seed=0)
    if rowsForARP !== nothing
        if X_T !== nothing
            X_T = X_T isa AbstractVector ? reshape(X_T[rowsForARP], :, 1) : X_T[rowsForARP, :]
        end
        sigma = sigma[rowsForARP, rowsForARP]
    end
    rng = MersenneTwister(seed)
    if X_T === nothing
        draws = _rmvnorm(rng, sigma, sims)                     # sims × M
        sdv = sqrt.(diag(sigma))
        xi = draws ./ reshape(sdv, 1, :)
        eta_vec = vec(maximum(xi; dims=2))
        return quantile(eta_vec, 1 - hybrid_kappa)
    else
        Xt = X_T isa AbstractVector ? reshape(X_T, :, 1) : X_T
        sdVec = sqrt.(diag(sigma))
        dimDelta = size(Xt, 2)
        f = vcat(1.0, zeros(dimDelta))
        C = -hcat(sdVec, Xt)
        draws = _rmvnorm(rng, sigma, sims)
        eta_vec = Float64[]
        for i in 1:sims
            b = -draws[i, :]
            push!(eta_vec, _compute_eta_lp(f, C, b))
        end
        eta_vec = filter(!isnan, eta_vec)
        return quantile(eta_vec, 1 - hybrid_kappa)
    end
end

function _compute_eta_lp(f, C, b)
    model = _lp_model()
    @variable(model, x[1:length(f)])
    @constraint(model, C * x .<= b)
    @objective(model, Min, f' * x)
    optimize!(model)
    return _is_optimal(model) ? objective_value(model) : NaN
end

# Robust multivariate-normal sampler (eigen method, handles PSD sigma).
function _rmvnorm(rng, sigma, n)
    S = Symmetric((Matrix(sigma) .+ Matrix(sigma)') ./ 2)
    E = eigen(S)
    L = E.vectors * Diagonal(sqrt.(max.(E.values, 0.0)))
    Z = randn(rng, n, size(sigma, 1))
    return Z * L'
end

# .lp_conditional_test_fn
function lp_conditional_test_fn(theta, y_T, X_T, sigma, alpha, hybrid_flag, hybrid_list, rowsForARP; ws=nothing)
    yv = _asvec(y_T)
    Xt = X_T isa AbstractVector ? reshape(X_T, :, 1) : X_T
    y_T_ARP = yv[rowsForARP]
    X_T_ARP = Xt[rowsForARP, :]
    sigma_ARP = sigma[rowsForARP, rowsForARP]

    M = size(sigma_ARP, 1)
    k = size(X_T_ARP, 2)

    linSoln = test_delta_lp_fn(y_T_ARP, X_T_ARP, sigma_ARP; ws=ws)
    if linSoln.error_flag > 0
        @warn "LP for eta did not converge properly. Not rejecting"
        return (reject=0, eta=linSoln.eta_star, delta=linSoln.delta_star)
    end

    if hybrid_flag == "LF"
        mod_size = (alpha - hybrid_list[:hybrid_kappa]) / (1 - hybrid_list[:hybrid_kappa])
        if linSoln.eta_star > hybrid_list[:lf_cv]
            return (reject=1, eta=linSoln.eta_star, delta=linSoln.delta_star)
        end
    elseif hybrid_flag == "FLCI"
        mod_size = (alpha - hybrid_list[:hybrid_kappa]) / (1 - hybrid_list[:hybrid_kappa])
        vbar = _asvec(hybrid_list[:vbar])
        VbarMat = vcat(reshape(vbar, 1, :), reshape(-vbar, 1, :))
        if maximum(VbarMat * yv - _asvec(hybrid_list[:dbar])) > 0
            return (reject=1, eta=linSoln.eta_star, delta=linSoln.delta_star)
        end
    elseif hybrid_flag == "ARP"
        mod_size = alpha
    else
        error("Hybrid flag must equal 'LF', 'FLCI' or 'ARP'")
    end

    tol_lambda = 1e-6
    lambda = linSoln.lambda
    degenerate_flag = (sum(lambda .> tol_lambda) != (k + 1))
    B_index = lambda .> tol_lambda
    Bc_index = .!B_index
    X_TB = reshape(X_T_ARP[B_index, :], :, size(X_T_ARP, 2))
    Xdim = minimum(size(X_TB))
    fullRank_flag = (Xdim == 0) ? false : (rank(X_TB) == Xdim)

    if !fullRank_flag || degenerate_flag
        # --- Dual approach ---
        lpDualSoln = lp_dual_fn(y_T_ARP, X_T_ARP, linSoln.eta_star, lambda, sigma_ARP; ws=ws)
        gamma_tilde = lpDualSoln.gamma_tilde
        sigma_B_dual2 = _scalar(gamma_tilde' * sigma_ARP * gamma_tilde)
        if abs(sigma_B_dual2) < eps(Float64)
            return (reject=(linSoln.eta_star > 0 ? 1 : 0), eta=linSoln.eta_star, delta=linSoln.delta_star, lambda=lambda)
        elseif sigma_B_dual2 < 0
            error(".vlo_vup_dual_fn returned a negative variance")
        end
        sigma_B_dual = sqrt(sigma_B_dual2)
        maxstat = lpDualSoln.eta / sigma_B_dual

        if hybrid_flag == "LF"
            zlo = lpDualSoln.vlo / sigma_B_dual
            zup = min(lpDualSoln.vup, hybrid_list[:lf_cv]) / sigma_B_dual
        elseif hybrid_flag == "FLCI"
            gamma_full = zeros(length(yv)); gamma_full[rowsForARP] = gamma_tilde
            sigma_gamma = (sigma * gamma_full) ./ _scalar(gamma_full' * sigma * gamma_full)
            S = yv - sigma_gamma .* _scalar(gamma_full' * yv)
            vFLCI = FLCI_computeVloVup(hybrid_list[:vbar], hybrid_list[:dbar], S, sigma_gamma)
            zlo = max(lpDualSoln.vlo, vFLCI.vlo) / sigma_B_dual
            zup = min(lpDualSoln.vup, vFLCI.vup) / sigma_B_dual
        else
            zlo = lpDualSoln.vlo / sigma_B_dual
            zup = lpDualSoln.vup / sigma_B_dual
        end

        if !(zlo <= maxstat <= zup)
            return (reject=0, eta=linSoln.eta_star, delta=linSoln.delta_star, lambda=lambda)
        else
            cval = max(0.0, norminvp_generalized(1 - mod_size, zlo, zup))
            return (reject=Int(maxstat > cval), eta=linSoln.eta_star, delta=linSoln.delta_star, lambda=lambda)
        end
    else
        # --- Primal approach ---
        size_B = sum(B_index)
        sdVec = sqrt.(diag(sigma_ARP))
        sdVec_B = sdVec[B_index]
        sdVec_Bc = sdVec[Bc_index]
        X_TBc = reshape(X_T_ARP[Bc_index, :], :, size(X_T_ARP, 2))
        Imat = Matrix{Float64}(I, M, M)
        S_B = Imat[B_index, :]
        S_Bc = Imat[Bc_index, :]

        WB = hcat(sdVec_B, X_TB)
        Gamma_B = hcat(sdVec_Bc, X_TBc) * inv(WB) * S_B - S_Bc
        e1 = vcat(1.0, zeros(size_B - 1))
        v_B = vec((e1' * inv(WB) * S_B))
        sigma2_B = _scalar(v_B' * sigma_ARP * v_B)
        sigma_B = sqrt(sigma2_B)
        rho = (Gamma_B * sigma_ARP * v_B) ./ sigma2_B
        maximand_or_minimand = (-(Gamma_B * y_T_ARP)) ./ rho .+ _scalar(v_B' * y_T_ARP)

        vlo = sum(rho .> 0) > 0 ? maximum(maximand_or_minimand[rho .> 0]) : -Inf
        vup = sum(rho .< 0) > 0 ? minimum(maximand_or_minimand[rho .< 0]) : Inf

        if hybrid_flag == "LF"
            zlo = vlo / sigma_B
            zup = min(vup, hybrid_list[:lf_cv]) / sigma_B
        elseif hybrid_flag == "FLCI"
            gamma_full = zeros(length(yv)); gamma_full[rowsForARP] = v_B
            sigma_gamma = (sigma * gamma_full) ./ _scalar(gamma_full' * sigma * gamma_full)
            S = yv - sigma_gamma .* _scalar(gamma_full' * yv)
            vFLCI = FLCI_computeVloVup(hybrid_list[:vbar], hybrid_list[:dbar], S, sigma_gamma)
            zlo = max(vlo, vFLCI.vlo) / sigma_B
            zup = min(vup, vFLCI.vup) / sigma_B
        else
            zlo = vlo / sigma_B
            zup = vup / sigma_B
        end

        maxstat = linSoln.eta_star / sigma_B
        if !(zlo <= maxstat <= zup)
            return (reject=0, eta=linSoln.eta_star, delta=linSoln.delta_star, lambda=lambda)
        else
            cval = max(0.0, norminvp_generalized(1 - mod_size, zlo, zup))
            return (reject=Int(maxstat > cval), eta=linSoln.eta_star, delta=linSoln.delta_star, lambda=lambda)
        end
    end
end

test_delta_lp_fn_wrapper(theta, y_T, X_T, sigma, alpha, hybrid_flag, hybrid_list, rowsForARP; ws=nothing) =
    lp_conditional_test_fn(theta, y_T, X_T, sigma, alpha, hybrid_flag, hybrid_list, rowsForARP; ws=ws).reject

# .ARP_computeCI
function ARP_computeCI(betahat, sigma, numPrePeriods, numPostPeriods, A, d, l_vec, alpha,
                       hybrid_flag, hybrid_list, returnLength, grid_lb, grid_ub, gridPoints, rowsForARP)
    bh = _asvec(betahat); dv = _asvec(d); lv = _asvec(l_vec)
    rows = rowsForARP === nothing ? collect(1:size(A, 1)) : collect(rowsForARP)
    thetaGrid = collect(range(grid_lb, grid_ub; length=gridPoints))

    Gamma = construct_Gamma(lv)
    Apost = A[:, (numPrePeriods+1):(numPrePeriods+numPostPeriods)]
    AGammaInv = Apost * inv(Gamma)
    AGammaInv_one = AGammaInv[:, 1]
    AGammaInv_minusOne = AGammaInv[:, 2:end]

    Y = A * bh - dv
    sigmaY = A * sigma * A'

    if hybrid_flag == "LF"
        hybrid_list[:lf_cv] = compute_least_favorable_cv(AGammaInv_minusOne, sigmaY,
                                                         hybrid_list[:hybrid_kappa]; rowsForARP=rows)
    end

    # Build the reusable LP models once (constraint structure is fixed across theta).
    ws = LPWorkspace(AGammaInv_minusOne[rows, :], sigmaY[rows, rows])

    accept = Float64[]
    for theta in thetaGrid
        if hybrid_flag == "FLCI"
            vbar = _asvec(hybrid_list[:vbar])
            hl = hybrid_list[:flci_halflength]
            hybrid_list[:dbar] = [hl - _scalar(vbar' * dv) + (1 - _scalar(vbar' * AGammaInv_one)) * theta,
                                  hl + _scalar(vbar' * dv) - (1 - _scalar(vbar' * AGammaInv_one)) * theta]
        end
        reject = test_delta_lp_fn_wrapper(theta, Y - AGammaInv_one .* theta, AGammaInv_minusOne,
                                          sigmaY, alpha, hybrid_flag, hybrid_list, rows; ws=ws)
        push!(accept, 1.0 - reject)
    end

    if (accept[1] == 1 || accept[end] == 1) && hybrid_flag != "FLCI"
        @warn "CI is open at one of the endpoints; CI length may not be accurate"
    end

    if returnLength
        gridLength = 0.5 .* (vcat(0.0, diff(thetaGrid)) .+ vcat(diff(thetaGrid), 0.0))
        return sum(accept .* gridLength)
    else
        return ConditionalCS(thetaGrid, accept)
    end
end
