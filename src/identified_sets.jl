# Identified-set bounds via LP, ported from the `.compute_IDset_Delta*`
# functions across the R delta files. Each solves
#   max / min  l'delta_post   s.t.  A delta <= d,  delta_pre = trueBeta_pre
# and returns (id.lb, id.ub) for theta = l'(beta_post - delta_post).

function _solve_idset_dir(fDelta, A, d, trueBeta, numPre, numPost; maximize::Bool)
    n = numPre + numPost
    model = _lp_model()
    @variable(model, x[1:n])
    if size(A, 1) > 0
        @constraint(model, A * x .<= d)
    end
    if numPre > 0
        @constraint(model, x[1:numPre] .== trueBeta[1:numPre])
    end
    if maximize
        @objective(model, Max, fDelta' * x)
    else
        @objective(model, Min, fDelta' * x)
    end
    optimize!(model)
    ok = _is_optimal(model)
    return (ok ? objective_value(model) : NaN, ok)
end

function _idset_bounds(A, d, M_unused, trueBeta, l_vec, numPre, numPost)
    lv = _asvec(l_vec)
    tb = _asvec(trueBeta)
    fDelta = vcat(zeros(numPre), lv)
    vmax, okmax = _solve_idset_dir(fDelta, A, d, tb, numPre, numPost; maximize=true)
    vmin, okmin = _solve_idset_dir(fDelta, A, d, tb, numPre, numPost; maximize=false)
    post = dot(lv, tb[(numPre+1):(numPre+numPost)])
    # R: warn + return point estimate only if BOTH directions fail.
    if !okmax && !okmin
        return (post, post)
    else
        return (post - vmax, post - vmin)  # (id.lb, id.ub)
    end
end

compute_IDset_DeltaSD(M, trueBeta, l_vec, numPre, numPost) =
    _idset_bounds(create_A_SD(numPre, numPost), create_d_SD(numPre, numPost, M),
                  M, trueBeta, l_vec, numPre, numPost)

compute_IDset_DeltaSDB(M, trueBeta, l_vec, numPre, numPost, biasDirection) =
    _idset_bounds(create_A_SDB(numPre, numPost; biasDirection=biasDirection),
                  create_d_SDB(numPre, numPost, M), M, trueBeta, l_vec, numPre, numPost)

compute_IDset_DeltaSDM(M, trueBeta, l_vec, numPre, numPost, monotonicityDirection) =
    _idset_bounds(create_A_SDM(numPre, numPost; monotonicityDirection=monotonicityDirection),
                  create_d_SDM(numPre, numPost, M), M, trueBeta, l_vec, numPre, numPost)

# RM / SDRM identified sets: union over the maximal period s and ± sign.
function compute_IDset_DeltaRM_fixedS(s, Mbar, max_positive, trueBeta, l_vec, numPre, numPost)
    A = create_A_RM(numPre, numPost; Mbar=Mbar, s=s, max_positive=max_positive)
    d = create_d_RM(numPre, numPost)
    return _idset_bounds(A, d, Mbar, trueBeta, l_vec, numPre, numPost)
end

function compute_IDset_DeltaRM(Mbar, trueBeta, l_vec, numPre, numPost)
    min_s = -(numPre - 1)
    lbs = Float64[]; ubs = Float64[]
    for s in min_s:0, mp in (true, false)
        lb, ub = compute_IDset_DeltaRM_fixedS(s, Mbar, mp, trueBeta, l_vec, numPre, numPost)
        push!(lbs, lb); push!(ubs, ub)
    end
    return (minimum(lbs), maximum(ubs))
end

function compute_IDset_DeltaSDRM_fixedS(s, Mbar, max_positive, trueBeta, l_vec, numPre, numPost)
    A = create_A_SDRM(numPre, numPost; Mbar=Mbar, s=s, max_positive=max_positive)
    d = create_d_SDRM(numPre, numPost)
    return _idset_bounds(A, d, Mbar, trueBeta, l_vec, numPre, numPost)
end

function compute_IDset_DeltaSDRM(Mbar, trueBeta, l_vec, numPre, numPost)
    min_s = -(numPre - 2)
    lbs = Float64[]; ubs = Float64[]
    for s in min_s:0, mp in (true, false)
        lb, ub = compute_IDset_DeltaSDRM_fixedS(s, Mbar, mp, trueBeta, l_vec, numPre, numPost)
        push!(lbs, lb); push!(ubs, ub)
    end
    return (minimum(lbs), maximum(ubs))
end
