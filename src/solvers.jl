# Thin JuMP wrappers. LPs use HiGHS (robust + exact dual/sensitivity info,
# matching R's lpSolveAPI / Rglpk); conic (SOCP/QP) problems use Clarabel
# (pure-Julia, replaces R's CVXR + ECOS).

function _lp_model()
    # direct_model avoids JuMP's MOI caching layer; presolve off is fastest for
    # the tiny LPs solved here (and they are rebuilt many times in the ARP loop).
    m = direct_model(HiGHS.Optimizer())
    set_silent(m)
    set_attribute(m, "presolve", "off")
    return m
end

# Conic (SOCP/QP) problems use ECOS - the same solver the R package uses via
# CVXR, which matches its behaviour at degenerate boundaries (e.g. h = hMin).
function _conic_model()
    m = Model(ECOS.Optimizer)
    set_silent(m)
    return m
end

# A usable primal solution exists (accept near-boundary points that stall at
# SLOW_PROGRESS, mirroring how R reads whatever ECOS returns).
_has_solution(m) = JuMP.primal_status(m) == MOI.FEASIBLE_POINT

# Status helpers mirroring R's interpretation of solver success.
_is_optimal(m) = termination_status(m) == MOI.OPTIMAL
function _status_string(m)
    st = termination_status(m)
    if st == MOI.OPTIMAL
        return "optimal"
    elseif st == MOI.ALMOST_OPTIMAL || st == MOI.LOCALLY_SOLVED
        return "optimal_inaccurate"
    else
        return lowercase(string(st))
    end
end
