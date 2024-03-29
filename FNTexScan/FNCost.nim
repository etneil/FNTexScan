import FNTheory
import math
import arraymancer
import gsl/gsl_min
import gsl/gsl_errno
import gsl/gsl_multimin

## Quark masses from 1112.3112
const Md_SM_ref = [2.4e-3, 0.049, 2.41]
const err_Md_SM_ref = [0.42e-3, 0.015, 0.14]
const Mu_SM_ref = [1.17e-3, 0.543, 148.1]
const err_Mu_SM_ref = [0.35e-3, 0.072, 1.3]

## CKM entries from PDG review
#     [V12, V23, V13]:
const CKM_SM_ref = [0.2245, 0.0410, 0.00382]
const err_CKM_SM_ref = [0.00067, 0.00085, 0.00011]
const J_SM_ref = 3.08e-5
const err_J_SM_ref = 0.14e-5


proc naturalness_cost_eps(eps: float, eps_max = 0.3): float =
    result = exp(10 * (eps - eps_max))

proc naturalness_cost*(thy: FNTheory, 
    flat_prior = false, sigma_weight = 0.691): float =

    var Cu_mag_sq, Cd_mag_sq: float
    var Cu_cost, Cd_cost = 0.0

    let sigma_fac = 1/(4*sigma_weight^2)

    for i in 0 ..< 9:
        # Log(|C|) should be distributed as chi^2 = Log(|C|)^2 / sigma^2
        # = log(|C|^2)^2 / (4*sigma^2)
        Cu_mag_sq = thy.C[i]^2 + thy.C[i+9]^2
        Cd_mag_sq = thy.C[i+18]^2 + thy.C[i+27]^2

        Cu_cost += ln(Cu_mag_sq)^2 * sigma_fac
        Cd_cost += ln(Cd_mag_sq)^2 * sigma_fac

    result = Cu_cost + Cd_cost

proc check_max_exp(current, x, x_ref: float): float =
    if x <= 0.0:
        return current + 1e6

    var trial_res = exp(abs(ln(x / x_ref)))

    if trial_res > current:
        result = trial_res
    else:
        result = current

proc score_SM_deviations*(thy: FNTheory, max_exp = false): float =
    # Score on: Mu, Md, |V12|, |V23|, |V13|, and J
    let SM_params = thy.getSMParams()

    if max_exp:
        result = -1.0

    var Mu, Md: array[3, float64]

    case thy.hierarchy
    of normal:
        Mu[0] = SM_params[0]
        Md[0] = SM_params[3]
    of inverted:
        Mu[0] = SM_params[3]
        Md[0] = SM_params[0]
    of degen:
        Mu[0] = (SM_params[0] + SM_params[3]) / 2
        Md[0] = (SM_params[0] + SM_params[3]) / 2

    Mu[1] = SM_params[1]
    Mu[2] = SM_params[2]
    Md[1] = SM_params[4]
    Md[2] = SM_params[5]

    for i in 0..<3:
        if max_exp:
            result = check_max_exp(result, Mu[i], Mu_SM_ref[i])
            result = check_max_exp(result, Md[i], Md_SM_ref[i])
        else:
            result += (Mu[i] - Mu_SM_ref[i])^2 / err_Mu_SM_ref[i]^2
            result += (Md[i] - Md_SM_ref[i])^2 / err_Md_SM_ref[i]^2

    let J = SM_params[24]
    if max_exp:
#        discard
        result = check_max_exp(result, J, J_SM_ref)
    else:
        result += (J - J_SM_ref)^2 / err_J_SM_ref^2

    let V12 = sqrt(SM_params[7]^2 + SM_params[16]^2)
    let V23 = sqrt(SM_params[11]^2 + SM_params[20]^2)
    let V13 = sqrt(SM_params[8]^2 + SM_params[17]^2)

    for i, V in [V12, V23, V13]:
        if max_exp:
            result = check_max_exp(result, V, CKM_SM_ref[i])
        else:
            result += (V - CKM_SM_ref[i])^2 / err_CKM_SM_ref[i]^2



proc optimize_eps*(thy: FNTheory, brute_force = true, 
                  max_exp = true, verbose = false): float =
    var optThy = FNTheory(tex:thy.tex, C:thy.C)

    if brute_force:
        let leps_range = linspace(-3.0, -0.6, 300)
        var min_cost, min_leps, cost: float

        min_leps = leps_range[0]
        min_cost = naturalness_cost_eps(exp(min_leps))
        min_cost += score_SM_deviations(optThy, max_exp=max_exp)

        for leps in leps_range:
            optThy.C[^1] = leps
            cost = naturalness_cost_eps(exp(leps)) 
            cost += score_SM_deviations(optThy, max_exp=max_exp)

            if cost < min_cost:
                min_cost = cost
                min_leps = leps

        result = exp(min_leps)

    else:
        # GSL minimization
        proc eps_opt_function(leps: cdouble, params: pointer): cdouble {.cdecl.} =
            var thyPtr = cast[ptr FNTheory](params)
            var thisThy : FNTheory = thyPtr[]
            thisThy.C[^1] = leps
            result = naturalness_cost_eps(exp(leps))
            result += score_SM_deviations(thisThy, max_exp=true)
        
        let f = gsl_min.gsl_function(
            function: eps_opt_function,
            params: addr optThy
        )

        discard gsl_errno.gsl_set_error_handler_off()

        var
            status = -2
            iter = 0
            max_iter = 500
            tol_abs = 1e-2
            lower = -8.0
            upper = 1.0
            leps_min = ln(0.02)

#        let T = gsl_min_fminimizer_brent
        let T = gsl_min.gsl_min_fminimizer_quad_golden
        var s = gsl_min_fminimizer_alloc(T)

        discard gsl_min_fminimizer_set(s, unsafeAddr f, leps_min, lower, upper)

        while status == -2 and iter <= max_iter:
            iter += 1
            discard gsl_min_fminimizer_iterate(s)
            leps_min = gsl_min_fminimizer_x_minimum(s)
            lower = gsl_min_fminimizer_x_lower(s)
            upper = gsl_min_fminimizer_x_upper(s)

            status = gsl_min.gsl_min_test_interval(lower, upper, tol_abs, 0.0)

        gsl_min_fminimizer_free(s)

        # Make sure minimization was successful
        if iter >= max_iter and verbose:
            echo "Warning: reached max iteration count [", max_iter, "] in optimize_eps."
        else:
            #doAssert status == 0
            discard

        result = exp(leps_min)

proc optimize_SM*(thy: FNTheory): FNTheory =
    result = FNTheory(tex:thy.tex, C:thy.C)

    proc SM_opt_function(v: ptr gsl_vector, params: pointer): cdouble {.cdecl.} =
        var thyPtr = cast[ptr FNTheory](params)
        var thisThy: FNTheory = thyPtr[]

        # Replace C array with entries in v
        for i in 0 .. 36:
            thisThy.C[i] = gsl_vector_get(v, csize_t(i))

        result = naturalness_cost_eps(exp(thisThy.C[^1]))
        result += naturalness_cost(thisThy)
        result += score_SM_deviations(thisThy)

    const dim = 37

    let f = gsl_multimin.gsl_multimin_function(
        n: dim,
        f: SM_opt_function,
        params: addr result
    )

    discard gsl_errno.gsl_set_error_handler_off()


    var ss, x: ptr gsl_vector

    ss = gsl_vector_alloc(dim)
    gsl_vector_set_all(ss, 1.0)

    # Starting guess
    x = gsl_vector_alloc(dim)
    for i in 0 .. 36:
        gsl_vector_set(x, csize_t(i), result.C[i])

    let T = gsl_multimin_fminimizer_nmsimplex2rand
    var s = gsl_multimin_fminimizer_alloc(T, dim)

    var status = -2
    var iter = 0
    var size: cdouble
    var max_iter = 10000

    discard gsl_multimin_fminimizer_set(s, unsafeAddr f, x, ss)

    while status == -2 and iter <= max_iter:
        iter += 1
        discard gsl_multimin_fminimizer_iterate(s)

        size = gsl_multimin_fminimizer_size(s)
        status = gsl_multimin_test_size(size, 1e-1)

    for i in 0 .. 36:
        result.C[i] = gsl_vector_get(s.x, csize_t(i))

    gsl_vector_free(x)
    gsl_vector_free(ss)
    gsl_multimin_fminimizer_free(s)

    # Make sure minimization was successful
    if iter >= max_iter:
        echo "Warning: reached max iteration count [", max_iter, "] in optimize_SM."
    else:
        doAssert status == 0    

