#=
Copyright 2021 BlackRock, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

abstract type TermCache end

TERM_COND_FLAG = Int64
CONVERGE_TERM_COND_FLAG = 1
FIRST_VARS_TERM_COND_FLAG = 2

"""
    get_term_cond_cache(params::AdmmParams, settings::Settings)

Return the TermCache corresponding to the `term_cond_type` specification in the ADMM Settings `settings`
"""
function get_term_cond_cache(params::AdmmParams, settings::Settings)
    if settings.term_cond_type == CONVERGE_TERM_COND_FLAG
        return ConvergeTermCache(params, settings)
    elseif settings.term_cond_type == FIRST_VARS_TERM_COND_FLAG
        return FirstVarsTermCache(params, settings)
    end
end

"""
    ConvergeTermCache

Cache state to helps determine whether or not to terminate

The fields represent:
- `x_last`: Previous x
- `z_last`: Previous z
- `w_last`: Previous w
- `y_last`: Previous y
"""
mutable struct ConvergeTermCache
    x_last::Vector{Float64}
    z_last::Vector{Float64}
    w_last::Vector{Float64}
    y_last::Vector{Float64}

    function ConvergeTermCache(num_vars::Int64, num_constrs::Int64)
        return new(zeros(num_vars), zeros(num_constrs), zeros(num_vars), zeros(num_constrs))
    end
end

"""
    ConvergeTermCache(params::AdmmParams, settings::Settings)

Construct a ConvergeTermCache given AdmmParams `params` and ADMM Settings `settings`
"""
function ConvergeTermCache(params::AdmmParams, settings::Settings)
    return ConvergeTermCache(get_num_vars(params), get_num_constrs(params))
end

"""
    check_term_conds!(tc::ConvergeTermCache, vars::Vars, params::AdmmParams, settings::Settings, iter::Int)

Check termination conditions
"""
function check_term_conds!(tc::ConvergeTermCache, vars::Vars, params::AdmmParams,
                           settings::Settings, iter::Int)

    # only check for convergence every so often
    iter_condition = iter % settings.term_cond_freq == 0
    if iter_condition
        converged = ((maximum(abs.(vars.x - vars.xt)) ≤ settings.ϵ) &&
                     (maximum(abs.(tc.x_last - vars.x)) ≤ settings.ϵ) &&
                     (maximum(abs.(tc.w_last - vars.w)) ≤ settings.ϵ))
        if get_num_constrs(vars) > 0
            converged &= ((maximum(abs.(tc.z_last - vars.z)) ≤ settings.ϵ) &&
                          (maximum(abs.(vars.z - vars.zt)) ≤ settings.ϵ) &&
                          (maximum(abs.(tc.y_last - vars.y)) ≤ settings.ϵ))
        end

        if converged
            return true
        end
    end

    # Store variables for next step.
    tc.x_last = vars.x
    tc.z_last = vars.z
    tc.w_last = vars.w
    tc.y_last = vars.y

    return false
end

"""
    FirstVarsTermCache

Cache state to helps determine whether or not to terminate

The fields represent:
- `obj_best`: The best encountered objective value
- `res_best`: The best encountered residual value
- `x_best`: The best encountered `x` iterate
- `z_best`: The best encountered `z` iterate
- `w_best`: The best encountered `w` (dual) iterate
- `y_best`: The best encountered `y` (dual) iterate
- `n1`: `n` - `m`
- `lb`: Lower bound
- `ub`: Upper bound
- `A1`: first `n1` columns of `A`
- `A2`: remaining columns of `A`
- `not_improved_count`: The number of iterations without improvement seeing improvement
- `not_improved_count_req`: The max allowed iterations without improvement

"""
mutable struct FirstVarsTermCache
    obj_best::Float64
    res_best::Float64
    x_best::Vector{Float64}
    z_best::Vector{Float64}
    w_best::Vector{Float64}
    y_best::Vector{Float64}
    n1::Int64
    lb::Vector{Float64}
    ub::Vector{Float64}
    A1::SparseMatrixCSC
    A2::SparseMatrixCSC
    not_improved_count::Int64
    not_improved_count_req::Int64

    function FirstVarsTermCache(params::AdmmParams, settings::Settings)
        n = get_num_vars(params)
        m = get_num_constrs(params)

        obj_best = Inf
        res_best = Inf
        x_best = zeros(n)
        z_best = zeros(m)
        w_best = zeros(n)
        y_best = zeros(m)
        n1 = n - m
        lb = [params.g[i, 1][1].lb for i in (n1 + 1):n]
        ub = [params.g[i, 1][end].ub for i in (n1 + 1):n]
        A1 = params.A[:, 1:n1]
        A2 = params.A[:, (n1 + 1):end]
        not_improved_count = 0
        not_improved_count_req = ceil(Int64,
                                      settings.non_improvement_iters /
                                      settings.term_cond_freq)
        return new(obj_best, res_best, x_best, z_best, w_best, y_best, n1, lb, ub, A1, A2,
                   not_improved_count, not_improved_count_req)
    end
end

"""
    check_term_conds!(tc::FirstVarsTermCache, vars::Vars, params::AdmmParams, settings::Settings, iter::Int)

Check termination conditions
"""
function check_term_conds!(tc::FirstVarsTermCache, vars::Vars, params::AdmmParams,
                           settings::Settings, iter::Int)

    # only check for convergence every so often
    iter_condition = iter % settings.term_cond_freq == 0
    if iter_condition
        obj, res = get_obj_and_res_from_first_vars(tc, vars, params)

        # in order to have "improved", the objective value must decrease, the
        # residuals must be sufficiently small
        improved = ((obj < tc.obj_best - settings.obj_tol) && (res < settings.res_tol))
        if improved
            tc.not_improved_count = 0
        else
            tc.not_improved_count += 1
        end

        # if found a new best objective/residual pair, save the iterate
        if (obj < tc.obj_best) && (res < settings.res_tol)
            tc.x_best = copy(vars.x)
            tc.z_best = copy(vars.z)
            tc.w_best = copy(vars.w)
            tc.y_best = copy(vars.y)
            tc.obj_best = obj
            tc.res_best = res
        end

        # if there hasn't been sufficient improvement over the best for a
        # certain number of iterations, set the vars to contain the best iterate
        converged = (tc.not_improved_count > tc.not_improved_count_req) &&
                    (tc.obj_best < Inf)
        if converged
            vars.x = copy(tc.x_best)
            vars.z = copy(tc.z_best)
            vars.w = copy(tc.w_best)
            vars.y = copy(tc.y_best)
        end
    else
        converged = false
    end

    return converged
end

"""
    get_obj_and_res_from_first_vars(tc::FirstVarsTermCache, vars::Vars, params::AdmmParams)

Get objective and residual from the FirstVarsTermCache `tc`
"""
function get_obj_and_res_from_first_vars(tc::FirstVarsTermCache, vars::Vars,
                                         params::AdmmParams)
    x = guess_point_from_first_vars(vars, params)
    n = get_num_vars(params)

    fx = 0.5 * x' * params.P * x + params.q' * x  # quad plus lin
    gx = 0.0  # sep
    residual = 0.0
    for i in 1:get_num_vars(params)
        if isfinite(params.g[i, 1](x[i]))
            gx += params.g[i, 1](x[i])
        end
        if i > tc.n1
            if x[i] < tc.lb[i - tc.n1]
                residual += tc.lb[i - tc.n1] - x[i]
            elseif x[i] > tc.ub[i - tc.n1]
                residual += x[i] - tc.ub[i - tc.n1]
            end
        end
    end

    objective = fx + gx
    return objective, residual
end

function guess_point_from_first_vars(vars::Vars, params::AdmmParams)
    n = get_num_vars(vars)
    m = get_num_constrs(vars)
    n1 = n - m
    x1 = vars.x[1:n1]
    A1 = params.A[:, 1:n1]
    A2 = params.A[:, (n1 + 1):end]
    x2 = A2 \ (params.b - A1 * x1)
    x = [x1; x2]
    return x
end
