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

"""
    optimize(params::AdmmParams[, settings::Settings[, vars::Vars]])

Carries out the ADMM algorithm by calling `admm_step` until termination conditions are met.
"""
function optimize(params::AdmmParams, settings::Settings, vars::Vars)
    assert_valid(params, settings, vars)
    kkt = factorize_kkt(params, settings)
    pc = ProxCache(params, settings)
    stats = Stats(settings)
    term_cache = get_term_cond_cache(params, settings)
    iter = 0
    while iter < settings.max_iters
        iter += 1
        admm_step!(vars, params, kkt, pc, settings, iter)
        if settings.compute_stats
            compute_stats(vars, params, stats, iter)
        end
        if check_term_conds!(term_cache, vars, params, settings, iter)
            break
        end
    end
    return vars, stats
end
optimize(params::AdmmParams, settings::Settings) = optimize(params, settings, Vars(params))
optimize(params::AdmmParams) = optimize(params, Settings(params))

@doc raw"""
    factorize_kkt(params::AdmmParams, settings::Settings)

Factorize the KKT matrix representing the optimzality equations of the provided problem.

# Extended help
Factorize the coefficient matrix of that KKT system representing the optimality equations
that will be used (and reused) for the constrained least squares problem:

```math
\begin{array}{ll}
\text{minimize} & \displaystyle
\frac12 x^T P x  + q^T x + \sum_{i=1}^n g_i(x_i) \\
\text{subject to} & Ax = b.
\end{array}
```

We can formulate this equivalently as

```math
\begin{array}{ll}
\text{minimize} & \displaystyle
\frac12 \tilde x^T P \tilde x  + q^T \tilde x + \mathcal{I}_\mathcal{A}(\tilde x, \tilde z)
        + g(x) + \mathcal{I}_{\mathcal B}(z) \\
\text{subject to} & \tilde x = x \\
                  & \tilde z = z.
\end{array}
```

where ``\mathcal I_{\mathcal A}`` is the indicator function over the affine
constraints of the original form, *i.e.*,
```math
\mathcal{I_{\mathcal A}}(x, z) = \begin{cases}
0 & Ax = z \\
\infty & {\rm otherwise}.
\end{cases}
```

and ``\mathcal I_{\mathcal B}`` is the indicator function over the singleton
set of vectors equal to ``b``, *i.e.*,

```math
\mathcal{I_{\mathcal B}}(x) = \begin{cases}
0 & x = b \\
\infty & {\rm otherwise}.
\end{cases}
```

We will denote the objective in this forumulation as ``f(x, \tilde x, z, \tilde z)``.

The augmented Lagrangian for this problem is

```math
L_{S,R} (x, \tilde x, z, \tilde z, w, y) =
f(x, \tilde x, z, \tilde z)
+ \frac12 \| \tilde x - x + S^{-1} w \|_S^2
+ \frac12 \| \tilde z - z + R^{-1} y \|_R^2,
```

where ``f`` is the previously mentioned objective function.
The diagonal, positive definite matrices ``S\in \bf S{++}^l`` and ``R\in \bf S_{++}^m``
are algorithm parameters.
We define ``\|\cdot \|_S`` is the ``S``-norm, *i.e.*,
``\|x\|_S = \sqrt{x^T S x}``.

Solving for the optimality conditions, we can find optimal ``\tilde x`` and ``\tilde z`` by solving the linear system

```math
\begin{bmatrix} P + S & A^T \\ A & -R^{-1} \end{bmatrix}
\begin{bmatrix} \tilde x \\ \nu \end{bmatrix}
=
\begin{bmatrix} Sx - q - w \\ z - R^{-1} y \end{bmatrix}.
```

where ``\tilde z`` can be recovered as ``\tilde z = z + R^{-1}(\nu - y)``

In this function, we factor the coefficient matrix on the left hand side so that we repeatedly
solve the linear system efficiently.
"""
function factorize_kkt(params::AdmmParams, settings::Settings)
    AT = sparse(params.A')
    Dσ = sparse(Diagonal(settings.σ))
    Dρ = sparse(Diagonal(1 ./ settings.ρ))

    # construct the block matrix K
    K = [params.P+Dσ AT;
         params.A -Dρ]
    return ldl(K)
end

"""
    solve_kkt!(vars::Vars, kkt::KKT, params::AdmmParams, settings::Settings)

Solve the KKT system set up in the description of `factorize_kkt` for a particular right hand side.
"""
function solve_kkt!(vars::Vars, kkt::KKT, params::AdmmParams, settings::Settings)
    m, n = size(params.A)
    rhs = Vector{Float64}(undef, n + m)
    for i in 1:n
        rhs[i] = settings.σ[i] * vars.x[i] - vars.w[i] - params.q[i]
    end
    for i in 1:m
        rhs[i + n] = vars.z[i] - (vars.y[i] / settings.ρ[i])
    end
    lhs = kkt \ rhs
    xt = lhs[1:n]
    ν = lhs[(n + 1):end]
    zt = vars.z + (ν - vars.y) ./ settings.ρ
    return xt, zt
end

"""
    admm_step!(vars::Vars, params::AdmmParams, kkt::KKT, pc::ProxCache, settings::Settings, iter::Int64; polish::Bool=false)

Carry out a single ADMM iteration.

# Steps:
1. Solve the KKT system to carry out the `xt` and `zt` updates.
2. Use a mixture of the result of (1) and the previous `x` iterate to solve `n` univariate
   subproblems by minimizing the Lagrangian in the description of `factorize_kkt` w.r.t.
   `x`. (Requires evaluation of a proximal operator.)
3. Update `w` and `y` using the residuals associated with the constraints `x = xt` and `z = zt`.

# Note:
* Note that for the entire run of the algorithm, the `z` update is just `z = b`, so we don't need
to do anything here.
"""
function admm_step!(vars::Vars, params::AdmmParams, kkt::KKT, pc::ProxCache,
                    settings::Settings, iter::Int64; polish::Bool=false)
    # update xt and zt by solving the KKT system
    xt, zt = solve_kkt!(vars, kkt, params, settings)

    # get the next x value by using the proximal operators of the separable function pieces
    prox_target = settings.α * xt + (1 - settings.α) * vars.x + vars.w ./ settings.σ
    x = prox_step(pc, params, settings, prox_target)
    @assert all(isfinite.(x))

    # update the dual variables
    z = params.b
    w = vars.w + settings.σ .* (settings.α * xt + (1 - settings.α) * vars.x - x)
    y = vars.y + settings.ρ .* (settings.α * zt + (1 - settings.α) * vars.z - z)

    # save in vars struct
    vars.xt = xt
    vars.zt = zt
    vars.x = x
    vars.w = w
    vars.z = z
    return vars.y = y
end

"""
    compute_stats(vars::Vars, params::AdmmParams, stats::Stats, iter::Int64)

Update the objective, residual, and iteration fields of the Stats object.
"""
function compute_stats(vars::Vars, params::AdmmParams, stats::Stats, iter::Int64)
    x = vars.x
    P = params.P
    q = params.q
    A = params.A
    b = params.b
    g = params.g

    fx = 0.5 * x' * P * x + q' * x
    gx = 0.0
    for i in 1:get_num_vars(params)
        if !isfinite(g[i, 1](x[i]))
        else
            gx += g[i, 1](x[i])
        end
    end
    stats.obj[iter] = fx + gx
    stats.iters = iter

    res = norm(A * x - b)
    return stats.res[iter] = res
end
