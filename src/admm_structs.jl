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

const KKT = LDLFactorizations.LDLFactorization{Float64,Int64}

"""
    AdmmParams

Parameters to the ADMM algorithm.

The fields represent:
- `P`: A sparse `n`x`n` positive semidefinite matrix
- `q`: an `n`-vector
- `A`: an `m`x`n` matrix
- `b`: an `m`-vector
- `g`: a vector of PiecewiseQuadratic functions
"""
mutable struct AdmmParams
    P::SparseMatrixCSC{Float64,Int64}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    g::Vector{PiecewiseQuadratic}

    function AdmmParams(P, q, A, b, g)
        m, n = size(A)
        @assert size(P) == (n, n)
        @assert size(q) == (n,)
        @assert size(A) == (m, n)
        @assert size(b) == (m,)
        @assert size(g) == (n,)
        return new(P, q, A, b, g)
    end
end

function Base.show(io::IO, params::AdmmParams)
    print(io, "AdmmParams")
    print(io, "  P: ")
    print(io, params.P)
    print(io, "  q: ")
    print(io, params.q)
    print(io, "  A: ")
    print(io, params.A)
    print(io, "  b: ")
    print(io, params.b)
    print(io, "  g: ")
    print(io, params.g)
    return
end

"""
    get_num_vars(params::AdmmParams)

Return the number of variables in the problem that `params` was initialized with.
"""
get_num_vars(params::AdmmParams) = size(params.A, 2)

"""
    get_num_constrs(params::AdmmParams)

Return the number of constraints in the problem that `params` was initialized with.
"""
get_num_constrs(params::AdmmParams) = size(params.A, 1)

"""
    Settings

Settings for the ADMM algorithm.

The fields represent:
- `ρ`: Augmented Lagrangian parameter for constraints (must be positive)
- `σ`: Augmented Lagrangian parameter for variables (must be positive)
- `α`:  Over-relaxation parameter (must be ∈ [0, 2], typically ∈ [1.0, 1.8])
- `max_iters`: Maximum allowed ADMM iterations (must be positive)
- `ϵ`: Threshold for determining convergence (must be nonnegative)
- `term_cond_freq`: Frequency in iterations with which to check for termination
- `compute_stats`: Whether to compute objectives and residuals on every iteration
- `term_cond_type`:
    - `1` uses the `ConvergeTermCache`
    - `2` uses the `FirstVarsTermCache`
- `obj_tol`: Minimum amount by which the objective needs to improve to be considered and improvement
- `res_tol`: Maximum allowable residual for an iterate to be considered an improvement
- `non_improvement_iters`: Maximum allowed non-improvement iterations before termination
"""
@with_kw mutable struct Settings
    ρ::Vector{Float64}
    σ::Vector{Float64}
    α::Float64 = 1.0
    max_iters::Int64 = 1000
    ϵ::Float64 = 1e-4
    term_cond_freq::Int = 10
    compute_stats::Bool = false
    term_cond_type::Int64 = 1
    obj_tol::Float64 = 1e-5
    res_tol::Float64 = 3e-4
    non_improvement_iters::Int64 = 50

    @assert (α > 0) & (α < 2)
    @assert all(ρ .> 0)
    @assert all(σ .> 0)
    @assert max_iters > 0
    @assert ϵ >= 0
    @assert obj_tol > 0
    @assert res_tol > 0
    @assert non_improvement_iters > 0
end

"""
    Settings(n:: Int, m:: Int)

Construct ADMM Settings with default values given the number of variables `n` and constraints `m`
"""
Settings(n::Int, m::Int) = Settings(; ρ=ones(m), σ=ones(n))

"""
    Settings(n:: Int, m:: Int)

Construct ADMM Settings with default values in accordance with set of ADMM Params `params`
"""
function Settings(params::AdmmParams)
    n = get_num_vars(params::AdmmParams)
    m = get_num_constrs(params::AdmmParams)
    return Settings(n::Int, m::Int)
end

"""
    get_num_vars(settings::Settings)

Return the number of variables in the problem that `settings` was initialized with.
"""
get_num_vars(settings::Settings) = length(settings.σ)

"""
    get_num_constrs(settings::Settings)

Return the number of constraints in the problem that `settings` was initialized with.
"""
get_num_constrs(settings::Settings) = length(settings.ρ)

"""
    Stats

Stats from the ADMM algorithm per iteration

The fields represent:
- `obj`: Objective values
- `res`: Residuals
- `iters`: Number of iterations completed
"""
@with_kw mutable struct Stats
    obj::Vector{Float64}
    res::Vector{Float64}
    iters::Int64

    function Stats(max_iters::Int64)
        obj = fill(NaN, max_iters)
        res = fill(NaN, max_iters)
        return new(obj, res, 0)
    end
end

Stats(settings::Settings) = Stats(settings.max_iters)

"""
    Vars

Store variables that need to be updated on each iteration of ADMM.

The fields represent:
- `x`: Previous value of x
- `z`: Previous value of z
- `w`: Previous value of w
- `y`: Previous value of y
- `xt`: Auxiliary variable
- `zt`: Auxiliary variable

# Notes:
* We keep the previous values of `x`, `w`, `z`, `y` in order to check for convergence.
* The actual problem we are solving here is detailed in [`SeparableOptimization.factorize_kkt`](@ref)
* `w` and `y` are used when minimizing the augmented Lagrangian of the original problem (and consequently
  in ADMM's update rules)
"""
@with_kw mutable struct Vars
    x::Vector{Float64}
    z::Vector{Float64}
    w::Vector{Float64}
    y::Vector{Float64}
    xt::Vector{Float64}
    zt::Vector{Float64}

    @assert length(x) == length(xt)
    @assert length(x) == length(w)
    @assert length(z) == length(zt)
    @assert length(z) == length(y)
end

"""
    Vars(x::Vector{Float64}, z::Vector{Float64})

Construct ADMM Vars with the appropriate shapes given `x` and `z`
"""
function Vars(x::Vector{Float64}, z::Vector{Float64})
    n = length(x)
    m = length(z)
    return Vars(x, z, zeros(n), zeros(m), x, z)
end
"""
    Vars(n::Int, m::Int)

Construct ADMM Vars given the number of variables `n` and constraints `m`
"""
Vars(n::Int, m::Int) = Vars(zeros(n), zeros(m))
"""
    Vars(params::AdmmParams)

Construct ADMM Vars with the appropriate shapes given AdmmParams `params`
"""
function Vars(params::AdmmParams)
    n = get_num_vars(params)
    m = get_num_constrs(params)
    return Vars(n, m)
end

"""
    get_num_vars(vars::Vars)

Return the number of variables in the problem that `vars` was initialized with.
"""
get_num_vars(vars::Vars) = length(vars.x)

"""
    get_num_constrs(vars::Vars)

Return the number of constraints in the problem that `vars` was initialized with.
"""
get_num_constrs(vars::Vars) = length(vars.z)

"""
    assert_valid(params::AdmmParams, settings::Settings, vars::Vars)

Assert that the number of variables and constraints are consistent across a `params`-`settings`-`vars` combination
"""
function assert_valid(params::AdmmParams, settings::Settings, vars::Vars)
    n = get_num_vars(params)
    m = get_num_constrs(params)
    @assert n == get_num_vars(settings)
    @assert n == get_num_vars(vars)
    @assert m == get_num_constrs(settings)
    @assert m == get_num_constrs(vars)
end
