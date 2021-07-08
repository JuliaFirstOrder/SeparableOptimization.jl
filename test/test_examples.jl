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

using Convex
using ECOS
using PiecewiseQuadratics
using LinearAlgebra

@testset "Only quadratic" begin
    P = ones(1, 1)
    q = -2 * ones(1)
    A = zeros(0, 1)
    b = zeros(0)
    g = zeros(PiecewiseQuadratic, 1)

    ρ = zeros(0)
    σ = ones(1)

    params = AdmmParams(P, q, A, b, g)
    settings = Settings(; ρ=ρ, σ=σ)
    vars, _ = optimize(params, settings)

    @test abs(vars.x[1] - 2) < 1e-3
    @test abs(vars.xt[1] - 2) < 1e-3
end

@testset "Only separable" begin
    P = zeros(1, 1)
    q = zeros(1)
    A = zeros(0, 1)
    b = zeros(0)
    g = [PiecewiseQuadratic([BoundedQuadratic(2.0, 3.0, 0.0, -1.0, 10.0),
                             BoundedQuadratic(3.0, 5.0, 0.0, 1.0, 4.0)])]

    ρ = zeros(0)
    σ = ones(1)

    params = AdmmParams(P, q, A, b, g)
    settings = Settings(; ρ=ρ, σ=σ)
    vars, _ = optimize(params, settings)

    display(vars)
    @test abs(vars.x[1] - 3) < 1e-3
    @test abs(vars.xt[1] - 3) < 1e-3
end

@testset "Random quadratic program" begin
    # Original generating code:
    # P = randn(5, 5)
    # P = P'*P
    # A = randn(2, 5)

    P = [1.29585 2.90775 -0.421317 1.52657 0.209435
         2.90775 10.5604 0.298516 4.64058 -2.92845
         -0.421317 0.298516 1.94608 -1.61207 -0.198397
         1.52657 4.64058 -1.61207 4.06981 -2.02592
         0.209435 -2.92845 -0.198397 -2.02592 4.66317]
    q = -10 * ones(5)
    A = [1.33694 -0.745464 -0.0531773 -2.11537 1.22246
         -0.0731486 -1.22006 -0.165136 -0.066768 0.567695]
    b = A * ones(5)
    g = repeat([indicator(0, Inf)], 5)

    x = Variable(5)
    problem = Convex.minimize(0.5 * quadform(x, P) + q' * x, [x >= 0, A * x == b])
    solve!(problem, ECOS.Optimizer)

    ρ = ones(2)
    σ = ones(5)
    params = AdmmParams(P, q, A, b, g)
    settings = Settings(; ρ=ρ, σ=σ, term_cond_freq=1000, max_iters=1000)
    vars, stats = optimize(params, settings)

    @test norm(vars.x - x.value) < 1e-3
    @test norm(vars.xt - x.value) < 1e-3
end

@testset "Random Lasso problem" begin
    n = 2  # num features
    m = 0  # num constraints

    X = [1.33694 -0.745464
         -0.0531773 -2.11537
         1.22246 -0.0731486
         -1.22006 -0.165136
         -0.066768 0.567695]
    y = X * ones(n)
    λ = 1.0

    β = Variable(n)
    problem = Convex.minimize(0.5 * sumsquares(X * β - y) + λ * sum(abs(β)))
    solve!(problem, ECOS.Optimizer)

    ρ = ones(m)
    σ = ones(n)
    P = X' * X
    q = -X' * y
    A = zeros(m, n)
    b = zeros(m)
    g = [PiecewiseQuadratic([BoundedQuadratic(-Inf, 0.0, 0.0, -λ, 0.0),
                             BoundedQuadratic(0.0, Inf, 0.0, λ, 0.0)]) for i in 1:n]
    params = AdmmParams(P, q, A, b, g)
    settings = Settings(; ρ=ρ, σ=σ, term_cond_freq=1000, max_iters=1000)
    vars, _ = optimize(params, settings)

    @test norm(vars.x - β.value) < 1e-3
    @test norm(vars.xt - β.value) < 1e-3

    # first vars term cache
    settings = Settings(; ρ=ρ, σ=σ, term_cond_freq=1000, max_iters=1000,
                        term_cond_type=FIRST_VARS_TERM_COND_FLAG, compute_stats=true)
    vars, stats = optimize(params, settings)

    @test norm(vars.x - β.value) < 1e-3
    @test norm(vars.xt - β.value) < 1e-3
end

@testset "small example" begin
    # NOTE: example from README.md

    # construct problem data (ensuring the problem is feasible)
    x0 = [0.0036165677461501566, 0.9751164348173793, 0.4825907494313493,
          0.47578842808561417]
    A = [0.477388 0.796388 0.509418 0.426954;
         0.772124 0.381946 0.412414 0.815745]
    b = A * x0
    X = [0.444226 0.331317 0.870612 0.576343;
         0.236767 0.392398 0.000545117 0.0703828;
         0.413522 0.893505 0.366298 0.654319;
         0.149137 0.0738117 0.737134 0.403253]
    P = X'X  # ensure P is positive definite
    @assert isposdef(P)
    q = [0.007228400155782522, 0.9979590358695292, 0.6284406861482683, 0.08716574428367818]

    # x1 has to be in union([-1, 2], [2.5, 3.5]) and has a quadratic penalty if
    # it lies in [-1, 2] and a linear penalty if it lies in [2.5, 3.5]
    g1 = PiecewiseQuadratic([BoundedQuadratic(-1, 2, 1, 0, 0),
                             BoundedQuadratic(2.5, 3.5, 0, 1, 0)])
    # x2 has to be between -20 and 10
    g2 = indicator(-20, 10)

    # x3 has to be between -5 and 10
    g3 = indicator(-5, 10)

    # x4 has to be exactly 1.2318
    g4 = indicator(1.2318, 1.2318)

    g = [g1, g2, g3, g4]

    # solve
    params = AdmmParams(P, q, A, b, g)

    # display
    display(params)

    vars, _ = optimize(params)

    expected = [-0.39769, 1.73815, -0.96784, 1.2318]
    @test norm(vars.x - expected) < 1e-4
end
