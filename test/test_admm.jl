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

@testset "Identity constrained QP" begin
    P = Diagonal(trues(4))  # === I(4) in 1.2+
    q = ones(4)
    A = spzeros(1, 4)
    b = [1.0]
    g = repeat([zero(PiecewiseQuadratic)], 4)
    params = AdmmParams(P, q, A, b, g)

    vars, _ = optimize(params)

    @test norm(vars.x - (-1) * ones(4)) < 1e-5
end

@testset "Diag constrained QP" begin
    P = spdiagm(0 => 1:4)
    q = Vector{Float64}(1:4)
    A = spzeros(1, 4) .+ 1
    b = [1.0]
    g = repeat([zero(PiecewiseQuadratic)], 4)
    params = AdmmParams(P, q, A, b, g)

    vars, _ = optimize(params)

    @test norm(vars.x - [1.4, 0.2, -0.2, -0.4]) < 1e-5
end

@testset "Test with ρ and σ" begin
    P = [65.0 76.0; 76.0 89.0]
    q = [7.0, 3.0]
    A = [7.0 8.0; 5.0 4.0]
    b = [8.0, 8.0]
    g = repeat([zero(PiecewiseQuadratic)], 2)
    params = AdmmParams(P, q, A, b, g)

    settings = Settings(; ρ=[0.2, 2.0], σ=[2.0, 2.0], ϵ=1e-5)

    vars, _ = optimize(params, settings)

    @test norm(vars.x - [2.66666651, -1.33333319]) < 1e-5
end

@testset "Test with g_i intervals containing feasible point" begin
    P = [5.0 3.0; 3.0 2.0]
    q = [1.0, 2.0]
    A = [1.0 1.0]
    b = [-3.0]
    g = repeat([zero(PiecewiseQuadratic)], 2)
    params = AdmmParams(P, q, A, b, g)

    settings = Settings(; ρ=[0.2], σ=[2.0, 2.0], ϵ=1e-6)

    vars, _ = optimize(params, settings)

    @test norm(vars.x - [4.0, -7.0]) < 1e-5
end

@testset "Test with g_i quadratic" begin
    P = [5.0 3.0; 3.0 2.0]
    q = [1.0, 2.0]
    A = [1.0 1.0; 1.0 2.0]
    b = [1.0, 1.0]
    g = [PiecewiseQuadratic(BoundedQuadratic(1, 0, 0)), zero(PiecewiseQuadratic)]
    params = AdmmParams(P, q, A, b, g)

    settings = Settings(; ρ=[5.0, 5.0], σ=[5.0, 2.0], ϵ=1e-6)

    vars, _ = optimize(params, settings)

    @test norm(vars.x - [1.0, 0.0]) < 1e-5
end
