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

@testset "Settings" begin
    # should not error
    Settings(; ρ=[1], σ=[1])

    # proper consntruction
    @test_throws AssertionError Settings(; ρ=[1], σ=[1], α=3)
    @test_throws AssertionError Settings(; ρ=[-1], σ=[1])
    @test_throws AssertionError Settings(; ρ=[1], σ=[-1])
    @test_throws AssertionError Settings(; ρ=[1], σ=[1], max_iters=-1)
    @test_throws AssertionError Settings(; ρ=[1], σ=[1], obj_tol=-1)
    @test_throws AssertionError Settings(; ρ=[1], σ=[1], res_tol=-1)
    @test_throws AssertionError Settings(; ρ=[1], σ=[1], non_improvement_iters=-1)

    s1 = Settings(3, 3)
    s2 = Settings(; ρ=[1, 1, 1], σ=[1, 1, 1])

    @test s1.ρ == s2.ρ
    @test s1.σ == s2.σ
end

@testset "Params" begin
    # should not error
    P = [1.0 0.0; 0.0 2.0]
    q = [1.0, 2.0]
    A = [1.0 2.0; 0.0 1.0; 1.0 2.0]
    b = [1.0, 2.0, 3.0]
    g = repeat([zero(PiecewiseQuadratic)], 2)
    AdmmParams(P, q, A, b, g)

    @test_throws AssertionError AdmmParams([1.0 0.0 3.0; 0.0 2.0 2.0], q, A, b, g)
    @test_throws AssertionError AdmmParams(P, [1.0, 2.0, 3.0], A, b, g)
    @test_throws AssertionError AdmmParams(P, q, A, [1.0, 2.0], g)
    @test_throws AssertionError AdmmParams(P, q, A, b,
                                           repeat([zero(PiecewiseQuadratic)], 3))
end

@testset "Vars" begin
    # should not error
    Vars(5, 1)
    Vars(zeros(5), zeros(1))
    Vars(zeros(5), zeros(1), zeros(5), zeros(1), zeros(5), zeros(1))

    @test_throws AssertionError Vars(zeros(5), zeros(2), zeros(5), zeros(1), zeros(5),
                                     zeros(1))
    @test_throws AssertionError Vars(zeros(5), zeros(1), zeros(2), zeros(1), zeros(5),
                                     zeros(1))
    @test_throws AssertionError Vars(zeros(5), zeros(1), zeros(5), zeros(2), zeros(5),
                                     zeros(1))
    @test_throws AssertionError Vars(zeros(5), zeros(1), zeros(5), zeros(1), zeros(2),
                                     zeros(1))
    @test_throws AssertionError Vars(zeros(5), zeros(1), zeros(5), zeros(1), zeros(5),
                                     zeros(2))
end

@testset "Validity" begin
    P = [1.0 0.0; 0.0 2.0]
    q = [1.0, 2.0]
    A = [1.0 2.0; 0.0 1.0; 1.0 2.0]
    b = [1.0, 2.0, 3.0]
    g = repeat([zero(PiecewiseQuadratic)], 2)
    params = AdmmParams(P, q, A, b, g)

    vars = Vars(2, 3)

    settings = Settings(; ρ=[1.0, 2.0, 3.0], σ=[1.0, 2.0])

    # should not error
    assert_valid(params, settings, vars)
end
