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

using SeparableOptimization
using PiecewiseQuadratics
using LinearAlgebra
using SparseArrays
using Documenter
using Test

@testset "Doctests" begin
    DocMeta.setdocmeta!(SeparableOptimization, :DocTestSetup,
                        :(using SeparableOptimization, PiecewiseQuadratics); recursive=true)
    doctest(SeparableOptimization)
end

@testset "ADMM Structs" begin
    include("test_admm_structs.jl")
end

@testset "ADMM" begin
    include("test_admm.jl")
end

@testset "Examples" begin
    include("test_examples.jl")
end
