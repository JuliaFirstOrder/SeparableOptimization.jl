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

module SeparableOptimization

using LinearAlgebra
using LDLFactorizations
using SparseArrays
using Parameters
using PiecewiseQuadratics

import Base: show

# Input data structures:
export AdmmParams, Settings, Stats, Vars

# Optimization function:
export optimize, factorize_kkt

# Utility functions:
export get_num_vars, get_num_constrs, assert_valid, show

# Termination condition types:
export CONVERGE_TERM_COND_FLAG, FIRST_VARS_TERM_COND_FLAG

include("admm_structs.jl")
include("term_conds.jl")
include("prox_cache.jl")
include("admm.jl")

end
