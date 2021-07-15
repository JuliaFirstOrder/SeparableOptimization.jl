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

push!(LOAD_PATH, "../src/")

using Documenter
using LCSO
using PiecewiseQuadratics


DocMeta.setdocmeta!(LCSO, :DocTestSetup, :(using LCSO, PiecewiseQuadratics); recursive=true)


makedocs(
    sitename = "LCSO.jl",
    authors = "Nick Moehle, Ellis Brown, Mykel Kochenderfer",
    repo = "github.com/JuliaFirstOrder/LCSO.jl.git",
    modules = [LCSO],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "index.md",
        "api.md"
    ]
)

deploydocs(;
    repo = "github.com/JuliaFirstOrder/LCSO.jl.git",
    devbranch = "main",
    push_preview = true,
)
