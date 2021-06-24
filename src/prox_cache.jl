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

@doc raw"""
    ProxCache

Store the proximal operators corresponding to each element of a vector of `PiecewiseQuadratic`s

# Note
The proximal operator of ``f``, ``\rho``, denoted
```math
\text{prox}_{f,\rho}(u) =
\argmin_{x \in \text{dom}(f)} f(x) + \frac{\rho}{2}  \|x - u \|_2^2
```
"""
struct ProxCache
    f_plus_quad::Vector{PiecewiseQuadratic}
end

"""
    ProxCache(params::AdmmParams, settings::Settings)

Construct a ProxCache given AdmmParams `params` and ADMM Settings `settings`
"""
function ProxCache(params::AdmmParams, settings::Settings)
    f_plus_quad = Vector{PiecewiseQuadratic}(undef, length(params.g))
    work = PwqSumWorkspace(2)
    for i in 1:length(params.g)
        σi = settings.σ[i]
        cm = _sum([params.g[i], PiecewiseQuadratic(BoundedQuadratic(σi / 2.0, 0.0, 0.0))],
                  work)
        env = simplify(envelope(cm))
        f_plus_quad[i] = env
    end

    return ProxCache(f_plus_quad)
end

"""
    prox_step(pc::ProxCache, params::AdmmParams, settings::Settings, u::Vector{Float64})

Return the next `x` by using the proximal operators of the separable function pieces
"""
function prox_step(pc::ProxCache, params::AdmmParams, settings::Settings,
                   u::Vector{Float64})
    x = zeros(length(u))
    for i in 1:length(u)
        x[i] = prox(pc.f_plus_quad[i], u[i], settings.σ[i]; use_quadratic=false)
    end
    return x
end
