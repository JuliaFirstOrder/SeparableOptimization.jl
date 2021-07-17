# LCSO.jl

**LCSO.jl** is a [Julia](http://julialang.org) package that solves **L**inearly **C**onstrained **S**eparable **O**ptimization Problems.


The package currently supports quadratic-plus-separable problems of the form

```math
\begin{array}{ll}
\text{minimize} &
\displaystyle \frac12 x^T P x + q^T x + \sum_{i=1}^n g_i(x_i) \\
\text{subject to} & Ax = b.
\end{array}
```

where:
* ``x \in \bf R^n`` is the decision variable
* ``A \in \bf R^{m\times n}`` is a sparse matrix
* ``b \in \bf R^m``
* ``P \in \bf S_{+}^n`` is symmetric positive semidefinite matrix
* ``q \in \bf R^n``
* ``g_i: \bf R \to \bf R`` is a piecewise quadratic function, as specified via [PiecewiseQuadratics.jl](https://github.com/JuliaFirstOrder/PiecewiseQuadratics.jl).


The algorithm used is the alternating direction method of multipliers (ADMM).  This method reaches moderate accuracy very quickly, but often requires some tuning, which may need to be done by hand.  This package is therefore best used by someone looking to solve a family of similar optimization problems with excellent performance, even when the function $g_i$ is very complicated.


## Contents

```@contents
Pages = ["index.md", "api.md"]
Depth = 2
```


## Installation
Use Julia's builtin package manager [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) to install.
From a Julia REPL:
```Julia
] add LCSO
```

## Example
Let's use LCSO to solve an example problem.

```@example 1
using LCSO
using PiecewiseQuadratics
using LinearAlgebra

n = 4 # num features
m = 2 # num constraints

# construct problem data (ensuring the problem is feasible)
x0 = rand(n)
A = rand(m, n)
b = A * x0
X = rand(n,n)
P = X'X  # ensure P is positive definite
@assert isposdef(P)
q = rand(n)

# x1 has to be in union([-1, 2], [2.5, 3.5]) and has a quadratic penalty if
# it lies in [-1, 2] and a linear penalty if it lies in [2.5, 3.5]
g1 = PiecewiseQuadratic([BoundedQuadratic(-1, 2, 1, 0, 0),
                        BoundedQuadratic(2.5, 3.5, 0, 1, 0)])
# x2 has to be between -20 and 10
g2 = indicator(-20, 10)

# x3 has to be between -5 and 10
g3 = indicator(-5, 10)

# x4 has to be exactly 1.2318
g4 = indicator(1.2318, 1.2318);

g = [g1,g2,g3,g4]

# solve
params = AdmmParams(P, q, A, b, g)
settings = Settings(; ρ=ones(m), σ=ones(n), compute_stats=true)

vars, stats = optimize(params, settings)

println("optimal x: ", vars.x)
println("final obj: ", stats.obj[stats.iters])
println("final res: ", stats.res[stats.iters])
```

## Authors
This package and [PiecewiseQuadratics.jl](https://github.com/JuliaFirstOrder/PiecewiseQuadratics.jl) were originally developed by [Nicholas Moehle](https://www.nicholasmoehle.com/), [Ellis Brown](http://ellisbrown.github.io), and [Mykel Kochenderfer](https://mykel.kochenderfer.com/) at BlackRock AI Labs.  They were developed to produce the results in the following paper: [arXiv:2103.05455](https://arxiv.org/abs/2103.05455).

## Reference
```@contents
Pages = ["api.md"]
Depth = 3
```
