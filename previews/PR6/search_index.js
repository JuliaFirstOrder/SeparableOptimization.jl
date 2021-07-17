var documenterSearchIndex = {"docs":
[{"location":"api/#API-Documentation","page":"API Documentation","title":"API Documentation","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"Docstrings for SeparableOptimization.jl interface members can be accessed through Julia's built-in documentation system or in the list below.","category":"page"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"CurrentModule = SeparableOptimization","category":"page"},{"location":"api/#Contents","page":"API Documentation","title":"Contents","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"Pages = [\"api.md\"]\nDepth = 3","category":"page"},{"location":"api/#Index","page":"API Documentation","title":"Index","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Types","page":"API Documentation","title":"Types","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"AdmmParams\nSettings\nStats\nVars","category":"page"},{"location":"api/#SeparableOptimization.AdmmParams","page":"API Documentation","title":"SeparableOptimization.AdmmParams","text":"AdmmParams\n\nParameters to the ADMM algorithm.\n\nThe fields represent:\n\nP: A sparse nxn positive semidefinite matrix\nq: an n-vector\nA: an mxn matrix\nb: an m-vector\ng: a vector of PiecewiseQuadratic functions\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.Settings","page":"API Documentation","title":"SeparableOptimization.Settings","text":"Settings\n\nSettings for the ADMM algorithm.\n\nThe fields represent:\n\nρ: Augmented Lagrangian parameter for constraints (must be positive)\nσ: Augmented Lagrangian parameter for variables (must be positive)\nα:  Over-relaxation parameter (must be ∈ [0, 2], typically ∈ [1.0, 1.8])\nmax_iters: Maximum allowed ADMM iterations (must be positive)\nϵ: Threshold for determining convergence (must be nonnegative)\nterm_cond_freq: Frequency in iterations with which to check for termination\ncompute_stats: Whether to compute objectives and residuals on every iteration\nterm_cond_type:\n1 uses the ConvergeTermCache\n2 uses the FirstVarsTermCache\nobj_tol: Minimum amount by which the objective needs to improve to be considered and improvement\nres_tol: Maximum allowable residual for an iterate to be considered an improvement\nnon_improvement_iters: Maximum allowed non-improvement iterations before termination\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.Stats","page":"API Documentation","title":"SeparableOptimization.Stats","text":"Stats\n\nStats from the ADMM algorithm per iteration\n\nThe fields represent:\n\nobj: Objective values\nres: Residuals\niters: Number of iterations completed\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.Vars","page":"API Documentation","title":"SeparableOptimization.Vars","text":"Vars\n\nStore variables that need to be updated on each iteration of ADMM.\n\nThe fields represent:\n\nx: Previous value of x\nz: Previous value of z\nw: Previous value of w\ny: Previous value of y\nxt: Auxiliary variable\nzt: Auxiliary variable\n\nNotes:\n\nWe keep the previous values of x, w, z, y in order to check for convergence.\nThe actual problem we are solving here is detailed in SeparableOptimization.factorize_kkt\nw and y are used when minimizing the augmented Lagrangian of the original problem (and consequently in ADMM's update rules)\n\n\n\n\n\n","category":"type"},{"location":"api/#Functions","page":"API Documentation","title":"Functions","text":"","category":"section"},{"location":"api/#Utility","page":"API Documentation","title":"Utility","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"get_num_vars\nget_num_constrs\nassert_valid","category":"page"},{"location":"api/#SeparableOptimization.get_num_vars","page":"API Documentation","title":"SeparableOptimization.get_num_vars","text":"get_num_vars(params::AdmmParams)\n\nReturn the number of variables in the problem that params was initialized with.\n\n\n\n\n\nget_num_vars(settings::Settings)\n\nReturn the number of variables in the problem that settings was initialized with.\n\n\n\n\n\nget_num_vars(vars::Vars)\n\nReturn the number of variables in the problem that vars was initialized with.\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.get_num_constrs","page":"API Documentation","title":"SeparableOptimization.get_num_constrs","text":"get_num_constrs(params::AdmmParams)\n\nReturn the number of constraints in the problem that params was initialized with.\n\n\n\n\n\nget_num_constrs(settings::Settings)\n\nReturn the number of constraints in the problem that settings was initialized with.\n\n\n\n\n\nget_num_constrs(vars::Vars)\n\nReturn the number of constraints in the problem that vars was initialized with.\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.assert_valid","page":"API Documentation","title":"SeparableOptimization.assert_valid","text":"assert_valid(params::AdmmParams, settings::Settings, vars::Vars)\n\nAssert that the number of variables and constraints are consistent across a params-settings-vars combination\n\n\n\n\n\n","category":"function"},{"location":"api/#Optimization","page":"API Documentation","title":"Optimization","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"optimize","category":"page"},{"location":"api/#SeparableOptimization.optimize","page":"API Documentation","title":"SeparableOptimization.optimize","text":"optimize(params::AdmmParams[, settings::Settings[, vars::Vars]])\n\nCarries out the ADMM algorithm by calling admm_step until termination conditions are met.\n\n\n\n\n\n","category":"function"},{"location":"api/#Internals","page":"API Documentation","title":"Internals","text":"","category":"section"},{"location":"api/#ADMM","page":"API Documentation","title":"ADMM","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"SeparableOptimization.factorize_kkt\nSeparableOptimization.solve_kkt!\nSeparableOptimization.admm_step!\nSeparableOptimization.compute_stats","category":"page"},{"location":"api/#SeparableOptimization.factorize_kkt","page":"API Documentation","title":"SeparableOptimization.factorize_kkt","text":"factorize_kkt(params::AdmmParams, settings::Settings)\n\nFactorize the KKT matrix representing the optimzality equations of the provided problem.\n\nExtended help\n\nFactorize the coefficient matrix of that KKT system representing the optimality equations that will be used (and reused) for the constrained least squares problem:\n\nbeginarrayll\ntextminimize  displaystyle\nfrac12 x^T P x  + q^T x + sum_i=1^n g_i(x_i) \ntextsubject to  Ax = b\nendarray\n\nWe can formulate this equivalently as\n\nbeginarrayll\ntextminimize  displaystyle\nfrac12 tilde x^T P tilde x  + q^T tilde x + mathcalI_mathcalA(tilde x tilde z)\n        + g(x) + mathcalI_mathcal B(z) \ntextsubject to  tilde x = x \n                   tilde z = z\nendarray\n\nwhere mathcal I_mathcal A is the indicator function over the affine constraints of the original form, i.e.,\n\nmathcalI_mathcal A(x z) = begincases\n0  Ax = z \ninfty  rm otherwise\nendcases\n\nand mathcal I_mathcal B is the indicator function over the singleton set of vectors equal to b, i.e.,\n\nmathcalI_mathcal B(x) = begincases\n0  x = b \ninfty  rm otherwise\nendcases\n\nWe will denote the objective in this forumulation as f(x tilde x z tilde z).\n\nThe augmented Lagrangian for this problem is\n\nL_SR (x tilde x z tilde z w y) =\nf(x tilde x z tilde z)\n+ frac12  tilde x - x + S^-1 w _S^2\n+ frac12  tilde z - z + R^-1 y _R^2\n\nwhere f is the previously mentioned objective function. The diagonal, positive definite matrices Sin bf S++^l and Rin bf S_++^m are algorithm parameters. We define cdot _S is the S-norm, i.e., x_S = sqrtx^T S x.\n\nSolving for the optimality conditions, we can find optimal tilde x and tilde z by solving the linear system\n\nbeginbmatrix P + S  A^T  A  -R^-1 endbmatrix\nbeginbmatrix tilde x  nu endbmatrix\n=\nbeginbmatrix Sx - q - w  z - R^-1 y endbmatrix\n\nwhere tilde z can be recovered as tilde z = z + R^-1(nu - y)\n\nIn this function, we factor the coefficient matrix on the left hand side so that we repeatedly solve the linear system efficiently.\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.solve_kkt!","page":"API Documentation","title":"SeparableOptimization.solve_kkt!","text":"solve_kkt!(vars::Vars, kkt::KKT, params::AdmmParams, settings::Settings)\n\nSolve the KKT system set up in the description of factorize_kkt for a particular right hand side.\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.admm_step!","page":"API Documentation","title":"SeparableOptimization.admm_step!","text":"admm_step!(vars::Vars, params::AdmmParams, kkt::KKT, pc::ProxCache, settings::Settings, iter::Int64; polish::Bool=false)\n\nCarry out a single ADMM iteration.\n\nSteps:\n\nSolve the KKT system to carry out the xt and zt updates.\nUse a mixture of the result of (1) and the previous x iterate to solve n univariate subproblems by minimizing the Lagrangian in the description of factorize_kkt w.r.t. x. (Requires evaluation of a proximal operator.)\nUpdate w and y using the residuals associated with the constraints x = xt and z = zt.\n\nNote:\n\nNote that for the entire run of the algorithm, the z update is just z = b, so we don't need\n\nto do anything here.\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.compute_stats","page":"API Documentation","title":"SeparableOptimization.compute_stats","text":"compute_stats(vars::Vars, params::AdmmParams, stats::Stats, iter::Int64)\n\nUpdate the objective, residual, and iteration fields of the Stats object.\n\n\n\n\n\n","category":"function"},{"location":"api/#Prox-cache","page":"API Documentation","title":"Prox cache","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"SeparableOptimization.ProxCache\nSeparableOptimization.prox_step","category":"page"},{"location":"api/#SeparableOptimization.ProxCache","page":"API Documentation","title":"SeparableOptimization.ProxCache","text":"ProxCache\n\nStore the proximal operators corresponding to each element of a vector of PiecewiseQuadratics\n\nNote\n\nThe proximal operator of f, rho, denoted\n\ntextprox_frho(u) =\nargmin_x in textdom(f) f(x) + fracrho2  x - u _2^2\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.prox_step","page":"API Documentation","title":"SeparableOptimization.prox_step","text":"prox_step(pc::ProxCache, params::AdmmParams, settings::Settings, u::Vector{Float64})\n\nReturn the next x by using the proximal operators of the separable function pieces\n\n\n\n\n\n","category":"function"},{"location":"api/#Convergence","page":"API Documentation","title":"Convergence","text":"","category":"section"},{"location":"api/","page":"API Documentation","title":"API Documentation","text":"SeparableOptimization.get_term_cond_cache\nSeparableOptimization.check_term_conds!\nSeparableOptimization.ConvergeTermCache\nSeparableOptimization.FirstVarsTermCache\nSeparableOptimization.get_obj_and_res_from_first_vars","category":"page"},{"location":"api/#SeparableOptimization.get_term_cond_cache","page":"API Documentation","title":"SeparableOptimization.get_term_cond_cache","text":"get_term_cond_cache(params::AdmmParams, settings::Settings)\n\nReturn the TermCache corresponding to the term_cond_type specification in the ADMM Settings settings\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.check_term_conds!","page":"API Documentation","title":"SeparableOptimization.check_term_conds!","text":"check_term_conds!(tc::ConvergeTermCache, vars::Vars, params::AdmmParams, settings::Settings, iter::Int)\n\nCheck termination conditions\n\n\n\n\n\ncheck_term_conds!(tc::FirstVarsTermCache, vars::Vars, params::AdmmParams, settings::Settings, iter::Int)\n\nCheck termination conditions\n\n\n\n\n\n","category":"function"},{"location":"api/#SeparableOptimization.ConvergeTermCache","page":"API Documentation","title":"SeparableOptimization.ConvergeTermCache","text":"ConvergeTermCache\n\nCache state to helps determine whether or not to terminate\n\nThe fields represent:\n\nx_last: Previous x\nz_last: Previous z\nw_last: Previous w\ny_last: Previous y\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.FirstVarsTermCache","page":"API Documentation","title":"SeparableOptimization.FirstVarsTermCache","text":"FirstVarsTermCache\n\nCache state to helps determine whether or not to terminate\n\nThe fields represent:\n\nobj_best: The best encountered objective value\nres_best: The best encountered residual value\nx_best: The best encountered x iterate\nz_best: The best encountered z iterate\nw_best: The best encountered w (dual) iterate\ny_best: The best encountered y (dual) iterate\nn1: n - m\nlb: Lower bound\nub: Upper bound\nA1: first n1 columns of A\nA2: remaining columns of A\nnot_improved_count: The number of iterations without improvement seeing improvement\nnot_improved_count_req: The max allowed iterations without improvement\n\n\n\n\n\n","category":"type"},{"location":"api/#SeparableOptimization.get_obj_and_res_from_first_vars","page":"API Documentation","title":"SeparableOptimization.get_obj_and_res_from_first_vars","text":"get_obj_and_res_from_first_vars(tc::FirstVarsTermCache, vars::Vars, params::AdmmParams)\n\nGet objective and residual from the FirstVarsTermCache tc\n\n\n\n\n\n","category":"function"},{"location":"#SeparableOptimization.jl","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"SeparableOptimization.jl is a Julia package that solves Linearly Constrained Separable Optimization Problems.","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"The package currently supports quadratic-plus-separable problems of the form","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"beginarrayll\ntextminimize \ndisplaystyle frac12 x^T P x + q^T x + sum_i=1^n g_i(x_i) \ntextsubject to  Ax = b\nendarray","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"where:","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"x in bf R^n is the decision variable\nA in bf R^mtimes n is a sparse matrix\nb in bf R^m\nP in bf S_+^n is symmetric positive semidefinite matrix\nq in bf R^n\ng_i bf R to bf R is a piecewise quadratic function, as specified via PiecewiseQuadratics.jl.","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"The algorithm used is the alternating direction method of multipliers (ADMM).  This method reaches moderate accuracy very quickly, but often requires some tuning, which may need to be done by hand.  This package is therefore best used by someone looking to solve a family of similar optimization problems with excellent performance, even when the function g_i is very complicated.","category":"page"},{"location":"#Contents","page":"SeparableOptimization.jl","title":"Contents","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"Pages = [\"index.md\", \"api.md\"]\nDepth = 2","category":"page"},{"location":"#Installation","page":"SeparableOptimization.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"Use Julia's builtin package manager Pkg to install. From a Julia REPL:","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"] add SeparableOptimization","category":"page"},{"location":"#Example","page":"SeparableOptimization.jl","title":"Example","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"Let's use SeparableOptimization to solve an example problem.","category":"page"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"using SeparableOptimization\nusing PiecewiseQuadratics\nusing LinearAlgebra\n\nn = 4 # num features\nm = 2 # num constraints\n\n# construct problem data (ensuring the problem is feasible)\nx0 = rand(n)\nA = rand(m, n)\nb = A * x0\nX = rand(n,n)\nP = X'X  # ensure P is positive definite\n@assert isposdef(P)\nq = rand(n)\n\n# x1 has to be in union([-1, 2], [2.5, 3.5]) and has a quadratic penalty if\n# it lies in [-1, 2] and a linear penalty if it lies in [2.5, 3.5]\ng1 = PiecewiseQuadratic([BoundedQuadratic(-1, 2, 1, 0, 0),\n                        BoundedQuadratic(2.5, 3.5, 0, 1, 0)])\n# x2 has to be between -20 and 10\ng2 = indicator(-20, 10)\n\n# x3 has to be between -5 and 10\ng3 = indicator(-5, 10)\n\n# x4 has to be exactly 1.2318\ng4 = indicator(1.2318, 1.2318);\n\ng = [g1,g2,g3,g4]\n\n# solve\nparams = AdmmParams(P, q, A, b, g)\nsettings = Settings(; ρ=ones(m), σ=ones(n), compute_stats=true)\n\nvars, stats = optimize(params, settings)\n\nprintln(\"optimal x: \", vars.x)\nprintln(\"final obj: \", stats.obj[stats.iters])\nprintln(\"final res: \", stats.res[stats.iters])","category":"page"},{"location":"#Authors","page":"SeparableOptimization.jl","title":"Authors","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"This package and PiecewiseQuadratics.jl were originally developed by Nicholas Moehle, Ellis Brown, and Mykel Kochenderfer at BlackRock AI Labs.  They were developed to produce the results in the following paper: arXiv:2103.05455.","category":"page"},{"location":"#Reference","page":"SeparableOptimization.jl","title":"Reference","text":"","category":"section"},{"location":"","page":"SeparableOptimization.jl","title":"SeparableOptimization.jl","text":"Pages = [\"api.md\"]\nDepth = 3","category":"page"}]
}
