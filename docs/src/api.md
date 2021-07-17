# API Documentation

Docstrings for SeparableOptimization.jl interface members can be [accessed through Julia's built-in documentation system](https://docs.julialang.org/en/v1/manual/documentation/#Accessing-Documentation) or in the list below.

```@meta
CurrentModule = SeparableOptimization
```

## Contents

```@contents
Pages = ["api.md"]
Depth = 3
```

## Index

```@index
Pages = ["api.md"]
```


## Types

```@docs
AdmmParams
Settings
Stats
Vars
```

## Functions

### Utility

```@docs
get_num_vars
get_num_constrs
assert_valid
```

### Optimization

```@docs
optimize
```

# Internals

## ADMM

```@docs
SeparableOptimization.factorize_kkt
SeparableOptimization.solve_kkt!
SeparableOptimization.admm_step!
SeparableOptimization.compute_stats
```

## Prox cache

```@docs
SeparableOptimization.ProxCache
SeparableOptimization.prox_step
```

## Convergence

```@docs
SeparableOptimization.get_term_cond_cache
SeparableOptimization.check_term_conds!
SeparableOptimization.ConvergeTermCache
SeparableOptimization.FirstVarsTermCache
SeparableOptimization.get_obj_and_res_from_first_vars
```
