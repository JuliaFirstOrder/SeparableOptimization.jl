# API Documentation

Docstrings for LCSO.jl interface members can be [accessed through Julia's built-in documentation system](https://docs.julialang.org/en/v1/manual/documentation/#Accessing-Documentation) or in the list below.

```@meta
CurrentModule = LCSO
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
LCSO.factorize_kkt
LCSO.solve_kkt!
LCSO.admm_step!
LCSO.compute_stats
```

## Prox cache

```@docs
LCSO.ProxCache
LCSO.prox_step
```

## Convergence

```@docs
LCSO.get_term_cond_cache
LCSO.check_term_conds!
LCSO.ConvergeTermCache
LCSO.FirstVarsTermCache
LCSO.get_obj_and_res_from_first_vars
```
