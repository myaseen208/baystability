# Deprecated names of baystability 0.2.0

These functions of baystability 0.2.0 still work in 1.0.0, with a
warning, and will be removed in a later version. Use the new names,
which follow those of the package stability and its Bayesian version
bayesstability:

|  |  |
|----|----|
| old name | new name |
| `ge_mean()` | [`bys_ls_means()`](https://myaseen208.com/baystability/reference/bys_ls_means.md) |
| `ge_var()` | [`bys_ls_var()`](https://myaseen208.com/baystability/reference/bys_ls_var.md) |
| `g_eff()` | [`bys_ls_gen_effects()`](https://myaseen208.com/baystability/reference/bys_ls_gen_effects.md) |
| `e_eff()` | [`bys_ls_env_effects()`](https://myaseen208.com/baystability/reference/bys_ls_env_effects.md) |
| `ge_eff()` | [`bys_ls_ge_effects()`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md) |
| `ge_model()` | [`bys_ls_model()`](https://myaseen208.com/baystability/reference/bys_ls_model.md) |
| `ge_ammi()` | [`bys_ls_ammi()`](https://myaseen208.com/baystability/reference/bys_ls_ammi.md) |
| `matrix_k()` | [`bys_contrasts()`](https://myaseen208.com/baystability/reference/bys_contrasts.md) |
| `orthnorm()` | [`bys_orthonormal()`](https://myaseen208.com/baystability/reference/bys_orthonormal.md) |

The arguments and results are unchanged.

## Usage

``` r
ge_mean(.data, .y, .gen, .env)

# Default S3 method
ge_mean(.data, .y, .gen, .env)

ge_var(.data, .y, .gen, .env)

# Default S3 method
ge_var(.data, .y, .gen, .env)

g_eff(.data, .y, .gen, .env)

# Default S3 method
g_eff(.data, .y, .gen, .env)

e_eff(.data, .y, .gen, .env)

# Default S3 method
e_eff(.data, .y, .gen, .env)

ge_eff(.data, .y, .gen, .env)

# Default S3 method
ge_eff(.data, .y, .gen, .env)

ge_model(.data, .y, .gen, .env, .rep)

# Default S3 method
ge_model(.data, .y, .gen, .env, .rep)

ge_ammi(.data, .y, .gen, .env, .rep)

# Default S3 method
ge_ammi(.data, .y, .gen, .env, .rep)

matrix_k(n)

# Default S3 method
matrix_k(n)

orthnorm(u = NULL, basis = TRUE, norm = TRUE)

# Default S3 method
orthnorm(u = NULL, basis = TRUE, norm = TRUE)
```

## Arguments

- .data, .y, .gen, .env, .rep, n, u, basis, norm:

  the arguments of the earlier version, passed on to the new function.

## Value

The result of the new function.
