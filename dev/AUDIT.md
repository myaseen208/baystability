# Audit of baystability 0.2.0

Run on cultivo2009 with cultivo2008 as prior, with the code of 0.2.0 (rstiefel replaced by an equivalent von
Mises-Fisher draw where it could not be installed).

| # | Finding | Evidence | Fixed in 1.0.0 |
|---|---|---|---|
| 1 | `bayes_ammi()` reads `n0`, `c0`, `mu0`, `tau0`, `tao0`, `delta0`, `lambdas0`, `alphas0`, `gammas0`, `ge_means0` from the user's workspace | "object 'n0' not found" in a clean session; the help example that creates them is `\dontrun` | `.prior_data` argument |
| 2 | Singular values never sampled | SD of every `lambdas1` column over the draws = 0 | truncated normal draw every iteration (eq. 2.15) |
| 3 | `tau1` is a precision (about 1e-7) drawn from a residual matrix `E` computed once before the loop | code: `E` outside the loop | sigma^2 from current residuals |
| 4 | Scores drawn without constraints, then orthogonalised (`orthnorm`) | code | vMF on the constrained spheres |
| 5 | All min(G, E) - 1 = 11 terms; with 12 genotypes a column-wise sampler cannot move the genotype scores (one free dimension per column) | R-hat of the genotype scores > 1e13 with column moves only | exact pair moves (von Mises angle) |
| 6 | `delta1` columns named `delta1..delta12` for 25 environments | `colnames(delta1) <- paste0("delta", 1:g)` | `delta1..delta25` |
| 7 | `cultivo2009` help page describes 1320 rows and 6 columns | data: 900 x 4 | corrected |
| 8 | Paper reports two replicates; data have three | Table B.1 = average of the two years' least-squares effects either way | documented |

Helpers (`ge_mean`, `ge_var`, `g_eff`, `e_eff`, `ge_eff`, `ge_ammi`, `matrix_k`, `orthnorm`): identical results
(dev/legacy/helpers.rds, test-helpers.R); `ge_ammi()` now computes the error variance in closed form, equal to the
REML fit of `ge_model()`.
