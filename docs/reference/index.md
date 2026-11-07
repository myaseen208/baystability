# Package index

## Bayesian linear-bilinear models

The Gibbs sampler of Perez-Elizalde et al. (2012) with an earlier trial
as prior, and the biplot with HPD regions.

- [`bys_priors()`](https://myaseen208.com/baystability/reference/bys_priors.md)
  [`print(`*`<bys_priors>`*`)`](https://myaseen208.com/baystability/reference/bys_priors.md)
  : Prior distributions of the Bayesian AMMI model
- [`bayes_ammi()`](https://myaseen208.com/baystability/reference/bayes_ammi.md)
  : Bayesian estimation of linear-bilinear models: AMMI and SREG (GGE)
- [`bayes_gge()`](https://myaseen208.com/baystability/reference/bayes_gge.md)
  : Bayesian GGE (sites regression) model
- [`bayes_ammi_biplot()`](https://myaseen208.com/baystability/reference/bayes_ammi_biplot.md)
  [`bayes_gge_biplot()`](https://myaseen208.com/baystability/reference/bayes_ammi_biplot.md)
  : AMMI and GGE biplots with highest posterior density regions

## Model choice, interpretation and references

- [`bys_compare()`](https://myaseen208.com/baystability/reference/bys_compare.md)
  : Compare Bayesian AMMI and GGE models
- [`as_draws_df(`*`<bys_bayes_ammi>`*`)`](https://myaseen208.com/baystability/reference/bys-draws.md)
  [`as_draws(`*`<bys_bayes_ammi>`*`)`](https://myaseen208.com/baystability/reference/bys-draws.md)
  [`as.mcmc.list(`*`<bys_bayes_ammi>`*`)`](https://myaseen208.com/baystability/reference/bys-draws.md)
  : The draws of a fit for posterior, bayesplot and coda
- [`bys_interpret()`](https://myaseen208.com/baystability/reference/bys_interpret.md)
  [`print(`*`<bys_interpretation>`*`)`](https://myaseen208.com/baystability/reference/bys_interpret.md)
  [`format(`*`<bys_interpretation>`*`)`](https://myaseen208.com/baystability/reference/bys_interpret.md)
  : Interpretation of a Bayesian AMMI analysis
- [`bys_cite()`](https://myaseen208.com/baystability/reference/bys_cite.md)
  : References of a result or of the package

## The application

- [`bys_app()`](https://myaseen208.com/baystability/reference/bys_app.md)
  [`bys_app_object()`](https://myaseen208.com/baystability/reference/bys_app.md)
  : Bayesian AMMI and SREG analysis in a web browser

## Tables and least-squares AMMI

Cell means and variances, main and interaction effects, least-squares
AMMI and the mixed model.

- [`bys_ls_means()`](https://myaseen208.com/baystability/reference/bys_ls_means.md)
  : Genotype-by-environment means
- [`bys_ls_var()`](https://myaseen208.com/baystability/reference/bys_ls_var.md)
  : Genotype-by-environment variances
- [`bys_ls_gen_effects()`](https://myaseen208.com/baystability/reference/bys_ls_gen_effects.md)
  : Genotype effects
- [`bys_ls_env_effects()`](https://myaseen208.com/baystability/reference/bys_ls_env_effects.md)
  : Environment effects
- [`bys_ls_ge_effects()`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md)
  : Genotype-by-environment interaction effects
- [`bys_ls_ammi()`](https://myaseen208.com/baystability/reference/bys_ls_ammi.md)
  : Least-squares AMMI model
- [`bys_ls_model()`](https://myaseen208.com/baystability/reference/bys_ls_model.md)
  : Linear mixed model of a multi-environment trial

## Matrices

- [`bys_contrasts()`](https://myaseen208.com/baystability/reference/bys_contrasts.md)
  : Orthonormal contrast matrix
- [`bys_orthonormal()`](https://myaseen208.com/baystability/reference/bys_orthonormal.md)
  : Gram-Schmidt orthonormalisation

## Data

- [`cultivo2008`](https://myaseen208.com/baystability/reference/cultivo2008.md)
  [`cultivo2009`](https://myaseen208.com/baystability/reference/cultivo2008.md)
  : Maize trials of 2008 and 2009

## Results and utilities

- [`print(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`summary(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`print(`*`<summary.bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`format(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`as.data.table(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`as.data.frame(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  [`plot(`*`<bys_result>`*`)`](https://myaseen208.com/baystability/reference/bys-methods.md)
  : Methods for baystability results
- [`bys_theme()`](https://myaseen208.com/baystability/reference/bys_theme.md)
  : ggplot2 theme of baystability
- [`bys_save()`](https://myaseen208.com/baystability/reference/bys_save.md)
  : Save a figure (vector PDF by default)
- [`bys_notation()`](https://myaseen208.com/baystability/reference/bys_notation.md)
  : Notation of baystability
- [`bys_reference()`](https://myaseen208.com/baystability/reference/bys_reference.md)
  : Published values used by the tests and vignettes
- [`baystability`](https://myaseen208.com/baystability/reference/baystability-package.md)
  [`baystability-package`](https://myaseen208.com/baystability/reference/baystability-package.md)
  : baystability: Bayesian Stability Analysis of Genotype by Environment
  Interaction (GEI)
