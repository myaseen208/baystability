# Least-squares AMMI model

`bys_ls_ammi()` computes the least-squares estimates of the AMMI model
with all interaction axes: grand mean, genotype and environment effects,
singular values and vectors of the interaction, and the error variance.
They are the starting values of
[`bayes_ammi()`](https://myaseen208.com/baystability/reference/bayes_ammi.md)
and the classical estimates it is compared with.

## Usage

``` r
bys_ls_ammi(.data, .y, .gen, .env, .rep)

# Default S3 method
bys_ls_ammi(.data, .y, .gen, .env, .rep)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env, .rep:

  the columns of the response, genotypes, environments and replicates:
  bare names, strings or vectors.

## Value

A list (class `c("bys_ge_ammi", "bys_result", "list")`) with the
components of the earlier versions:

- `g`, `e`, `Rep`, `k`:

  numbers of genotypes, environments, replicates and axes.

- `ge_means`:

  the result of
  [`bys_ls_means()`](https://myaseen208.com/baystability/reference/bys_ls_means.md).

- `mu`, `sigma2`, `tau`:

  grand mean, error variance and its reciprocal (precision).

- `tao`, `delta`:

  genotype and environment effects
  ([`bys_ls_gen_effects()`](https://myaseen208.com/baystability/reference/bys_ls_gen_effects.md),
  [`bys_ls_env_effects()`](https://myaseen208.com/baystability/reference/bys_ls_env_effects.md)).

- `lambdas`, `alphas`, `gammas`:

  singular values and genotype and environment singular vectors of the
  interaction.

## Details

With the cell means \\\bar{y}\_{ij}\\ of \\B\\ replicates, \\\hat{\mu} =
\bar{y}\\, \\\hat{\tau}\_i = \bar{y}\_{i\cdot} - \bar{y}\\,
\\\hat{\delta}\_j = \bar{y}\_{\cdot j} - \bar{y}\\ and the singular
value decomposition of the interaction effects
([`bys_ls_ge_effects()`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md)),
\\K = \min(G, E) - 1\\ axes. The error variance is the residual variance
of the model with random replicates within environments
([`bys_ls_model()`](https://myaseen208.com/baystability/reference/bys_ls_model.md)),
in closed form for balanced data: \$\$\hat{\sigma}^2 =
\frac{\mathrm{SSW} - \mathrm{SS}\_{B(E)}}{E(B - 1)(G - 1)}\$\$ when the
replicate mean square is not smaller, otherwise (the replicate variance
component estimated as zero) \\\mathrm{SSW}/(GE(B - 1))\\; SSW is the
within-cell and \\\mathrm{SS}\_{B(E)}\\ the replicate-within-environment
sum of squares.

**Change from the earlier versions.** The error variance came from
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html); the closed
form gives the same value for balanced data without needing lme4.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

Gauch, H. G. (1992). *Statistical Analysis of Regional Yield Trials:
AMMI Analysis of Factorial Designs*. Elsevier, Amsterdam.

## See also

[`bayes_ammi`](https://myaseen208.com/baystability/reference/bayes_ammi.md),
[`bys_ls_ge_effects`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
fit <-
  bys_ls_ammi(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
fit$lambdas
#>  [1] 5204.1759 3885.7373 3740.4723 2777.3430 2704.3153 2106.1909 1881.7436
#>  [8] 1874.5398 1362.9746  973.6497  827.5602
fit$sigma2
#> [1] 589051.8
```
