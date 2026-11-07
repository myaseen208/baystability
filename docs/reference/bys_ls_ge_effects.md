# Genotype-by-environment interaction effects

`bys_ls_ge_effects()` computes the interaction effects (the cell means
centred on genotype and environment means) and their singular value
decomposition, the least-squares estimate of the multiplicative part of
the AMMI model.

## Usage

``` r
bys_ls_ge_effects(.data, .y, .gen, .env)

# Default S3 method
bys_ls_ge_effects(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bys_ge_eff", "bys_result", "list")`) with `ge_effects`
(the \\G \times E\\ matrix \\\mathbf{Z}\\) and `ge_svd`
([`svd()`](https://rdrr.io/r/base/svd.html) of it: `d`, `u`, `v`).

## Details

\$\$\mathbf{Z} = (z\_{ij}), \quad z\_{ij} = \bar{y}\_{ij} -
\bar{y}\_{i\cdot} - \bar{y}\_{\cdot j} + \bar{y}, \qquad \mathbf{Z} =
\sum_k \lambda_k\boldsymbol{\alpha}\_k\boldsymbol{\gamma}\_k^\top,\$\$
with singular values \\\lambda_1 \ge \lambda_2 \ge \cdots\\ and
orthonormal singular vectors orthogonal to \\\mathbf{1}\\.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

Gauch, H. G. (1992). *Statistical Analysis of Regional Yield Trials:
AMMI Analysis of Factorial Designs*. Elsevier, Amsterdam.

## See also

[`bys_ls_ammi`](https://myaseen208.com/baystability/reference/bys_ls_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
fit <-
  bys_ls_ge_effects(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
fit$ge_svd$d
#>  [1] 5.204176e+03 3.885737e+03 3.740472e+03 2.777343e+03 2.704315e+03
#>  [6] 2.106191e+03 1.881744e+03 1.874540e+03 1.362975e+03 9.736497e+02
#> [11] 8.275602e+02 5.623826e-12
```
