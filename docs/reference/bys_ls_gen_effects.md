# Genotype effects

`bys_ls_gen_effects()` computes the genotype main effects, the
deviations of the genotype means from the grand mean.

## Usage

``` r
bys_ls_gen_effects(.data, .y, .gen, .env)

# Default S3 method
bys_ls_gen_effects(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A \\G \times 1\\ matrix with column `GEffs` (row names: genotypes).

## Details

\$\$\hat{\tau}\_i = \bar{y}\_{i\cdot} - \bar{y}, \qquad
\bar{y}\_{i\cdot} = \frac{1}{E}\sum_j \bar{y}\_{ij},\$\$ the
least-squares estimates of the AMMI model (equal to the deviations of
the plot means when the data are balanced).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`bys_ls_env_effects`](https://myaseen208.com/baystability/reference/bys_ls_env_effects.md),
[`bys_ls_ge_effects`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
bys_ls_gen_effects(
    .data = cultivo2009
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#>          GEffs
#> 2   622.710632
#> 10  125.855653
#> 12   -9.319319
#> 13  117.406089
#> 14   53.786287
#> 18  283.405021
#> 19  293.819680
#> 20   11.972994
#> 23  199.565996
#> 38 -653.521319
#> 40 -278.263429
#> 42 -767.418284
```
