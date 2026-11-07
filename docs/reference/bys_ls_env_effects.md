# Environment effects

`bys_ls_env_effects()` computes the environment main effects, the
deviations of the environment means from the grand mean.

## Usage

``` r
bys_ls_env_effects(.data, .y, .gen, .env)

# Default S3 method
bys_ls_env_effects(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

An \\E \times 1\\ matrix with column `EEffs` (row names: environments).

## Details

\$\$\hat{\delta}\_j = \bar{y}\_{\cdot j} - \bar{y}, \qquad
\bar{y}\_{\cdot j} = \frac{1}{G}\sum_i \bar{y}\_{ij}.\$\$

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`bys_ls_gen_effects`](https://myaseen208.com/baystability/reference/bys_ls_gen_effects.md),
[`bys_ls_ge_effects`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
bys_ls_env_effects(
    .data = cultivo2009
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#>          EEffs
#> 1    867.53704
#> 2    -12.65547
#> 3   1398.26550
#> 4   -462.15834
#> 5   2185.12706
#> 6  -1759.08684
#> 7  -1681.29639
#> 8   -363.79748
#> 9  -1593.28134
#> 10 -1696.13893
#> 11 -1118.88245
#> 12 -2052.17221
#> 13    15.77075
#> 14   -52.07042
#> 15 -1741.30888
#> 16   656.17796
#> 17   460.24213
#> 18  1911.54259
#> 19  -182.96801
#> 20 -1465.68379
#> 21  -517.87235
#> 22  3912.49644
#> 23  2224.03219
#> 24  2113.39367
#> 25 -1045.21241
```
