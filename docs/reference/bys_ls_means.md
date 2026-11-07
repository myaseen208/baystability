# Genotype-by-environment means

`bys_ls_means()` computes the table of genotype-by-environment cell
means and the grand mean, the data of the AMMI model.

## Usage

``` r
bys_ls_means(.data, .y, .gen, .env)

# Default S3 method
bys_ls_means(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bys_ge_mean", "bys_result", "list")`) with `ge_means`,
the \\G \times E\\ matrix of means (row and column names: genotypes and
environments), and `grand_mean`.

## Details

With \\y\_{ijl}\\ the response of genotype \\i\\ in replicate \\l\\ of
environment \\j\\ and \\n\_{ij}\\ the number of replicates of the cell,
\$\$\bar{y}\_{ij} = \frac{1}{n\_{ij}}\sum_l y\_{ijl}, \qquad \bar{y} =
\frac{1}{GE}\sum\_{ij}\bar{y}\_{ij}.\$\$ Genotypes and environments are
in the order of their factor levels (sorted values when the columns are
not factors).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`bys_ls_var`](https://myaseen208.com/baystability/reference/bys_ls_var.md),
[`bys_ls_ge_effects`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md),
[`bys_ls_ammi`](https://myaseen208.com/baystability/reference/bys_ls_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
bys_ls_means(
    .data = cultivo2009
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#> <bys_ge_mean> Genotype-by-environment means
#> 
#> $ge_means
#>       1    2    3    4    5    6    7    8      9     10   11   12   13   14
#> 2  5091 4206 5346 3968 6048 2240 1805 3465 2170.8 1994.2 2909 1598 5128 4737
#> 10 4738 3009 4636 2527 5971 1402 1389 3432 1703.5  736.6 2263 1573 4285 4427
#> 12 4243 3205 4446 3470 5543 2300 2225 3135 3199.7 2429.8 2869 1534 3347 3420
#> 13 4662 4345 5001 3682 5965 3094 3407 2611 1250.6 1872.1 2828 1348 3379 3403
#> 14 4609 3698 5175 3585 5771 1275 2003 4154 1615.5 1698.0 2222 1549 3747 3553
#> 18 4936 4048 6758 3832 6409 1716 1179 3341  790.7 2390.8 2828 1703 3820 3405
#> 19 4747 4940 5478 4066 5240 1501 1793 3778 2711.4 2675.8 2788 1509 4025 4039
#> 20 4617 3964 5378 3102 5674 1717 2040 3375 2288.0 2001.0 2182 1742 2981 3416
#> 23 4169 3546 4301 2771 6131 2228 2058 3626 1846.2 2771.2 2384 1465 4166 3650
#> 38 3589 2648 4580 2120 5193 1798 1462 2199 3098.8 1262.9 2545 1931 2962 3228
#> 40 4625 3546 4296 3147 6410 1606 1981 3414 1247.7 1391.1 2182 1514 3517 3008
#> 42 3912 2220 4911 1713 5394 1540 2011 2633 2485.2 1950.4 2101 1436 2360 2616
#>      15   16   17   18   19   20   21   22   23   24   25
#> 2  2752 4010 5010 6214 4097 3394 4364 8157 7885 7191 2472
#> 10 1987 4725 3313 6254 5131 2182 3030 7980 7345 6865 2927
#> 12 1750 4852 4242 5207 3148 2141 1939 6440 6878 6218 2267
#> 13 2385 4495 4323 6034 2294 2303 3030 7429 6053 5774 2650
#> 14 2370 4034 3556 5867 3881 1818 3758 8525 5183 5649 2732
#> 18 1651 5076 4566 6212 3645 2545 3475 8516 5609 6031 3285
#> 19 2007 4966 4485 5844 2610 1899 3495 7437 6208 6893 2893
#> 20 1936 3283 4323 5102 3600 2182 4162 7948 5348 5819 2802
#> 23 1651 5430 4606 5696 4301 2263 3596 8037 6319 5367 3294
#> 38 1349 3615 3475 4638 3449 1535 1899 7197 3510 3810 1249
#> 40 1789 3079 3677 5372 2881 1859 2990 7182 5627 5050 2337
#> 42 1005 3838 3475 4026 2295 1818 1576 5630 4253 4221 2078
#> 
#> $grand_mean
#> [1] 3627
```
