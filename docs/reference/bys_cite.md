# References of a result or of the package

`bys_cite()` gives the references of the methods behind a result of the
package, or the references of every method in the package, as text or
BibTeX. `citation("baystability")` cites the package and its main
methods.

## Usage

``` r
bys_cite(x = NULL, style = c("text", "bibtex"))

# Default S3 method
bys_cite(x = NULL, style = c("text", "bibtex"))
```

## Arguments

- x:

  a result of a function of the package, or `NULL` (default) for every
  reference of the package.

- style:

  `"text"` (default) or `"bibtex"`.

## Value

An object of class `bibentry` (printed in the chosen style).

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bys_cite(
    x     = NULL
  , style = "text"
  )
#> Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive main effects and multiplicative interaction model. Crop Science, 51(4), 1458-1469. doi:10.2135/cropsci2010.06.0343
#> 
#> Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear models applied to plant breeding trials with genotype x environment interaction. Journal of Agricultural, Biological, and Environmental Statistics, 17(1), 15-37. doi:10.1007/s13253-011-0063-9
#> 
#> Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of multiplicative interaction in analysis of variance using Kullback-Leibler information. Journal of Statistical Planning and Inference, 84(1-2), 201-219.
#> 
#> Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with applications to multivariate and relational data. Journal of Computational and Graphical Statistics, 18(2), 438-456. doi:10.1198/jcgs.2009.07177
#> 
#> Yan, W., Hunt, L. A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar evaluation and mega-environment investigation based on the GGE biplot. Crop Science, 40(3), 597-605. doi:10.2135/cropsci2000.403597x
#> 
#> Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society B, 64(4), 583-639. doi:10.1111/1467-9868.00353
#> 
#> Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571-3594.
#> 
#> Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413-1432. doi:10.1007/s11222-016-9696-4
#> 
#> Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society A, 182(2), 389-402. doi:10.1111/rssa.12378 
data(cultivo2009)
g <-
  bys_ls_ammi(
      .data = cultivo2009
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
bys_cite(
    x     = g
  , style = "bibtex"
  )
#> @Article{PerezElizalde2012,
#>   title = {A General {Bayesian} Estimation Method of Linear-Bilinear Models Applied to Plant Breeding Trials with Genotype x Environment Interaction},
#>   author = {Sergio Perez-Elizalde and Diego Jarquin and Jose Crossa},
#>   journal = {Journal of Agricultural, Biological, and Environmental Statistics},
#>   year = {2012},
#>   volume = {17},
#>   number = {1},
#>   pages = {15--37},
#>   doi = {10.1007/s13253-011-0063-9},
#> }
```
