# Compare Bayesian AMMI and GGE models

`bys_compare()` compares fits of
[`bayes_ammi()`](https://myaseen208.com/baystability/reference/bayes_ammi.md)
that differ in the number of multiplicative terms, the model (AMMI or
GGE), the error variances or the priors, by the deviance information
criterion (DIC) and the widely applicable information criterion (WAIC).
Smaller values predict new data better.

## Usage

``` r
bys_compare(...)

# Default S3 method
bys_compare(...)
```

## Arguments

- ...:

  fits of
  [`bayes_ammi()`](https://myaseen208.com/baystability/reference/bayes_ammi.md)
  to the same data, named or not (names label the rows).

## Value

An object of class `c("bys_compare", "bys_result", "list")` with `table`
(`data.table` sorted by WAIC: `model`, `terms`, `variances`, `DIC`,
`pD`, `WAIC`, `p_waic`, `se` (standard error of WAIC),
`lambda_last_lower` (lower HPD limit of the last singular value),
`lambda_last_mean`, `dWAIC` (difference from the best) and `dWAIC_se`
(its standard error)). Methods:
[`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`bys_interpret()`](https://myaseen208.com/baystability/reference/bys_interpret.md).

## Details

With \\D(\boldsymbol{\theta}) = -2\log p(\mathbf{y} \mid
\boldsymbol{\theta})\\ over the observed plots, \$\$\mathrm{DIC} =
D(\bar{\boldsymbol{\theta}}) + 2p_D, \qquad p_D =
\overline{D(\boldsymbol{\theta})} - D(\bar{\boldsymbol{\theta}})\$\$
(Spiegelhalter et al. 2002), the deviance evaluated at the posterior
means of the cell means and variances; and \$\$\mathrm{WAIC} =
-2\left(\sum_o \log\frac{1}{S}\sum_s p(y_o \mid
\boldsymbol{\theta}^{(s)}) - \sum_o \mathrm{Var}\_s\log p(y_o \mid
\boldsymbol{\theta}^{(s)})\right)\$\$ (Watanabe 2010; Gelman et al.
2014), accumulated over all draws during sampling. All fits must use the
same data.

The standard error of a difference in WAIC comes from the plot-by-plot
differences of the contributions \\w_o\\ to WAIC: \\\mathrm{se} =
\sqrt{n\\\mathrm{Var}\_o(w_o^{(m)} - w_o^{(\mathrm{best})})}\\ (Vehtari
et al. 2017). Because the models are fitted to the same plots, it is
smaller than the standard error of either WAIC. A difference of less
than about two of its standard errors does not separate two models.

## References

Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
(2002). Bayesian measures of model complexity and fit. *Journal of the
Royal Statistical Society B*, **64**(4), 583–639.
[doi:10.1111/1467-9868.00353](https://doi.org/10.1111/1467-9868.00353)

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, **11**, 3571–3594.

Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. *Statistics and Computing*,
**24**(6), 997–1016.
[doi:10.1007/s11222-013-9416-2](https://doi.org/10.1007/s11222-013-9416-2)

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, **27**(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
fits <-
  lapply(
      X   = 1:2
    , FUN = function(K) {
        bayes_ammi(
            .data       = cultivo2009
          , .y          = y
          , .gen        = entry
          , .env        = site
          , .rep        = rep
          , .nIter      = 600
          , .prior_data = NULL
          , .n0         = NULL
          , .model      = "AMMI"
          , .k          = K
          , .burnin     = 100
          , .thin       = 1
          , .chains     = 1
          , .priors     = bys_priors(
                mu_mean     = 0
              , mu_var      = 1e+8
              , tau_mean    = 0
              , tau_var     = 1e+8
              , delta_mean  = 0
              , delta_var   = 1e+8
              , lambda_mean = 0
              , lambda_var  = 1e+8
              , nu          = 1
              , s2          = NULL
              )
          , .level      = 0.95
          , .seed       = 1
          , .align      = TRUE
          , .waic       = TRUE
          )
      }
    )
bys_compare(
    AMMI1 = fits[[1]]
  , AMMI2 = fits[[2]]
  )
#> Comparison of Bayesian models (smaller is better) 
#> 
#>     model terms  prior      DIC      pD     WAIC  p_waic      se
#>    <char> <int> <char>    <num>   <num>    <num>   <num>   <num>
#> 1:  AMMI1     1  vague 14766.74 66.8358 14769.34 64.8932 52.8212
#> 2:  AMMI2     2  vague 14764.29 85.5752 14770.03 84.0827 54.3560
#>    lambda_last_lower lambda_last_mean  dWAIC dWAIC_se
#>                <num>            <num>  <num>    <num>
#> 1:         3202.4996         4232.212 0.0000   0.0000
#> 2:          151.7392         2114.736 0.6854   5.0523
#> 
#> Sources: Spiegelhalter et al. (2002); Watanabe (2010); Vehtari et al. (2017).
#> Full references: bys_cite(x).
```
