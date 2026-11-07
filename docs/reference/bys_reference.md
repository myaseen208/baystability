# Published values used by the tests and vignettes

Reads `inst/extdata/reference.csv`: one row per published number, with
its source, the R call that recomputes it with this package (`call`),
the expression that extracts the number from the result (`extract`, in
terms of `x`), the printed value and the digits to compare. The same
file drives the test of published values and the vignette that
reproduces the published examples. Each distinct `call` is evaluated
once.

## Usage

``` r
bys_reference(source = NULL, evaluate = FALSE, level = 1)
```

## Arguments

- source:

  optional regular expression; keep rows whose `source` matches.

- evaluate:

  logical; evaluate the rows and add columns `computed` and `agrees`
  (default `FALSE`).

- level:

  largest `level` evaluated (1 = fast rows, 2 = slow rows; default 1);
  rows above it get `NA`.

## Value

A `data.table` with columns `id`, `source`, `quantity`, `call`,
`extract`, `value`, `digits`, `level`, `note` (and `computed`,
`agrees`).

## See also

Other baystability utilities:
[`bys_notation()`](https://myaseen208.com/baystability/reference/bys_notation.md),
[`bys_save()`](https://myaseen208.com/baystability/reference/bys_save.md),
[`bys_theme()`](https://myaseen208.com/baystability/reference/bys_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
ref <-
  bys_reference(
      source   = NULL
    , evaluate = FALSE
    , level    = 1
    )
ref
#>                  id                                            source
#>              <char>                                            <char>
#> 1:     perez2012_mu Perez-Elizalde et al. (2012), JABES 17, Table B.1
#> 2: perez2012_alpha7 Perez-Elizalde et al. (2012), JABES 17, Table B.1
#> 3:  perez2012_beta5 Perez-Elizalde et al. (2012), JABES 17, Table B.1
#>                                                                                                          quantity
#>                                                                                                            <char>
#> 1: posterior mean of mu with the first year as prior and equal weights: the average of the two yearly grand means
#> 2:                                          posterior mean of alpha_7: the average of the yearly genotype effects
#> 3:                                        posterior mean of beta_5: the average of the yearly environment effects
#>                                                                                                                                                        call
#>                                                                                                                                                      <char>
#> 1: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#> 2: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#> 3: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#>                                   extract value digits level
#>                                    <char> <int>  <int> <int>
#> 1:                           x$grand_mean  3848      0     1
#> 2: rowMeans(x$ge_means)[7] - x$grand_mean   499      0     1
#> 3: colMeans(x$ge_means)[5] - x$grand_mean  1536      0     1
#>                                                                                                   note
#>                                                                                                 <char>
#> 1:                                                                                                    
#> 2:                                                                                                    
#> 3: the published value is 1537; the posterior mean of a run is within 0.03 SD of it (test-perez2012.R)
# \donttest{
bys_reference(
    source   = "."
  , evaluate = TRUE
  , level    = 1
  )
#>                  id                                            source
#>              <char>                                            <char>
#> 1:     perez2012_mu Perez-Elizalde et al. (2012), JABES 17, Table B.1
#> 2: perez2012_alpha7 Perez-Elizalde et al. (2012), JABES 17, Table B.1
#> 3:  perez2012_beta5 Perez-Elizalde et al. (2012), JABES 17, Table B.1
#>                                                                                                          quantity
#>                                                                                                            <char>
#> 1: posterior mean of mu with the first year as prior and equal weights: the average of the two yearly grand means
#> 2:                                          posterior mean of alpha_7: the average of the yearly genotype effects
#> 3:                                        posterior mean of beta_5: the average of the yearly environment effects
#>                                                                                                                                                        call
#>                                                                                                                                                      <char>
#> 1: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#> 2: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#> 3: bys_ls_means(.data = rbind(cultivo2008[c('site', 'rep', 'entry', 'y')], cultivo2009[c('site', 'rep', 'entry', 'y')]), .y = y, .gen = entry, .env = site)
#>                                   extract value digits level
#>                                    <char> <int>  <int> <int>
#> 1:                           x$grand_mean  3848      0     1
#> 2: rowMeans(x$ge_means)[7] - x$grand_mean   499      0     1
#> 3: colMeans(x$ge_means)[5] - x$grand_mean  1536      0     1
#>                                                                                                   note
#>                                                                                                 <char>
#> 1:                                                                                                    
#> 2:                                                                                                    
#> 3: the published value is 1537; the posterior mean of a run is within 0.03 SD of it (test-perez2012.R)
#>     computed agrees
#>        <num> <lgcl>
#> 1: 3848.2749   TRUE
#> 2:  499.1293   TRUE
#> 3: 1535.5777   TRUE
# }
```
