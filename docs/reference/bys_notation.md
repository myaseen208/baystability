# Notation of baystability

The symbols used in the help pages and vignettes of baystability, with
their meaning and the argument or result component that holds them. The
notation follows one rule throughout: scalars are upper-case letters in
normal face (\\P\\, \\R\\), indices lower-case italic (\\i, j\\),
vectors lower-case bold (\\\mathbf{y}\\, \\\boldsymbol{\beta}\\) and
matrices upper-case bold (\\\mathbf{X}\\, \\\boldsymbol{\Sigma}\\).

## Usage

``` r
bys_notation(format = c("table", "path", "tex"))
```

## Arguments

- format:

  `"table"` (default: a `data.table` with columns `symbol` (LaTeX),
  `text` (plain text), `kind` (scalar, index, vector, matrix, set,
  function, operator), `meaning` and `argument`), `"path"` (the path of
  `notation.csv`) or `"tex"` (the path of `Notations.tex`, the notation
  appendix used in the articles).

## Value

See `format`.

## See also

Other baystability utilities:
[`bys_reference()`](https://myaseen208.com/baystability/reference/bys_reference.md),
[`bys_save()`](https://myaseen208.com/baystability/reference/bys_save.md),
[`bys_theme()`](https://myaseen208.com/baystability/reference/bys_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bys_notation(
    format = c("table", "path", "tex")[1]
  )
#>                                                                     symbol
#>                                                                     <char>
#>  1:                                                                      G
#>  2:                                                                      E
#>  3:                                                                      B
#>  4:                                                                      K
#>  5:                                                                      n
#>  6:                                                             i, j, l, k
#>  7:                                          \\mathbf{Y} = (\\bar{y}_{ij})
#>  8:                                                                   \\mu
#>  9:                                        \\boldsymbol{\\tau} = (\\tau_i)
#> 10:                                    \\boldsymbol{\\delta} = (\\delta_j)
#> 11:                                                             \\lambda_k
#> 12:                       \\boldsymbol{\\alpha}_k, \\boldsymbol{\\gamma}_k
#> 13:                                                              \\sigma^2
#> 14:                                                     \\mathbf{Y}_0, n_0
#> 15:                                  \\mathbf{U}, \\mathbf{D}, \\mathbf{V}
#> 16:                                                                      S
#> 17:                                          \\mathrm{DIC}, \\mathrm{WAIC}
#> 18:                                                          \\mathbf{K}_n
#> 19: \\sigma^2_\\mu, \\sigma^2_\\tau, \\sigma^2_\\delta, \\sigma^2_\\lambda
#> 20:                                                            \\nu, S^2_0
#> 21:                                                          \\mathrm{SSW}
#> 22:                            \\mathbf{A}^\\top, \\mathbf{1}, \\mathbf{I}
#>                                                                     symbol
#>                                                                     <char>
#>                            text     kind
#>                          <char>   <char>
#>  1:                           G   scalar
#>  2:                           E   scalar
#>  3:                           B   scalar
#>  4:                           K   scalar
#>  5:                           n   scalar
#>  6:                     i j l k    index
#>  7:                           Y   matrix
#>  8:                          mu   scalar
#>  9:                         tau   vector
#> 10:                       delta   vector
#> 11:                    lambda_k   scalar
#> 12:             alpha_k gamma_k   vector
#> 13:                     sigma^2   scalar
#> 14:                       Y0 n0   matrix
#> 15:                       U D V   matrix
#> 16:                           S   scalar
#> 17:                    DIC WAIC   scalar
#> 18:                         K_n   matrix
#> 19: v_mu v_tau v_delta v_lambda   scalar
#> 20:                     nu s0^2   scalar
#> 21:                         SSW   scalar
#> 22:                     A^T 1 I operator
#>                            text     kind
#>                          <char>   <char>
#>                                                                                          meaning
#>                                                                                           <char>
#>  1:                                                                          number of genotypes
#>  2:                                                                       number of environments
#>  3: number of replicates in each environment (largest number of plots in a cell when unbalanced)
#>  4:                                                                     number of bilinear terms
#>  5:                                                                   number of observations BGE
#>  6:                                              genotype, environment, replicate, bilinear term
#>  7:                                                       table of genotype-by-environment means
#>  8:                                                                                 overall mean
#>  9:                                                             genotype main effects (sum zero)
#> 10:                                                          environment main effects (sum zero)
#> 11:                                                            singular value of bilinear term k
#> 12:                               genotype and environment score vectors of term k (unit length)
#> 13:                                                                               error variance
#> 14:                        cell means and weight (replicates) of the earlier trial used as prior
#> 15:                genotype scores, singular values and environment scores of the bilinear terms
#> 16:                                                               number of kept posterior draws
#> 17:                                          deviance and widely applicable information criteria
#> 18:                                     orthonormal contrasts, (n - 1) x n, rows orthogonal to 1
#> 19:                                                                              prior variances
#> 20:                              degrees of freedom and scale of the prior of the error variance
#> 21:                                                                   within-cell sum of squares
#> 22:                                                   transpose, vector of ones, identity matrix
#>                                                                                          meaning
#>                                                                                           <char>
#>        argument
#>          <char>
#>  1:        .gen
#>  2:        .env
#>  3:        .rep
#>  4:          .k
#>  5:            
#>  6:            
#>  7:          .y
#>  8:            
#>  9:            
#> 10:            
#> 11:            
#> 12:            
#> 13:            
#> 14: .prior_data
#> 15:          .k
#> 16:      .nIter
#> 17:       .waic
#> 18:            
#> 19:     .priors
#> 20:     .priors
#> 21:            
#> 22:            
#>        argument
#>          <char>
```
