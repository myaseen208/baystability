# Bayesian AMMI and SREG analysis in a web browser

`bys_app()` starts the application of baystability: load a trial (the
example
[`cultivo2009`](https://myaseen208.com/baystability/reference/cultivo2008.md)
or a CSV file with one row per plot) and, optionally, an earlier trial
of the same genotypes and environments as the prior (the example
[`cultivo2008`](https://myaseen208.com/baystability/reference/cultivo2008.md)
or a CSV file), fit the Bayesian AMMI or SREG (GGE) model of
Perez-Elizalde et al. (2012), check convergence, draw the biplot with
HPD regions and compare models of different numbers of terms, with and
without the earlier trial. Every page shows an interpretation computed
from the results
([`bys_interpret()`](https://myaseen208.com/baystability/reference/bys_interpret.md)),
the output as R prints it, formatted tables, plots and the R code that
repeats the analysis, each with a Copy button.

## Usage

``` r
bys_app(launch.browser = interactive(), ...)

bys_app_object()
```

## Arguments

- launch.browser:

  passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html): `TRUE`
  (default in an interactive session) opens the application in the web
  browser.

- ...:

  further arguments to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), such
  as `port`.

## Value

`bys_app()` is called for its side effect and returns the value of
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) when the
application is closed; `bys_app_object()` returns the application object
(for [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html),
[`shiny::testServer()`](https://rdrr.io/pkg/shiny/man/testServer.html)
or deployment on a server).

## Details

The application needs the suggested packages shiny, bslib and DT. Plots
are sent to the browser as WebP images when magick is installed (smaller
than PNG; PNG otherwise) and can be downloaded as PNG, PDF or WebP. Text
is set in IBM Plex Sans and output in IBM Plex Mono when the browser can
reach Google Fonts, in the system fonts otherwise.

## References

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A General
Bayesian Estimation Method of Linear-Bilinear Models Applied to Plant
Breeding Trials With Genotype x Environment Interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, **17**(1),
15–37.
([doi:10.1007/s13253-011-0063-9](https://doi.org/10.1007/s13253-011-0063-9)
)

## See also

[`bayes_ammi`](https://myaseen208.com/baystability/reference/bayes_ammi.md),
[`bayes_ammi_biplot`](https://myaseen208.com/baystability/reference/bayes_ammi_biplot.md),
[`bys_compare`](https://myaseen208.com/baystability/reference/bys_compare.md),
[`bys_interpret`](https://myaseen208.com/baystability/reference/bys_interpret.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
if (interactive()) {
  bys_app(
      launch.browser = TRUE
    )
}
```
