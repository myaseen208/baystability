# Maize trials of 2008 and 2009

Grain yield of 12 maize hybrids in 25 environments in two consecutive
years, the multi-environment trial of Perez-Elizalde et al. (2012),
Section 3.2: `cultivo2008` (the first year, used as the prior) and
`cultivo2009` (the second year). Randomized complete blocks with 3
replicates in every environment.

## Format

Data frames with 900 rows (12 genotypes x 25 environments x 3
replicates) and 4 columns:

- `site`:

  environment, 1 to 25 (the same environments in both years).

- `rep`:

  replicate, 1 to 3.

- `entry`:

  genotype, 1 to 12.

- `y`:

  grain yield (kg/ha).

## Source

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general
Bayesian estimation method of linear-bilinear models applied to plant
breeding trials with genotype x environment interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, **17**(1),
15–37.
[doi:10.1007/s13253-011-0063-9](https://doi.org/10.1007/s13253-011-0063-9)

## Details

The paper describes two replicates; the data have three. The posterior
means of its Table B.1 are, to within a posterior standard deviation of
a hundredth, the averages of the least-squares estimates of the two
years, as the method gives with equal weights for both years (see the
vignette on the published analysis). The help page of version 0.2.0
described another data set.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2008)
data(cultivo2009)
str(cultivo2009)
#> tibble [900 × 4] (S3: tbl_df/tbl/data.frame)
#>  $ site : int [1:900] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ entry: int [1:900] 2 2 2 10 10 10 12 12 12 13 ...
#>  $ rep  : int [1:900] 1 2 3 1 2 3 1 2 3 1 ...
#>  $ y    : num [1:900] 5079 4899 5296 5353 4353 ...
#>  - attr(*, "spec")=List of 2
#>   ..$ cols   :List of 4
#>   .. ..$ site : list()
#>   .. .. ..- attr(*, "class")= chr [1:2] "collector_integer" "collector"
#>   .. ..$ entry: list()
#>   .. .. ..- attr(*, "class")= chr [1:2] "collector_integer" "collector"
#>   .. ..$ rep  : list()
#>   .. .. ..- attr(*, "class")= chr [1:2] "collector_integer" "collector"
#>   .. ..$ y    : list()
#>   .. .. ..- attr(*, "class")= chr [1:2] "collector_double" "collector"
#>   ..$ default: list()
#>   .. ..- attr(*, "class")= chr [1:2] "collector_guess" "collector"
#>   ..- attr(*, "class")= chr "col_spec"
bys_ls_means(
    .data = cultivo2009
  , .y    = y
  , .gen  = entry
  , .env  = site
  )$grand_mean
#> [1] 3627.3
```
