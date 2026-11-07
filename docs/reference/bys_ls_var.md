# Genotype-by-environment variances

`bys_ls_var()` computes the variance of the replicates of every
genotype-by-environment cell.

## Usage

``` r
bys_ls_var(.data, .y, .gen, .env)

# Default S3 method
bys_ls_var(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bys_ge_var", "bys_result", "list")`) with
`ge_variances`, the \\G \times E\\ matrix of variances.

## Details

\$\$s^2\_{ij} = \frac{1}{n\_{ij} - 1}\sum_l (y\_{ijl} -
\bar{y}\_{ij})^2,\$\$ `NA` for a cell with one replicate. The
within-cell sum of squares of the Bayesian AMMI model is
\\\sum\_{ij}(n\_{ij} - 1)s^2\_{ij}\\.

**Change from the earlier versions.**
[`bayes_ammi()`](https://myaseen208.com/baystability/reference/bayes_ammi.md)
used the sum of these variances in place of the within-cell sum of
squares.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`bys_ls_means`](https://myaseen208.com/baystability/reference/bys_ls_means.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(cultivo2009)
bys_ls_var(
    .data = cultivo2009
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#> <bys_ge_var> Genotype-by-environment variances
#> 
#> $ge_variances
#>          1      2        3      4       5       6       7       8      9
#> 2    39470 274034  2503325 453518   41830 2332382  959741  757317 784856
#> 10  289824 722742 10375176 310521  294905  409686  173121  304449 101951
#> 12  625859 909008  2112686 349242   83504 1225126 1288101  192972   2021
#> 13  608856 214628    36129 175845  653386 2858682  490624   94485 175958
#> 14  444901 347144  3910158 453039  463206  281084  198444  107160  42028
#> 18 2310584 551958   499273 688415    4765   54840  646333  331180  33711
#> 19  192938 884172  3785782  21475   11712   78995 2290000 1212527  51885
#> 20 1026582 506201  1059297 128378  104915  826646  780541  634708 131675
#> 23  885776   5008  3354400 475354 1901210  306905  536543  403222  92217
#> 38  411810 410029  3685413 393989   61360 1805145  138332   39535  11653
#> 40  349649  22194  3568365 185385 1359236 2027615 1151782  288038  13509
#> 42   11607 349160   442895 850647   67922  370990  824088  856355 556514
#>         10     11      12      13     14       15      16      17     18
#> 2   867453 176309   94701  949208  33100  40742.7  506377  151821 703773
#> 10  476984 651362   96167 1806196  21812 320322.9  664010  445669 224306
#> 12 2982145  34282  148879 1327817 431886  73901.2  769889  367309  55392
#> 13 1330566   4897   12268  154456  87146 536509.0   19224 1488828 495851
#> 14   10864 181206  144769   73226  28597  82999.5 8666355   19590 964905
#> 18 1121722 239976  105258   72402 165340   1497.2  792097   93052  75726
#> 19 2644768  44077  502340  521185 152962    513.6   27811  573003 230956
#> 20 1391869 191001 1009629   39909 124785 120173.4 1244444 4250995 445702
#> 23 2153537 181206   55232  194625   9457   1195.8   90796  455464  34918
#> 38   68092 573003  270714  499200 387217  19610.4  773173  445669 616459
#> 40  855526 132231   95044  145059  64132  89257.8 6052591  607285 140243
#> 42 1537852 210591   75565  240540 179038   7091.3  118774 1621059 132322
#>         19      20      21      22      23      24      25
#> 2   687988 1601469  102847  634458 1140635    9816  993103
#> 10  824939  573003  631772 1965501 1974687  357254  663457
#> 12  589368  357515 1337006  253693  200954  257821  384286
#> 13  828464  705234  826446   15575  296991  526396 1495489
#> 14  396977 2659320  411387   97857  618060 1023504  419849
#> 18  415958 1513315 2256504   41765  595486  433408  468408
#> 19 1920226  445669 1404346  382047  336106  573418 1547984
#> 20  467069  396694  181206  576807  137216  109872 1637276
#> 23   15980 1973676  831344  311796   27944  748192 1103218
#> 38 2418624 3310683  210591   78790   29552 4240196  475614
#> 40  302810  357515  559535   71769  313142   79716  360778
#> 42   43365  455464  179982   10722   99371   24938    1352
```
