# baystability: Bayesian Stability Analysis of Genotype by Environment Interaction (GEI)

baystability implements the general Bayesian estimation method of
linear-bilinear models of Perez-Elizalde, Jarquin and Crossa (2012) for
multi-environment trials: the AMMI model and the sites regression (SREG,
GGE) model, with the data of an earlier trial as the prior or with vague
priors. Every parameter has a posterior distribution, so the biplot
carries inference: highest posterior density (HPD) regions show which
genotypes and environments interact significantly. For Bayesian
stability analysis in general (Finlay-Wilkinson, Eberhart-Russell and
Shukla models, stability indices, probabilities of superior performance)
see the package bayesstability; for the Bayesian AMMI model of Crossa et
al. (2011) see bayesammi.

## Details

**Functions**, in the order of an analysis:

- [`bys_ls_means`](https://myaseen208.com/baystability/reference/bys_ls_means.md),
  [`bys_ls_var`](https://myaseen208.com/baystability/reference/bys_ls_var.md),
  [`bys_ls_gen_effects`](https://myaseen208.com/baystability/reference/bys_ls_gen_effects.md),
  [`bys_ls_env_effects`](https://myaseen208.com/baystability/reference/bys_ls_env_effects.md),
  [`bys_ls_ge_effects`](https://myaseen208.com/baystability/reference/bys_ls_ge_effects.md):
  cell means and variances, main effects, interaction effects and their
  singular value decomposition.

- [`bys_ls_ammi`](https://myaseen208.com/baystability/reference/bys_ls_ammi.md),
  [`bys_ls_model`](https://myaseen208.com/baystability/reference/bys_ls_model.md):
  least-squares AMMI and the mixed model of the error variance.

- [`bys_priors`](https://myaseen208.com/baystability/reference/bys_priors.md),
  [`bayes_ammi`](https://myaseen208.com/baystability/reference/bayes_ammi.md):
  vague priors and the Gibbs sampler, with an earlier trial as prior
  (`.prior_data`).

- [`bayes_ammi_biplot`](https://myaseen208.com/baystability/reference/bayes_ammi_biplot.md):
  the biplot with HPD regions.

- [`bys_compare`](https://myaseen208.com/baystability/reference/bys_compare.md):
  DIC and WAIC of several fits (for example numbers of terms).

- [`bys_interpret`](https://myaseen208.com/baystability/reference/bys_interpret.md),
  [`bys_cite`](https://myaseen208.com/baystability/reference/bys_cite.md):
  a plain-language interpretation and the references of a result.

- [`bys_contrasts`](https://myaseen208.com/baystability/reference/bys_contrasts.md),
  [`bys_orthonormal`](https://myaseen208.com/baystability/reference/bys_orthonormal.md):
  orthonormal contrasts and Gram-Schmidt orthonormalisation.

Every analysis returns the named list documented on its help page, with
a class and methods (`print`, `summary`, `plot`, `as.data.table`); see
[`bys-methods`](https://myaseen208.com/baystability/reference/bys-methods.md).
Utilities:
[`bys_theme`](https://myaseen208.com/baystability/reference/bys_theme.md),
[`bys_save`](https://myaseen208.com/baystability/reference/bys_save.md),
[`bys_notation`](https://myaseen208.com/baystability/reference/bys_notation.md),
[`bys_reference`](https://myaseen208.com/baystability/reference/bys_reference.md).

**Data.**
[`cultivo2008`](https://myaseen208.com/baystability/reference/cultivo2008.md),
[`cultivo2009`](https://myaseen208.com/baystability/reference/cultivo2008.md):
grain yield of 12 maize hybrids in 25 environments in two years
(Perez-Elizalde et al. 2012).

**Vignettes.** Start with
[`vignette("baystability_01_getting_started", package = "baystability")`](https://myaseen208.com/baystability/articles/baystability_01_getting_started.md);
`browseVignettes("baystability")` lists all of them, each as HTML and as
PDF.

## Notation

Scalars are upper-case letters in normal face, indices lower-case
italic, vectors lower-case bold and matrices upper-case bold. The full
table is
[`bys_notation()`](https://myaseen208.com/baystability/reference/bys_notation.md).

|  |  |
|----|----|
| \\G\\, \\E\\, \\B\\, \\K\\ | numbers of genotypes, environments, replicates, bilinear terms |
| \\\mathbf{Y}\\, \\\mathbf{Y}\_0\\ | tables of cell means of the trial and of the earlier trial |
| \\n_0\\ | weight (replicates) of the earlier trial |
| \\\mu\\, \\\boldsymbol{\tau}\\, \\\boldsymbol{\delta}\\ | overall mean, genotype and environment effects |
| \\\lambda_k\\, \\\boldsymbol{\alpha}\_k\\, \\\boldsymbol{\gamma}\_k\\ | singular value and genotype and environment score vectors of term \\k\\ |
| \\\sigma^2\\ | error variance |

## References

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general
Bayesian estimation method of linear-bilinear models applied to plant
breeding trials with genotype x environment interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, **17**(1),
15–37.
[doi:10.1007/s13253-011-0063-9](https://doi.org/10.1007/s13253-011-0063-9)

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

## See also

[`bys-methods`](https://myaseen208.com/baystability/reference/bys-methods.md),
[`bys_notation`](https://myaseen208.com/baystability/reference/bys_notation.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>), ORCID 0000-0002-5923-1714.
