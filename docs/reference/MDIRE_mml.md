# Fit a Multidimensional DIRE Latent Regression

Fits a multidimensional DIRE model with one subtest per latent
dimension. Main covariates enter the latent regression explicitly, while
residual or raw principal components summarize the remaining
high-dimensional covariate information.

## Usage

``` r
MDIRE_mml(
  X,
  Y,
  parTab,
  main,
  subtest_names = NULL,
  test_name = "comp",
  pca_type = c("cov", "cor"),
  use_residual = TRUE,
  var_threshold = 0.9,
  calcCor = TRUE,
  ...
)
```

## Arguments

- X:

  \`N\` by \`p\` numeric covariate matrix.

- Y:

  \`N\` by \`J_total\` binary response matrix.

- parTab:

  Item parameter data frame with one row per item. For the
  simple-structure format used by \`sim_m1\`, it must contain \`a\` (or
  \`slope\`), \`d\` (or \`b\`/\`difficulty\`), \`c\` (or \`guessing\`),
  and \`dimension\` (or \`subtest\`).

- main:

  Character names or numeric column indices of the main covariates
  retained explicitly in the latent regression.

- subtest_names:

  Optional character vector of \`D\` subtest names.

- test_name:

  Name of the composite test passed to DIRE.

- pca_type:

  Use covariance (\`"cov"\`) or correlation (\`"cor"\`) PCA.

- use_residual:

  If \`TRUE\`, perform PCA on covariate residuals after conditioning on
  the main variables.

- var_threshold:

  Cumulative explained-variance threshold used to retain principal
  components.

- calcCor:

  Logical passed to \[Dire::mml()\].

- ...:

  Additional arguments passed to \[Dire::mml()\].

## Value

An object of class \`MDIREfit\` containing the fitted DIRE object and
all data structures needed by \[MDIRE_drawPVs()\].

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- MDIRE_mml(
  X = sim_m1$X,
  Y = sim_m1$Y,
  parTab = sim_m1$parTab,
  main = sim_m1$main
)
} # }
```
