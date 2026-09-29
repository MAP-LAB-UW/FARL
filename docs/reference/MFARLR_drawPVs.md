# Draw Multidimensional Plausible Values from an MFARLR-GVEM Model

Draws Monte Carlo samples from each respondent's fitted multidimensional
latent-regression prior, weights the samples by the multidimensional 2PL
item-response likelihood, and computes respondent-specific posterior EAP
means and covariance matrices. Plausible values are then drawn from the
multivariate distribution determined by these posterior moments.

## Usage

``` r
MFARLR_drawPVs(
  mmlcomp,
  npv = 10L,
  n_mc = 20000L,
  theta = NULL,
  seed = NULL,
  covariance_floor = 1e-08,
  progress = TRUE,
  verbose = TRUE,
  return_details = TRUE
)
```

## Arguments

- mmlcomp:

  A fitted object returned by \[MFARLR_GVEM()\]. The object retains
  \`X\`, \`Y\`, \`parTab\`, the multidimensional loading matrix, and the
  item intercepts, so these inputs do not need to be supplied again.

- npv:

  Positive integer. Number of plausible values per respondent.

- n_mc:

  Positive integer. Number of prior Monte Carlo samples used for each
  respondent. The original implementation used \`20000\`.

- theta:

  Optional true \`N\` by \`D\` latent-trait matrix used only to report
  the difference between the mean EAP estimate and the mean true score.

- seed:

  Optional random-number seed.

- covariance_floor:

  Small positive eigenvalue used to stabilize posterior covariance
  matrices before drawing plausible values.

- progress:

  Logical; display a respondent-level progress bar.

- verbose:

  Logical; print a completion message and, when \`theta\` is supplied,
  the EAP mean differences.

- return_details:

  Logical. If \`TRUE\`, return posterior summaries along with the
  plausible values. If \`FALSE\`, return only the wide PV data frame.

## Value

When \`return_details = TRUE\`, a list containing:

- data:

  An \`N\` by \`npv\` by \`D\` plausible-value array.

- datPVs:

  A wide data frame with one row per respondent.

- rawPV:

  An \`(N \* npv)\` by \`D\` stacked plausible-value matrix.

- EAP_estimates:

  An \`N\` by \`D\` matrix of posterior EAP means.

- EAP_covariance:

  An \`N\` by \`D\` by \`D\` array of posterior covariance matrices.

- effective_sample_size:

  The importance-sampling effective sample size for every respondent.

- prior_mean:

  The fitted latent-regression prior means.

- sigma:

  The fitted latent residual covariance matrix.

- mean_difference:

  Dimension-specific EAP minus true-score mean differences when
  \`theta\` is supplied.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- MFARLR_GVEM(
  X = sim_m1$X,
  Y = sim_m1$Y,
  parTab = sim_m1$parTab,
  farlr_args = list(main = sim_m1$main)
)

PVs <- MFARLR_drawPVs(
  mmlcomp = fit,
  npv = 10L
)
} # }
```
