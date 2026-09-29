# Fit a high-dimensional multidimensional FARLR model

Each latent dimension is first screened separately using the
\`FARLR_Debias\` method in \`Farlr_mml()\`. The union of the
dimension-specific nonzero patterns is stored in a \`P\` by \`D\` index
matrix. Conditional on this fixed sparsity pattern, multidimensional
GVEM iterations jointly update the regression coefficients and latent
residual covariance matrix.

## Usage

``` r
MFARLR_mml(
  X,
  Y,
  parTab,
  farlr_args = list(),
  selection_tol = 0,
  sigma_start = NULL,
  xi_start = NULL,
  max_iter = 1000L,
  threshold = 1e-04,
  ridge = 1e-06,
  progress = TRUE,
  verbose = TRUE
)
```

## Arguments

- X:

  \`N\` by \`p\` covariate matrix.

- Y:

  \`N\` by \`J\` binary response matrix.

- parTab:

  Item parameter data frame with one row per item. For the
  simple-structure format used by \`sim_m1\`, it must contain \`a\` (or
  \`slope\`), \`d\` (or \`b\`/\`difficulty\`), and \`dimension\` (or
  \`subtest\`). A full loading matrix can alternatively be supplied
  through \`a1\`, ..., \`aD\` columns. The current GVEM implementation
  supports multidimensional 2PL items, so \`c\`/\`guessing\`, when
  present, must be zero.

- farlr_args:

  Named list of additional arguments passed to each \`Farlr_mml()\`
  call. Do not include \`X\`, \`Y\`, \`parTab\`, or \`method\`.

- selection_tol:

  Coefficients with absolute value greater than this value are selected
  in \`index\`.

- sigma_start:

  Optional positive-definite \`D\` by \`D\` starting covariance.

- xi_start:

  Optional positive \`N\` by \`J\` starting matrix for \`xi\`.

- max_iter:

  Maximum number of joint GVEM iterations.

- threshold:

  Convergence threshold for the maximum absolute parameter change in
  \`sigma_hat\` and \`beta_hat\`.

- ridge:

  Nonnegative ridge constant used in the selected regressions.

- progress:

  Logical; display a progress bar for the joint GVEM iterations when
  \`TRUE\`.

- verbose:

  Logical; report dimension-specific screening results and the final
  convergence status when \`TRUE\`. Iteration differences are printed
  when \`verbose = TRUE\` and \`progress = FALSE\`.

## Value

A list containing the dimension-specific FARLR fits, selection
\`index\`, joint coefficient matrix, residual covariance, variational
moments, convergence information, the design matrix \`Z\`, and the
original \`X\`, \`Y\`, and \`parTab\` inputs required for
plausible-value generation.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- MFARLR_mml(
  X = sim_m1$X,
  Y = sim_m1$Y,
  parTab = sim_m1$parTab,
  farlr_args = list(main = sim_m1$main)
)
} # }
```
