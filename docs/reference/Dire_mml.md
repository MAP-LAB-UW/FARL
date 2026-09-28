# DIRE Marginal Maximum Likelihood Estimation with Optional PCA

Fits a DIRE model with \`Dire::mml()\`. When \`X\` and \`main_vars\` are
supplied, the function first extracts PCs with
\`PCA_extraction_merge()\`, merges the PC scores into \`stuDat\`,
appends the PC terms to \`formula\`, and then generates the fitted
\`mmlcomp\` object.

## Usage

``` r
Dire_mml(
  formula,
  stuItems,
  stuDat,
  idVar,
  dichotParamTab = NULL,
  polyParamTab = NULL,
  testScale = NULL,
  Q = 30,
  minNode = -4,
  maxNode = 4,
  polyModel = c("GPCM", "GRM"),
  weightVar = NULL,
  multiCore = FALSE,
  bobyqaControl = NULL,
  composite = TRUE,
  strataVar = NULL,
  PSUVar = NULL,
  fast = TRUE,
  calcCor = TRUE,
  verbose = 0,
  X = NULL,
  main_vars = NULL,
  pca_type = c("cov", "cor"),
  use_residual = TRUE,
  var_threshold = 0.9
)
```

## Arguments

- formula:

  Formula passed to \`Dire::mml()\`. It should already contain the main
  covariates that are to remain explicit in the latent regression.

- stuItems:

  Long-format item response data.

- stuDat:

  Person-level data. Its row order must match the row order of \`X\`
  when PCA is requested.

- idVar:

  Character name of the individual identifier.

- dichotParamTab:

  Dichotomous item parameter table or \`NULL\`.

- polyParamTab:

  Polytomous item parameter table or \`NULL\`.

- testScale:

  Optional test-scale table.

- Q:

  Number of quadrature nodes.

- minNode, maxNode:

  Lower and upper quadrature bounds.

- polyModel:

  Polytomous model, \`"GPCM"\` or \`"GRM"\`.

- weightVar:

  Optional sampling-weight variable.

- multiCore:

  Whether to use supported parallel calculations.

- bobyqaControl:

  Optional optimizer controls.

- composite:

  Whether to use composite likelihood.

- strataVar, PSUVar:

  Optional complex-survey variables.

- fast:

  Whether to use DIRE computational shortcuts.

- calcCor:

  Whether to calculate the coefficient correlation matrix.

- verbose:

  DIRE verbosity level.

- X:

  Optional numeric covariate matrix/data frame used for PCA. Leave
  \`NULL\` to fit the original non-PCA DIRE model.

- main_vars:

  Character vector naming the main variables in \`X\`.

- pca_type:

  PCA based on covariance (\`"cov"\`) or correlation (\`"cor"\`).

- use_residual:

  Whether to perform residual PCA.

- var_threshold:

  Cumulative explained-variance threshold for selecting PCs.

## Value

The fitted DIRE \`mmlcomp\` object. When PCA is used, the returned
object also contains \`pca_object\` and \`formula_before_pca\`
components; its standard \`formula\` and \`stuDat\` components contain
the augmented versions.
