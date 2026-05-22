# DIRE Marginal Maximum Likelihood Estimation

Fits a Direct Item Regression Effects (DIRE) model using a marginal
maximum likelihood (MML) framework. The DIRE model extends standard item
response theory by allowing item parameters to depend directly on
person-level covariates, enabling the assessment of covariate-induced
differential item functioning (DIF). This function provides a unified
estimation interface for dichotomous and polytomous item responses,
supports complex survey designs, and allows flexible numerical
integration and optimization options.

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
  verbose = 0
)
```

## Arguments

- formula:

  Formula. A model formula specifying person-level covariates entering
  the DIRE model. The left-hand side is ignored; the right-hand side
  defines covariates with potential direct effects on item parameters.

- stuItems:

  Data frame. Long-format item response data containing item identifiers
  and response values for each individual.

- stuDat:

  Data frame. Person-level data frame containing covariates referenced
  in `formula`.

- idVar:

  Character. Name of the variable in `stuItems` identifying individuals.

- dichotParamTab:

  Data frame or `NULL`. Item parameter table for dichotomous items. If
  `NULL`, dichotomous items are not modeled.

- polyParamTab:

  Data frame or `NULL`. Item parameter table for polytomous items. If
  `NULL`, polytomous items are not modeled.

- testScale:

  Character vector or `NULL`. Optional specification of test or scale
  membership for items.

- Q:

  Integer. Number of quadrature nodes used for numerical integration.
  Defaults to `30`.

- minNode:

  Numeric. Lower bound of the quadrature nodes. Defaults to `-4`.

- maxNode:

  Numeric. Upper bound of the quadrature nodes. Defaults to `4`.

- polyModel:

  Character. Polytomous item response model to be used. Supported
  options include `"GPCM"` and `"GRM"`.

- weightVar:

  Character or `NULL`. Optional sampling weight variable in `stuDat`
  used for weighted likelihood estimation.

- multiCore:

  Logical. If `TRUE`, enables parallel computation for likelihood
  evaluation where supported. Defaults to `FALSE`.

- bobyqaControl:

  List or `NULL`. Optional control parameters passed to the `bobyqa`
  optimizer for numerical optimization.

- composite:

  Logical. If `TRUE`, uses a composite likelihood approximation for
  estimation. Defaults to `TRUE`.

- strataVar:

  Character or `NULL`. Optional stratification variable for complex
  survey designs.

- PSUVar:

  Character or `NULL`. Optional primary sampling unit (PSU) variable for
  complex survey designs.

- fast:

  Logical. If `TRUE`, uses computational shortcuts to accelerate
  estimation. Defaults to `TRUE`.

- calcCor:

  Logical. If `TRUE`, computes correlation matrices for estimated item
  effects. Defaults to `TRUE`.

- verbose:

  Integer. Verbosity level controlling diagnostic output. `0` suppresses
  output; larger values produce more detailed messages.

## Value

A list containing estimation results from the DIRE model, typically
including:

- `item.par`:

  Estimated baseline item parameters.

- `dire.coef`:

  Estimated direct item regression effect coefficients associated with
  person-level covariates.

- `vcov`:

  Estimated variance–covariance matrix of parameter estimates.

- `LogLik`:

  Maximized marginal (or composite) log-likelihood value.

- `convergence`:

  Indicator of convergence status of the numerical optimization.

## Details

The DIRE model allows person-level covariates to enter item parameter
models directly, providing a flexible framework for assessing
covariate-related DIF without requiring anchor items. Latent variables
are integrated out using Gaussian quadrature with `Q` nodes over the
interval \[`minNode`, `maxNode`\]. For large-scale assessments or
complex survey data, composite likelihood and survey design adjustments
can be employed to improve computational feasibility.

## See also

[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html),
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`optim`](https://rdrr.io/r/stats/optim.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- Dire_mml(
  formula = ~ gender + ses,
  stuItems = stuItems,
  stuDat = stuDat,
  idVar = "student_id",
  dichotParamTab = dichotTab,
  polyParamTab = polyTab,
  polyModel = "GPCM",
  Q = 30,
  verbose = 1
)
} # }
```
