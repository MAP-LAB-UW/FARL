# Draw Plausible Values from a Fitted DIRE Model

Generates plausible values (PVs) of latent traits from a fitted Direct
Item Regression Effects (DIRE) model. Plausible values are random draws
from the posterior distribution of the latent variable given observed
responses and estimated model parameters, and are commonly used for
secondary analyses to properly account for measurement uncertainty.

## Usage

``` r
Dire_drawPVs(x, npv, pvVariableNameSuffix = "_dire", ...)
```

## Arguments

- x:

  Object. A fitted DIRE model object returned by
  [`Dire_mml`](https://yijunchenguw.github.io/FARL/reference/Dire_mml.md).

- npv:

  Integer. Number of plausible values to draw for each individual.

- pvVariableNameSuffix:

  Character. Suffix appended to the names of generated plausible value
  variables. Defaults to `"_dire"`.

- ...:

  Additional arguments passed to internal sampling routines.

## Value

A data frame containing the generated plausible values. Each row
corresponds to an individual, and each column corresponds to a plausible
value draw. Column names follow the pattern `PV1<pvVariableNameSuffix>`,
`PV2<pvVariableNameSuffix>`, …, `PVnpv<pvVariableNameSuffix>`.

## Details

Plausible values are drawn by sampling from the posterior distribution
of the latent trait implied by the fitted DIRE model, conditional on
observed item responses and estimated item and regression parameters.
The resulting plausible values can be used in downstream analyses (e.g.,
regression, group comparisons, or secondary modeling) by combining
results across draws using Rubin’s rules or other multiple imputation
techniques.

This function is intended for post-estimation use and does not refit the
DIRE model.

## See also

[`Dire_mml`](https://yijunchenguw.github.io/FARL/reference/Dire_mml.md),
[`fscores`](https://philchalmers.github.io/mirt/reference/fscores.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- Dire_mml(
  formula = ~ gender + ses,
  stuItems = stuItems,
  stuDat = stuDat,
  idVar = "student_id",
  dichotParamTab = dichotTab
)

pv <- Dire_drawPVs(
  x = fit,
  npv = 5
)
} # }
```
