# FARL: An Applied Guide for Unidimensional Assessments

## Overview

The `FARL` package fits factor-augmented regularized latent regression
models for unidimensional assessments. It is designed for applications
with many correlated background variables. The common variation in those
variables is summarized by factor scores, while regularization
identifies additional covariates that are useful for predicting the
latent trait.

The applied workflow has two main steps:

1.  fit the population model with
    [`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md);
    and
2.  generate plausible values with
    [`Farlr_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md).

Two estimation methods are available.

| Method | When to use it |
|:---|:---|
| `FARLR_EMM` | Use when the goal is a sparse model estimated through the iterative regularized procedure. |
| `FARLR_Debias` | Use when the goal is variable selection followed by a correction of regularization shrinkage. This is a useful default for reporting coefficient estimates. |

The examples below focus on data preparation, function calls, output
interpretation, and plausible-value analysis rather than algorithmic
details.

## Load FARL

``` r
library(FARL)
```

## Prepare the data

[`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md)
uses three required data objects.

| Object | Required organization | Purpose |
|:---|:---|:---|
| `X` | Respondents by covariates | Background variables used in the population model. All columns must be numeric, finite, and have nonzero variance. |
| `Y` | Respondents by items | Item responses. Rows must correspond exactly to the rows of `X`. Missing responses may be coded `NA`. |
| `parTab` | One row per item | Fixed item parameters. Its row order must match the column order of `Y`. |

The package supports binary 1PL, 2PL, and 3PL items, as well as
three-category GPCM items. Binary responses must be coded 0/1, and GPCM
responses must be coded 0/1/2. The package determines item type from the
supplied item parameters.

### Example data

`sim_a1` contains 300 respondents, 100 background variables, and six
binary 2PL items. The true latent trait is retained only for simulation
diagnostics.

``` r
data(sim_a1)

dim(sim_a1$X)
# [1] 300 100

dim(sim_a1$Y)
# [1] 300   6

head(sim_a1$Y)
head(sim_a1$X[, 1:6])
sim_a1$parTab
```

Before fitting the model, verify row alignment, item alignment, response
codes, and missingness.

``` r
stopifnot(nrow(sim_a1$X) == nrow(sim_a1$Y))
stopifnot(nrow(sim_a1$parTab) == ncol(sim_a1$Y))

mean(is.na(sim_a1$Y))
table(sim_a1$Y, useNA = "ifany")

apply(sim_a1$X, 2, sd) |>
  summary()
```

### Choose primary covariates

The `main` argument identifies covariates that must remain in the
population model. These variables are not removed by penalized
selection. Supply their column indices in `X`.

``` r
main <- c(1, 2, 15, 29, 45)

main_names <- if (is.null(colnames(sim_a1$X))) {
  paste0("X", main)
} else {
  colnames(sim_a1$X)[main]
}

main_names
```

Use `main = integer(0)` when every original covariate may be selected or
excluded by regularization.

## Fit FARL

### Main fitting arguments

| Argument | Applied meaning |
|:---|:---|
| `X`, `Y`, `parTab` | Data objects described above. |
| `method` | Selects `"FARLR_EMM"` or `"FARLR_Debias"`. |
| `main` | Column indices of primary covariates that remain in the model. |
| `lambda` | Candidate regularization values. A wider grid is useful when the selected value falls at an endpoint. |
| `n_sam` | Monte Carlo sample size used during fitting. Larger values reduce simulation noise but increase runtime. |
| `n_sam_final` | Monte Carlo sample size used in the final evaluation or debiased update. |
| `K_hat` | Optional number of common factors. When omitted, FARL estimates it using parallel analysis. |
| `delta.criteria` | Convergence tolerance for iterative updates. |
| `iter.max` | Maximum number of iterations for each candidate regularization value. |
| `window.size` | Number of recent updates used for stabilization. |
| `proposal_inflation` | Controls the spread of the latent-trait proposal distribution used during Monte Carlo integration. |
| `seed` | Makes factor selection and Monte Carlo calculations reproducible. |
| `progress`, `verbose` | Control progress bars and printed fitting information. |

The defaults provide a convenient starting point. For a final analysis,
inspect convergence and check whether results remain stable with larger
Monte Carlo sample sizes and a broader `lambda` grid.

### Fit FARLR-Debias

``` r
fit_debias <- with(
  sim_a1,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_Debias",
    main = main,
    lambda = seq(0.05, 0.50, length.out = 10),
    n_sam = 25L,
    n_sam_final = 60L,
    seed = 2026,
    progress = TRUE
  )
)
```

### Fit FARLR-EMM

The input interface is the same; only `method` changes.

``` r
fit_emm <- with(
  sim_a1,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_EMM",
    main = main,
    lambda = seq(0.05, 0.50, length.out = 10),
    n_sam = 25L,
    n_sam_final = 60L,
    seed = 2026,
    progress = TRUE
  )
)
```

## Examine a fitted FARL model

The fitted object contains the coefficient estimates together with the
factor-adjusted design, item information, scaling values, and fitting
diagnostics needed for later plausible-value generation.

| Component | Interpretation |
|:---|:---|
| `coefficients` | Named regression coefficients. Names beginning with `FS` refer to common-factor scores; names beginning with `X` refer to original covariates. |
| `sigma` | Residual standard deviation of the latent trait after conditioning on the fitted predictors. |
| `lambda` | Selected regularization value. |
| `LogLik` | Objective or likelihood summary returned for the selected model. |
| `Convergence` | Method-specific convergence indicator. |
| `factor_scores` | Estimated common-factor scores for each respondent. |
| `K_hat` | Number of common factors used in the fitted model. |
| `hatU` | Factor-adjusted covariates used by `FARLR_Debias`. |
| `Fan` | Final factor-plus-adjusted-covariate design used by `FARLR_Debias`. |
| `X`, `X_center`, `X_scale` | Standardized covariates and the values used to standardize them. |
| `item_params`, `item_type` | Item parameters and detected response-model type. |
| `stuDat`, `stuItems` | Subject-level design information and long-format item responses retained for compatibility and PV generation. |

### Check method and convergence

``` r
fit_debias$method
fit_debias$K_hat
fit_debias$lambda
fit_debias$LogLik
fit_debias$Convergence

names(fit_debias)
```

If the selected `lambda` is the smallest or largest candidate, extend
the grid and refit. If the model does not converge, increase `iter.max`,
review the input variables, and assess stability across seeds and Monte
Carlo settings.

### Summarize coefficients

Printing every coefficient is difficult when `X` is wide. A compact
table of the largest coefficients is usually more informative.

``` r
coefficient_table <- data.frame(
  predictor = names(fit_debias$coefficients),
  coefficient = as.numeric(fit_debias$coefficients),
  row.names = NULL
)

coefficient_table$source <- ifelse(
  grepl("^FS", coefficient_table$predictor),
  "Common factor",
  "Original covariate"
)

coefficient_table <- coefficient_table[
  order(abs(coefficient_table$coefficient), decreasing = TRUE),
]

head(coefficient_table, 15L)
```

A positive coefficient indicates that higher values of the predictor are
associated with a higher population mean for the latent trait, holding
the other included predictors constant. Coefficients describe
associations and should not be interpreted causally.

To inspect variables retained by the fitted model, apply a small
numerical tolerance rather than comparing floating-point values directly
with zero.

``` r
selection_tolerance <- 1e-8

selected_coefficients <- subset(
  coefficient_table,
  abs(coefficient) > selection_tolerance
)

selected_coefficients
```

The common-factor scores summarize shared variation across `X`.
Covariate coefficients represent information that remains useful after
the factor adjustment. Primary variables supplied through `main` should
be interpreted alongside the selected covariates rather than as a
separate model.

### Interpret residual variation

``` r
fit_debias$sigma
fit_emm$sigma
```

`sigma` is a residual standard deviation, not a variance. A smaller
value indicates that the fitted population model leaves less
latent-trait variation unexplained. Comparisons are most useful when the
same data, item parameters, and latent scale are used.

### Inspect the factor-adjusted design

``` r
dim(fit_debias$factor_scores)
head(fit_debias$factor_scores)

dim(fit_debias$hatU)
dim(fit_debias$Fan)
head(colnames(fit_debias$Fan))

fit_debias$X_center[1:6]
fit_debias$X_scale[1:6]
```

`FARLR_Debias` uses `Fan`, which combines common-factor scores and
factor-adjusted covariates. These stored objects are especially
important when applying the model to new data or diagnosing unexpected
coefficients.

### Compare the two estimation methods

``` r
method_comparison <- data.frame(
  method = c(fit_emm$method, fit_debias$method),
  lambda = c(fit_emm$lambda, fit_debias$lambda),
  sigma = c(fit_emm$sigma, fit_debias$sigma),
  nonzero_coefficients = c(
    sum(abs(fit_emm$coefficients) > selection_tolerance),
    sum(abs(fit_debias$coefficients) > selection_tolerance)
  ),
  row.names = NULL
)

method_comparison
```

Different coefficient counts are expected because the two methods
perform estimation and shrinkage correction differently. Method choice
should be made before inspecting the final results and should be checked
through sensitivity analyses when substantive conclusions depend on a
small set of variables.

## Generate plausible values

[`Farlr_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md)
combines the fitted population model with each respondent’s observed
item responses. It returns plausible values together with posterior
summaries that are useful for checking the generated draws. Plausible
values are intended for population and subgroup analyses, not for
reporting an individual’s single score.

### Plausible-value arguments

| Argument | Applied meaning |
|:---|:---|
| `mmlcomp` | Object returned by [`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md). |
| `npv` | Number of plausible values generated for each respondent. |
| `n_quad` | Number of points used to approximate each respondent’s posterior distribution. |
| `theta_range` | Lower and upper limits of the initial latent-score grid. |
| `adaptive_range` | When `TRUE`, allows the grid to expand for respondents whose fitted prior lies near or beyond the initial limits. |
| `range_sd` | Controls how far the adaptive grid extends around the fitted prior mean. |
| `design_matrix` | Optional replacement design matrix. Normally omitted so that FARL reconstructs the fitted design automatically. |
| `theta` | Optional true trait values used only for simulation diagnostics. Omit this argument for real data. |
| `seed` | Makes plausible-value generation reproducible. |
| `return_details` | When `TRUE`, returns the PV matrix and posterior summaries in a list. |

### Draw plausible values

``` r
farl_pvs <- Farlr_drawPVs(
  mmlcomp = fit_debias,
  npv = 10L,
  n_quad = 30L,
  theta_range = c(-4, 4),
  adaptive_range = TRUE,
  range_sd = 4,
  seed = 2026,
  progress = TRUE,
  return_details = TRUE
)
```

For simulation studies, `theta = sim_a1$theta` may be added to report
the difference between the average posterior estimate and the average
true trait. It should not be supplied in an empirical analysis.

### Examine plausible-value output

| Component | Interpretation |
|:---|:---|
| `data` | Numeric matrix containing only the plausible values. |
| `datPVs` | Data frame containing respondent IDs and columns `_farl1`, `_farl2`, and so on. |
| `EAP_estimates` | Respondent-level posterior means. |
| `EAP_variance` | Respondent-level posterior variances. |
| `posterior_mode` | Most likely grid value for each respondent. |
| `prior_mean` | Latent-regression mean predicted by the fitted population model. |
| `sigma` | Residual standard deviation copied from the fitted model. |
| `mean_difference` | Simulation-only mean recovery diagnostic; `NA` when `theta` is omitted. |

``` r
names(farl_pvs)
dim(farl_pvs$data)

head(farl_pvs$datPVs)
summary(farl_pvs$EAP_estimates)
summary(farl_pvs$EAP_variance)

cor(farl_pvs$EAP_estimates, farl_pvs$prior_mean)
```

Large posterior variances identify respondents whose response patterns
provide less information about the latent trait. Differences between
`prior_mean` and `EAP_estimates` show how the observed item responses
update predictions from the background variables.

## Use plausible values in secondary analyses

Run the intended analysis separately for every plausible-value column
and then combine the resulting estimates and uncertainty. Do not average
the plausible values first and analyze that single average as if it were
an observed score.

The following helpers summarize the overall spread across a set of
plausible values. `variance.sample()` averages the respondent-level
variance across PV datasets, while `variance.cal()` also accounts for
variation among the PV dataset means.

``` r
variance.sample <- function(theta_matrix) {
  theta_matrix <- as.matrix(theta_matrix)
  mean(apply(theta_matrix, 2L, stats::var))
}

variance.cal <- function(theta_matrix) {
  theta_matrix <- as.matrix(theta_matrix)
  number_of_pvs <- ncol(theta_matrix)
  within_variance <- variance.sample(theta_matrix)
  between_variance <- stats::var(colMeans(theta_matrix))

  within_variance +
    (1 + 1 / number_of_pvs) * between_variance
}

pv_matrix <- farl_pvs$data

variance.sample(pv_matrix)
variance.cal(pv_matrix)
```

For a substantive regression, repeat the model for each plausible value
and store the coefficient and its sampling variance.

``` r
analysis_data <- data.frame(
  outcome = sim_a1$X[, 1],
  group = sim_a1$X[, 2]
)

pv_results <- lapply(
  seq_len(ncol(farl_pvs$data)),
  function(pv_index) {
    model_data <- transform(
      analysis_data,
      plausible_value = farl_pvs$data[, pv_index]
    )

    model <- lm(plausible_value ~ group, data = model_data)

    data.frame(
      pv = pv_index,
      estimate = coef(model)["group"],
      sampling_variance = vcov(model)["group", "group"],
      row.names = NULL
    )
  }
)

pv_results <- do.call(rbind, pv_results)
pv_results
```

These stored estimates can then be combined using the
multiple-imputation pooling procedure appropriate for the target
analysis.

## Mixed 3PL and GPCM responses

The same FARL workflow applies when the assessment contains both binary
3PL items and three-category GPCM items. `sim_a2` demonstrates this
case.

``` r
data(sim_a2)

dim(sim_a2$X)
dim(sim_a2$Y)
table(sim_a2$itemtype)

sim_a2$parTab[, c(
  "ItemID", "itemtype", "a", "b", "c", "b1", "b2"
)]
```

Fit the model and generate PVs using the same public functions.

``` r
fit_mixed <- with(
  sim_a2,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_Debias",
    main = main,
    seed = 2026
  )
)

table(fit_mixed$item_type)

pvs_mixed <- Farlr_drawPVs(
  mmlcomp = fit_mixed,
  npv = 10L,
  adaptive_range = TRUE,
  seed = 2026
)

head(pvs_mixed$datPVs)
```

For mixed-format assessments, check that every `parTab` row contains the
parameters required by its item type. Missing responses remain coded as
`NA` and are ignored for the corresponding item during posterior
calculation.

## Practical checks

Before reporting a FARL analysis, verify the following:

- rows of `X` and `Y` refer to the same respondents;
- rows of `parTab` match the item-column order in `Y`;
- covariates are numeric and have nonzero variance;
- response codes agree with their item types;
- primary covariates in `main` are the intended columns of `X`;
- the selected `lambda` is not constrained by an inadequate grid;
- fitting results are stable across Monte Carlo seeds and larger sample
  sizes;
- `sigma` is positive and convergence diagnostics are acceptable; and
- all plausible values, rather than their respondent-level average, are
  used in secondary analyses.

## Main FARL outputs

| Function | Main output |
|:---|:---|
| [`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md) | Coefficients, residual SD, selected regularization value, factor scores, factor-adjusted designs, item information, and fitting diagnostics. |
| [`Farlr_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md) | Plausible-value matrix, ID-based PV data frame, EAP estimates, posterior variances, posterior modes, and fitted prior means. |
