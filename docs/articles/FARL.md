# FARL: Use Factor-Augmented Regularized Latent Regression for Large-Scale Assessment

## Introduction

The `FARL` package implements factor-augmented regularized latent
regression (FARLR) for large-scale assessments. FARLR is designed for
latent regression with high-dimensional, strongly dependent background
variables. It decomposes their variation into common latent factors and
idiosyncratic components, includes both components in the latent
regression, and uses sparse regularization to select a parsimonious set
of relevant covariates. This structure improves numerical stability
while retaining interpretable background variables and supporting
plausible-value generation for group-level inference.

The current package provides a unified fitting interface for the FARLR
and FARLR-Debias estimators described in the accompanying paper. It
supports dichotomous 2PL and 3PL items, three-category generalized
partial credit model (GPCM) items, high-dimensional covariates, and
plausible-value generation.

The development version can be installed with

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("MAP-LAB-UW/FARL", build_vignettes = TRUE)
torch::install_torch()
```

``` r
library(FARL)
```

[`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md)
is the public interface for both FARLR estimators. The option
`method = "FARLR_EMM"` implements FARLR through an
importance-sampling-based expectation-maximization-maximization
algorithm, whereas `method = "FARLR_Debias"` applies a bias correction
to the regularized coefficients.
[`Dire_mml()`](https://map-lab-uw.github.io/FARL/reference/Dire_mml.md)
provides a PCA-based latent regression comparison.

## Data Input

Data required for the analyses are summarized below.

| Analysis | Item Responses | Item Parameters | Background Covariates | Main Covariates | Formula |
|:--:|:--:|:--:|:--:|:--:|:--:|
| `Farlr_mml(..., method = "FARLR_EMM")` | \checkmark | \checkmark | \checkmark | \checkmark |  |
| `Farlr_mml(..., method = "FARLR_Debias")` | \checkmark | \checkmark | \checkmark | \checkmark |  |
| [`Dire_mml()`](https://map-lab-uw.github.io/FARL/reference/Dire_mml.md) | \checkmark | \checkmark | \checkmark | \checkmark | \checkmark |

Here we first use `sim_a1`, a simulated dataset for one-dimensional 2PL
analysis. It contains N=3000 respondents, J=10 items, and P=60
background covariates.

### Item responses

Item responses should be an N by J numeric matrix. Dichotomous responses
are coded as 0 or 1, and missing responses may be coded as `NA`.

``` r
data(sim_a1)

dim(sim_a1$Y)
# [1] 3000   10
head(sim_a1$Y)
#      i001 i002 i003 i004 i005 i006 i007 i008 i009 i010
# [1,]    1    1    1    1    1    1    0    1    1    1
# [2,]    0    0    1    0    1    0    0    1    1    0
# [3,]    1    1    1    0    0    0    0    0    1    0
# [4,]    0    0    1    0    1    0    0    0    0    0
# [5,]    1    0    1    0    1    1    0    1    0    0
# [6,]    0    0    1    0    1    1    0    1    1    0
```

### Background covariates

Background covariates should be an N by P numeric matrix. Covariates may
be binary or continuous, but each column must have positive variance.

``` r
dim(sim_a1$X)
# [1] 3000   60
head(sim_a1$X[, 1:10])
#      [,1] [,2]       [,3] [,4] [,5] [,6] [,7]       [,8] [,9]      [,10]
# [1,]    1    1  2.0679032    1    1    1    1  1.3037993    1 1.09682610
# [2,]    1    0  0.1913680    0    1    1    0  0.2450981    1 0.63345311
# [3,]    1    1  0.9718647    1    1    1    1  0.1302181    0 0.36494637
# [4,]    0    0 -0.4845058    1    0    0    0  0.1197165    0 0.21706507
# [5,]    1    1  0.8427349    1    1    1    0  0.1999399    1 0.47819559
# [6,]    0    1  0.2229768    1    1    1    1 -0.4628130    1 0.04376757
```

The `main` argument contains the column indices of primary covariates.
These variables are not penalized in the regularized regression. In this
example, the five main covariates are

``` r
main <- c(1, 2, 15, 29, 45)
```

### Item parameters

The item parameter table must follow the same item order as the columns
of the response matrix. For a 2PL item, `a` is the discrimination, `b`
is the difficulty, and `c` is zero.

``` r
sim_a1$parTab

# 1 -0.1280856  1.0005257   0   item1   comp    main    1.0005257   -0.1280856  
# 2 0.4959430   1.1258010   0   item2   comp    main    1.1258010   0.4959430   
# 3 -0.8648743  0.8435464   0   item3   comp    main    0.8435464   -0.8648743  
# 4 0.6850886   1.8747028   0   item4   comp    main    1.8747028   0.6850886   
# 5 -1.4793505  0.8534336   0   item5   comp    main    0.8534336   -1.4793505  
# 6 -0.3336796  1.0371741   0   item6   comp    main    1.0371741   -0.3336796  
# 7 1.8529386   1.4060425   0   item7   comp    main    1.4060425   1.8529386   
# 8 -0.4022468  1.4886061   0   item8   comp    main    1.4886061   -0.4022468  
# 9 0.3655995   1.0044162   0   item9   comp    main    1.0044162   0.3655995   
# 10    0.6513511   1.2750785   0   item10  comp    main    1.2750785   0.6513511   
```

``` text
# Paste the printed sim_a1 dimensions, response rows, covariate rows,
# and parameter table here.
```

### `sim_a2`: mixed 3PL and GPCM responses

The package also includes `sim_a2`, which contains mixed 3PL and GPCM
items. For a 3PL item, `c` is the lower asymptote. A GPCM item is
identified by finite `b1` and `b2` step parameters.

``` r
data(sim_a2)

dim(sim_a2$X)
# [1] 3000   60
dim(sim_a2$Y)
# [1] 3000   10
table(sim_a2$itemtype)
#  3PL gpcm 
#   5    5 
head(sim_a2$Y)
#      i001 i002 i003 i004 i005 i006 i007 i008 i009 i010
# [1,]    1    2    1    2    1    1    1    2    1    1
# [2,]    0    1    0    2    1    2    0    2    1    1
# [3,]    1    0    1    2    1    1    0    0    1    2
# [4,]    0    0    0    2    1    1    0    0    1    0
# [5,]    1    1    1    1    1    0    0    1    1    0
# [6,]    1    1    1    2    1    0    0    2    1    0

sim_a2$parTab[, c("ItemID", "itemtype", "a", "b", "c", "b1", "b2")]
# item1 3PL 1.0005257   -0.1280856  0.1624870   NA  NA
# item2 gpcm    1.1258010   NA  0.0000000   -0.3738297  0.6175084
# item3 3PL 0.8435464   0.4959430   0.2444852   NA  NA
# item4 gpcm    1.8747028   NA  0.0000000   -1.4734702  -0.3204836
# item5 3PL 0.8534336   -0.8648743  0.1599157   NA  NA
# item6 gpcm    1.0371741   NA  0.0000000   0.1658423   1.0450782
# item7 3PL 1.4060425   0.6850886   0.1887100   NA  NA
# item8 gpcm    1.4886061   NA  0.0000000   -0.7241589  0.4505988
# item9 3PL 1.0044162   -1.4793505  0.1994257   NA  NA
# item10    gpcm    1.2750785   NA  0.0000000   -0.2383690  0.8670740
```

``` text
# Paste the sim_a2 item-type table and mixed-item parameter table here.
```

## Data Output

The main model and plausible-value functions return the following
objects.

| Function | Main output |
|:--:|:---|
| [`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md) | Regression coefficients, residual standard deviation, selected tuning parameter, factor scores, and method-specific design matrices |
| [`Farlr_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md) | Plausible values, EAP estimates, posterior variances, posterior modes, and prior means |
| [`Dire_mml()`](https://map-lab-uw.github.io/FARL/reference/Dire_mml.md) | A fitted `mmlMeans` object with DIRE coefficients and the residual-PCA object |
| [`Dire_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Dire_drawPVs.md) | A data frame containing respondent IDs and plausible values |

## Factor-Augmented Regularized Latent Regression

[`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md)
is the unified interface for the two FARLR estimators. In both methods,
estimated factors capture shared dependence among covariates, while
regularization identifies a sparse set of relevant idiosyncratic
predictors. Variables listed in `main` remain unpenalized so that key
reporting variables are retained in the population model.

| Method | Description |
|:--:|:---|
| `FARLR_EMM` | Importance-sampling E-step, regularized variable-selection M-step, and an additional unpenalized M-step for the active coefficients |
| `FARLR_Debias` | Regularized estimation using factors and idiosyncratic components, followed by correction of LASSO shrinkage bias |

When `K_hat` is not supplied, the number of factors is estimated by
parallel analysis. If the number of factors is known, supplying `K_hat`
avoids this additional step.

### FARLR via ISEMM

This method combines the estimated factor scores and observed covariates
in the latent regression. At each iteration, importance sampling
approximates the posterior expectations of latent proficiency, a LASSO
update selects relevant covariates, and a second unpenalized M-step
refines the active coefficients.

``` r
mmlcomp_emm <- with(
  sim_a1,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_EMM",
    main = main,
    K_hat = NULL,
    seed = 2026
  )
)

mmlcomp_emm$coefficients
mmlcomp_emm$sigma
```

### FARLR Debias

The debias method separates the common factor structure from the
idiosyncratic component of each covariate. It first obtains regularized
coefficient estimates and then corrects part of the LASSO shrinkage
bias. This decorrelated parameterization is intended to stabilize
selection when background variables are strongly dependent.

``` r
mmlcomp_debias <- with(
  sim_a1,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_Debias",
    main = main,
    K_hat = 2,
    seed = 2026
  )
)

mmlcomp_debias$coefficients
mmlcomp_debias$sigma
```

### Mixed 3PL and GPCM items

The same interface can analyze the mixed-format responses in `sim_a2`.
[`Farlr_mml()`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md)
determines item types from `parTab`: finite `b1` and `b2` identify GPCM
items, and a nonzero `c` identifies 3PL items.

``` r
mmlcomp_mixed <- with(
  sim_a2,
  Farlr_mml(
    X = X,
    Y = Y,
    parTab = parTab,
    method = "FARLR_Debias",
    main = main,
    K_hat = 2,
    seed = 2026
  )
)

mmlcomp_mixed$coefficients
mmlcomp_mixed$sigma
```

## Plausible Values for FARLR

[`Farlr_drawPVs()`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md)
combines the fitted latent-regression prior with each respondent’s
item-response likelihood to obtain a posterior latent-trait
distribution. The default normal method draws from a normal
approximation based on the EAP mean and variance. The grid method
samples directly from the discretized posterior.

``` r
PVs_debias <- Farlr_drawPVs(
  mmlcomp = mmlcomp_debias,
  npv = 10L,
  draw_method = "normal",
  seed = 2026
)

head(PVs_debias$datPVs)
head(PVs_debias$EAP_estimates)
```

``` text
# Paste head(PVs_debias$datPVs) and head(PVs_debias$EAP_estimates) here.
```

For direct grid-based posterior draws, use

``` r
PVs_grid <- Farlr_drawPVs(
  mmlcomp = mmlcomp_debias,
  npv = 10L,
  draw_method = "grid",
  seed = 2026
)
```

### Plausible values for `sim_a2`

Plausible values for the mixed 3PL/GPCM data are drawn with the same
function. The fitted model retains the item-specific response model
supplied through `sim_a2$parTab`.

``` r
PVs_mixed <- Farlr_drawPVs(
  mmlcomp = mmlcomp_mixed,
  npv = 10L,
  draw_method = "normal",
  seed = 2026
)

head(PVs_mixed$datPVs)
head(PVs_mixed$EAP_estimates)
```

``` text
# Paste head(PVs_mixed$datPVs) and head(PVs_mixed$EAP_estimates) here.
```

## DIRE Marginal Maximum Likelihood

The PCA-based latent regression commonly used in large-scale assessments
reduces high-dimensional background variables to a selected set of
principal components. In FARL,
[`Dire_mml()`](https://map-lab-uw.github.io/FARL/reference/Dire_mml.md)
implements this comparison through the DIRE package. It can apply PCA to
the original covariates or to non-main covariates after residualizing
them on the key predictors, then adds the retained PC scores to the
latent regression formula.

### Prepare the DIRE data

``` r
library(Dire)

X <- as.matrix(sim_a1$X)
resp <- as.matrix(sim_a1$Y)
n <- nrow(X)
J <- ncol(resp)

colnames(X) <- paste0("X", seq_len(ncol(X)))
main_vars <- paste0("X", main)

subject <- factor(seq_len(n))
stuDat <- data.frame(subject = subject, X)

item_names <- paste0("item", seq_len(J))
stuItems <- data.frame(
  subject = rep(subject, times = J),
  key = factor(rep(item_names, each = n), levels = item_names),
  score = as.vector(resp)
)

parTab <- sim_a1$parTab[, c(
  "ItemID", "test", "subtest", "slope",
  "difficulty", "guessing", "D"
)]
parTab$ItemID <- item_names

testDat <- data.frame(
  test = "comp",
  subtest = "main",
  location = 0,
  scale = 1
)
```

### Fit the model

The `pca_type` argument chooses covariance or correlation PCA, and
`var_threshold` specifies the cumulative proportion of variance
retained.

``` r
mmlcomp_dire <- Dire_mml(
  formula = reformulate(main_vars, response = "comp"),
  stuItems = stuItems,
  stuDat = stuDat,
  idVar = "subject",
  dichotParamTab = parTab,
  testScale = testDat,
  X = X,
  main_vars = main_vars,
  pca_type = "cov",
  use_residual = TRUE,
  var_threshold = 0.90
)

mmlcomp_dire$coefficients
mmlcomp_dire$pca_object$num_pcs
```

### Draw DIRE plausible values

``` r
PVs_dire <- Dire_drawPVs(
  x = mmlcomp_dire,
  npv = 10L,
  pvVariableNameSuffix = "_dire"
)$data

# DIRE may sort character IDs internally. Restore the original respondent order.
pv_order <- match(
  as.character(seq_len(n)),
  as.character(PVs_dire$id)
)

PVs_dire <- PVs_dire[pv_order, , drop = FALSE]
rownames(PVs_dire) <- NULL

head(PVs_dire)
```

``` text
# Paste head(PVs_dire) here.
```

## Package Evaluation

The following short functions can be used to test the primary FARL
workflows. Each function fits a model and returns the plausible-value
result.

### FARLR Debias test

``` r
mml_test_debias <- function() {
  mmlcomp <- with(
    sim_a1,
    Farlr_mml(
      X,
      Y,
      parTab,
      method = "FARLR_Debias",
      main = c(1, 2, 15, 29, 45),
      K_hat = 2,
      seed = 2026
    )
  )

  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp,
    npv = 10L,
    seed = 2026
  )

  return(PVs)
}
```

### FARLR EMM test

``` r
mml_test_emm <- function() {
  mmlcomp <- with(
    sim_a1,
    Farlr_mml(
      X,
      Y,
      parTab,
      method = "FARLR_EMM",
      main = c(1, 2, 15, 29, 45),
      K_hat = 2,
      seed = 2026
    )
  )

  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp,
    npv = 10L,
    seed = 2026
  )

  return(PVs)
}
```

### Mixed-item test

``` r
mml_test_mixed <- function() {
  mmlcomp <- with(
    sim_a2,
    Farlr_mml(
      X,
      Y,
      parTab,
      method = "FARLR_Debias",
      main = c(1, 2, 15, 29, 45),
      K_hat = 2,
      seed = 2026
    )
  )

  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp,
    npv = 10L,
    seed = 2026
  )

  return(PVs)
}
```

## Notes

- FARLR is intended to improve compatibility with common secondary
  analyses; it does not guarantee congeniality in the strict statistical
  sense.
- The current FARLR estimators are designed for a unidimensional latent
  proficiency model.
- Rows of `X` and `Y` must correspond to the same respondents.
- Rows of `parTab` must follow the item-column order in `Y`.
- `main` uses numeric column indices from `X`.
- Set `seed` for reproducible FARLR estimation and plausible values.
- Set `K_hat` when the number of factors is known; otherwise parallel
  analysis is used.
- For long-running package vignettes, prebuild the document or keep the
  model chunks unevaluated during routine package checks.
