library(Dire)


PCA_extraction_merge <- function(main_vars, X,
                                 pca_type = c("cov", "cor"),
                                 use_residual = TRUE,
                                 var_threshold = 0.90) {
  pca_type <- match.arg(pca_type)

  if (!is.numeric(var_threshold) || length(var_threshold) != 1L ||
      is.na(var_threshold) || var_threshold <= 0 || var_threshold > 1) {
    stop("`var_threshold` must be one number in (0, 1].", call. = FALSE)
  }

  original_names <- colnames(X)
  X <- as.data.frame(X)
  if (is.null(original_names) || any(original_names == "")) {
    names(X) <- paste0("X", seq_len(ncol(X)))
  }

  if (anyDuplicated(names(X))) {
    stop("`X` must have unique column names.", call. = FALSE)
  }
  if (anyNA(X) || any(!vapply(X, is.numeric, logical(1)))) {
    stop("`X` must contain only numeric values without missing values.",
         call. = FALSE)
  }
  if (!length(main_vars) || anyDuplicated(main_vars) ||
      !all(main_vars %in% names(X))) {
    stop("Every value in `main_vars` must uniquely name a column of `X`.",
         call. = FALSE)
  }

  scale_flag <- identical(pca_type, "cor")

  if (use_residual) {
    remaining_vars <- setdiff(names(X), main_vars)
    if (!length(remaining_vars)) {
      stop("Residual PCA requires at least one non-main variable.",
           call. = FALSE)
    }

    residual_formula <- reformulate(main_vars)
    design_matrix <- model.matrix(residual_formula, data = X)
    residual_fit <- lm.fit(
      x = design_matrix,
      y = as.matrix(X[, remaining_vars, drop = FALSE])
    )

    residuals_matrix <- as.matrix(residual_fit$residuals)
    coefficients_matrix <- as.matrix(residual_fit$coefficients)
    colnames(residuals_matrix) <- remaining_vars
    colnames(coefficients_matrix) <- remaining_vars
    rownames(coefficients_matrix) <- colnames(design_matrix)
    pca_input <- residuals_matrix
  } else {
    remaining_vars <- names(X)
    pca_input <- as.matrix(X)
    residuals_matrix <- NULL
    coefficients_matrix <- NULL
  }

  finite_nonconstant <- vapply(
    seq_len(ncol(pca_input)),
    function(j) {
      values <- pca_input[, j]
      all(is.finite(values)) && stats::sd(values) > sqrt(.Machine$double.eps)
    },
    logical(1)
  )
  dropped_vars <- colnames(pca_input)[!finite_nonconstant]
  if (length(dropped_vars)) {
    warning(
      "Dropping non-finite or effectively constant PCA columns: ",
      paste(dropped_vars, collapse = ", "),
      call. = FALSE
    )
    pca_input <- pca_input[, finite_nonconstant, drop = FALSE]
  }
  if (!ncol(pca_input)) {
    stop("No usable columns remain for PCA.", call. = FALSE)
  }

  pca_res <- stats::prcomp(
    pca_input,
    center = TRUE,
    scale. = scale_flag
  )
  explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cumulative_var <- cumsum(explained_var)
  num_pcs <- which(cumulative_var >= var_threshold)[1L]
  pc_scores <- pca_res$x[, seq_len(num_pcs), drop = FALSE]
  colnames(pc_scores) <- paste0("PC", seq_len(num_pcs))
  names(explained_var) <- paste0("PC", seq_along(explained_var))
  names(cumulative_var) <- names(explained_var)

  list(
    pc_scores = pc_scores,
    num_pcs = num_pcs,
    explained_var = explained_var,
    cumulative_var = cumulative_var,
    residuals_matrix = residuals_matrix,
    coefficients_matrix = coefficients_matrix,
    pca_input = pca_input,
    pca_res = pca_res,
    dropped_vars = dropped_vars,
    main_vars = main_vars,
    pca_type = pca_type,
    use_residual = use_residual,
    var_threshold = var_threshold
  )
}


#' DIRE Marginal Maximum Likelihood Estimation with Optional PCA
#'
#' Fits a DIRE model with `Dire::mml()`. When `X` and `main_vars` are supplied,
#' the function first extracts PCs with `PCA_extraction_merge()`, merges the PC
#' scores into `stuDat`, appends the PC terms to `formula`, and then generates
#' the fitted `mmlcomp` object.
#'
#' @param formula Formula passed to `Dire::mml()`. It should already contain
#'   the main covariates that are to remain explicit in the latent regression.
#' @param stuItems Long-format item response data.
#' @param stuDat Person-level data. Its row order must match the row order of
#'   `X` when PCA is requested.
#' @param idVar Character name of the individual identifier.
#' @param dichotParamTab Dichotomous item parameter table or `NULL`.
#' @param polyParamTab Polytomous item parameter table or `NULL`.
#' @param testScale Optional test-scale table.
#' @param Q Number of quadrature nodes.
#' @param minNode,maxNode Lower and upper quadrature bounds.
#' @param polyModel Polytomous model, `"GPCM"` or `"GRM"`.
#' @param weightVar Optional sampling-weight variable.
#' @param multiCore Whether to use supported parallel calculations.
#' @param bobyqaControl Optional optimizer controls.
#' @param composite Whether to use composite likelihood.
#' @param strataVar,PSUVar Optional complex-survey variables.
#' @param fast Whether to use DIRE computational shortcuts.
#' @param calcCor Whether to calculate the coefficient correlation matrix.
#' @param verbose DIRE verbosity level.
#' @param X Optional numeric covariate matrix/data frame used for PCA. Leave
#'   `NULL` to fit the original non-PCA DIRE model.
#' @param main_vars Character vector naming the main variables in `X`.
#' @param pca_type PCA based on covariance (`"cov"`) or correlation (`"cor"`).
#' @param use_residual Whether to perform residual PCA.
#' @param var_threshold Cumulative explained-variance threshold for selecting
#'   PCs.
#'
#' @return The fitted DIRE `mmlcomp` object. When PCA is used, the returned
#'   object also contains `pca_object` and `formula_before_pca` components; its
#'   standard `formula` and `stuDat` components contain the augmented versions.
#' @export
Dire_mml <- function(formula, stuItems, stuDat, idVar,
                     dichotParamTab = NULL, polyParamTab = NULL,
                     testScale = NULL, Q = 30, minNode = -4, maxNode = 4,
                     polyModel = c("GPCM", "GRM"), weightVar = NULL,
                     multiCore = FALSE, bobyqaControl = NULL,
                     composite = TRUE, strataVar = NULL, PSUVar = NULL,
                     fast = TRUE, calcCor = TRUE, verbose = 0,
                     X = NULL, main_vars = NULL,
                     pca_type = c("cov", "cor"),
                     use_residual = TRUE, var_threshold = 0.90) {
  polyModel <- match.arg(polyModel)
  pca_type <- match.arg(pca_type)
  formula_before_pca <- formula
  pca_object <- NULL

  if (is.null(X) != is.null(main_vars)) {
    stop("Supply both `X` and `main_vars` to use PCA, or leave both `NULL`.",
         call. = FALSE)
  }

  if (!is.null(X)) {
    if (nrow(X) != nrow(stuDat)) {
      stop("`X` and `stuDat` must contain the same individuals in the same order.",
           call. = FALSE)
    }

    pca_object <- PCA_extraction_merge(
      main_vars = main_vars,
      X = X,
      pca_type = pca_type,
      use_residual = use_residual,
      var_threshold = var_threshold
    )
    pc_scores <- as.data.frame(pca_object$pc_scores)

    duplicate_pc_names <- intersect(names(pc_scores), names(stuDat))
    if (length(duplicate_pc_names)) {
      stop(
        "PC columns already exist in `stuDat`: ",
        paste(duplicate_pc_names, collapse = ", "),
        call. = FALSE
      )
    }
    stuDat <- cbind(stuDat, pc_scores)

    pc_formula <- paste(names(pc_scores), collapse = " + ")
    formula <- stats::update.formula(
      formula,
      stats::as.formula(paste(". ~ . +", pc_formula))
    )

    if (verbose >= 1) {
      message(
        "PCA retained ", pca_object$num_pcs, " component(s), explaining ",
        round(100 * pca_object$cumulative_var[pca_object$num_pcs], 2),
        "% of PCA-input variance."
      )
    }
  }

  mmlcomp <- Dire::mml(
    formula = formula,
    stuItems = stuItems,
    stuDat = stuDat,
    idVar = idVar,
    dichotParamTab = dichotParamTab,
    polyParamTab = polyParamTab,
    testScale = testScale,
    Q = Q,
    minNode = minNode,
    maxNode = maxNode,
    polyModel = polyModel,
    weightVar = weightVar,
    multiCore = multiCore,
    bobyqaControl = bobyqaControl,
    composite = composite,
    strataVar = strataVar,
    PSUVar = PSUVar,
    fast = fast,
    calcCor = calcCor,
    verbose = verbose
  )

  if (!is.null(pca_object)) {
    mmlcomp$pca_object <- pca_object
    mmlcomp$formula_before_pca <- formula_before_pca
  }

  mmlcomp
}


#' Draw Plausible Values from a Fitted DIRE Model
#'
#' This is the second step of the workflow. It accepts the `mmlcomp` object
#' generated by `Dire_mml()` and delegates plausible-value generation to
#' `Dire::drawPVs()`.
#'
#' @param x A fitted DIRE `mmlcomp` object returned by `Dire_mml()`.
#' @param npv Positive integer number of plausible values per individual.
#' @param pvVariableNameSuffix Suffix appended to plausible-value column names.
#' @param ... Additional arguments passed to `Dire::drawPVs()`.
#'
#' @return The object returned by `Dire::drawPVs()`.
#' @export
Dire_drawPVs <- function(x, npv,
                         pvVariableNameSuffix = "_dire", ...) {
  if (!inherits(x, "mmlMeans")) {
    stop("`x` must be the fitted mmlcomp object returned by `Dire_mml()`.",
         call. = FALSE)
  }
  if (length(npv) != 1L || is.na(npv) || npv < 1 || npv != as.integer(npv)) {
    stop("`npv` must be a positive integer.", call. = FALSE)
  }

  Dire::drawPVs(
    x = x,
    npv = as.integer(npv),
    pvVariableNameSuffix = pvVariableNameSuffix,
    ...
  )
}


# Example workflow:
#
# mmlcomp <- Dire_mml(
#   formula = comp ~ X1 + X2 + X29 + X15 + X45,
#   stuItems = stuItems,
#   stuDat = stuDat,
#   idVar = "subject",
#   dichotParamTab = parTab,
#   testScale = testDat,
#   X = X,
#   main_vars = main_vars,
#   pca_type = "cov",
#   use_residual = TRUE,
#   var_threshold = 0.90,
#   calcCor = TRUE
# )
#
# datPVs <- Dire_drawPVs(mmlcomp, npv = 10L)$data
