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

#' Fit a Multidimensional DIRE Latent Regression
#'
#' Fits a multidimensional DIRE model with one subtest per latent dimension.
#' Main covariates enter the latent regression explicitly, while residual or
#' raw principal components summarize the remaining high-dimensional
#' covariate information.
#'
#' @param X `N` by `p` numeric covariate matrix.
#' @param Y `N` by `J_total` binary response matrix.
#' @param parTab Item parameter data frame with one row per item. For the
#'   simple-structure format used by `sim_m1`, it must contain `a` (or `slope`),
#'   `d` (or `b`/`difficulty`), `c` (or `guessing`), and `dimension` (or
#'   `subtest`).
#' @param main Character names or numeric column indices of the main
#'   covariates retained explicitly in the latent regression.
#' @param subtest_names Optional character vector of `D` subtest names.
#' @param test_name Name of the composite test passed to DIRE.
#' @param pca_type Use covariance (`"cov"`) or correlation (`"cor"`) PCA.
#' @param use_residual If `TRUE`, perform PCA on covariate residuals after
#'   conditioning on the main variables.
#' @param var_threshold Cumulative explained-variance threshold used to retain
#'   principal components.
#' @param calcCor Logical passed to [Dire::mml()].
#' @param ... Additional arguments passed to [Dire::mml()].
#'
#' @return An object of class `MDIREfit` containing the fitted DIRE object and
#'   all data structures needed by [MDIRE_drawPVs()].
#'
#' @examples
#' \dontrun{
#' fit <- MDIRE_mml(
#'   X = sim_m1$X,
#'   Y = sim_m1$Y,
#'   parTab = sim_m1$parTab,
#'   main = sim_m1$main
#' )
#' }
#'
#' @export
MDIRE_mml <- function(
    X,
    Y,
    parTab,
    main,
    subtest_names = NULL,
    test_name = "comp",
    pca_type = c("cov", "cor"),
    use_residual = TRUE,
    var_threshold = 0.90,
    calcCor = TRUE,
    ...
) {
  pca_type <- match.arg(pca_type)

  if (!requireNamespace("Dire", quietly = TRUE)) {
    stop("Package 'Dire' is required.", call. = FALSE)
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  input_parTab <- parTab

  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.", call. = FALSE)
  }
  suppressWarnings(storage.mode(Y) <- "numeric")
  if (any(!is.finite(Y[!is.na(Y)]))) {
    stop("Observed responses must be finite numeric values.", call. = FALSE)
  }

  N <- nrow(X)
  p <- ncol(X)
  J_total <- ncol(Y)

  if (nrow(Y) != N) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }

  parameter_function <- get0(
    "multidim_item_parameters",
    mode = "function",
    inherits = TRUE
  )
  if (is.null(parameter_function)) {
    stop("Internal function multidim_item_parameters() was not found.", call. = FALSE)
  }

  item_parameters <- parameter_function(
    parTab = parTab,
    J = J_total,
    allow_guessing = TRUE
  )
  a <- item_parameters$a
  d <- item_parameters$d
  c <- item_parameters$c
  item_index <- item_parameters$item_index
  D_latent <- ncol(a)

  observed_response <- Y[!is.na(Y)]
  if (any(!observed_response %in% c(0, 1))) {
    stop("Y must contain only 0, 1, or NA.", call. = FALSE)
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  if (anyDuplicated(colnames(X))) {
    stop("X column names must be unique.", call. = FALSE)
  }

  if (is.numeric(main)) {
    main_indices <- unique(as.integer(main))
    if (anyNA(main_indices) ||
        any(main_indices < 1L | main_indices > p)) {
      stop("Numeric main must index columns of X.", call. = FALSE)
    }
    main_vars <- colnames(X)[main_indices]
  } else {
    main_vars <- unique(as.character(main))
    if (length(main_vars) == 0L || any(!main_vars %in% colnames(X))) {
      stop("main must name columns of X.", call. = FALSE)
    }
  }

  # -----------------------------------------------------------------------
  # Item-to-dimension assignment
  # -----------------------------------------------------------------------
  nonzero_loading <- a != 0
  if (any(rowSums(nonzero_loading) != 1L)) {
    stop(
      "MDIRE_mml() requires each item in parTab to load on exactly one ",
      "dimension.",
      call. = FALSE
    )
  }
  item_dimension <- max.col(nonzero_loading, ties.method = "first")

  item_slope <- a[cbind(seq_len(J_total), item_dimension)]
  if (any(!is.finite(item_slope)) || any(item_slope == 0)) {
    stop("Every assigned item must have a finite nonzero slope.", call. = FALSE)
  }

  b <- -d / item_slope

  if (is.null(subtest_names)) {
    subtest_names <- if (D_latent <= length(LETTERS)) {
      LETTERS[seq_len(D_latent)]
    } else {
      paste0("Dimension", seq_len(D_latent))
    }
  }
  subtest_names <- as.character(subtest_names)
  if (length(subtest_names) != D_latent ||
      anyNA(subtest_names) || anyDuplicated(subtest_names)) {
    stop("subtest_names must contain D unique names.", call. = FALSE)
  }
  names(item_index) <- subtest_names

  # -----------------------------------------------------------------------
  # PCA covariates
  # -----------------------------------------------------------------------
  pca_function <- get0(
    "PCA_extraction_merge",
    mode = "function",
    inherits = TRUE
  )
  if (is.null(pca_function)) {
    stop(
      "Internal function PCA_extraction_merge() was not found.",
      call. = FALSE
    )
  }

  pca_object <- pca_function(
    main_vars = main_vars,
    X = X,
    pca_type = pca_type,
    use_residual = use_residual,
    var_threshold = var_threshold
  )

  pc_scores <- pca_object$pc_scores
  if (is.null(pc_scores)) {
    pc_scores <- pca_object$pc_score
  }
  if (is.null(pc_scores)) {
    stop("PCA_extraction_merge() did not return PC scores.", call. = FALSE)
  }
  pc_scores <- as.matrix(pc_scores)
  if (nrow(pc_scores) != N) {
    stop("PCA scores and X have different row counts.", call. = FALSE)
  }
  if (ncol(pc_scores) > 0L) {
    colnames(pc_scores) <- paste0("PC", seq_len(ncol(pc_scores)))
  }

  # -----------------------------------------------------------------------
  # DIRE data objects
  # -----------------------------------------------------------------------
  subject <- factor(seq_len(N))
  stuDat <- data.frame(
    subject = subject,
    X,
    check.names = FALSE
  )
  if (ncol(pc_scores) > 0L) {
    stuDat <- cbind(stuDat, as.data.frame(pc_scores, check.names = FALSE))
  }

  item_within_dimension <- integer(J_total)
  for (dimension in seq_len(D_latent)) {
    item_within_dimension[item_index[[dimension]]] <-
      seq_along(item_index[[dimension]])
  }
  item_ids <- paste0(
    "item", subtest_names[item_dimension], item_within_dimension
  )

  stuItems <- data.frame(
    subject = rep(subject, times = J_total),
    key = factor(rep(item_ids, each = N), levels = item_ids),
    score = as.vector(Y),
    row.names = NULL
  )

  testDat <- data.frame(
    test = rep(test_name, D_latent),
    subtest = subtest_names,
    location = rep(0, D_latent),
    scale = rep(1, D_latent),
    subtestWeight = rep(1 / D_latent, D_latent),
    stringsAsFactors = FALSE
  )

  dire_parTab <- data.frame(
    ItemID = item_ids,
    test = rep(test_name, J_total),
    subtest = subtest_names[item_dimension],
    slope = item_slope,
    difficulty = b,
    guessing = c,
    D = 1,
    stringsAsFactors = FALSE
  )

  predictors <- c(main_vars, colnames(pc_scores))
  predictors <- predictors[nzchar(predictors)]
  if (length(predictors) == 0L) {
    stop("The latent regression requires at least one predictor.", call. = FALSE)
  }

  formula <- stats::reformulate(
    termlabels = predictors,
    response = test_name,
    intercept = FALSE
  )

  start_time <- Sys.time()
  mmlcomp <- Dire::mml(
    formula,
    stuItems = stuItems,
    stuDat = stuDat,
    idVar = "subject",
    dichotParamTab = dire_parTab,
    testScale = testDat,
    calcCor = calcCor,
    ...
  )
  elapsed <- Sys.time() - start_time

  result <- list(
    mmlcomp = mmlcomp,
    coefficients = mmlcomp$coefficients,
    sigma = mmlcomp$SubscaleVC,
    formula = formula,
    pca_object = pca_object,
    pc_scores = pc_scores,
    stuDat = stuDat,
    stuItems = stuItems,
    parTab = input_parTab,
    dire_parTab = dire_parTab,
    testDat = testDat,
    X = X,
    Y = Y,
    response = Y,
    a = a,
    b = b,
    d = d,
    c = c,
    main_vars = main_vars,
    item_index = item_index,
    subtest_names = subtest_names,
    test_name = test_name,
    D = D_latent,
    elapsed = elapsed,
    call = match.call()
  )
  class(result) <- c("MDIREfit", "list")
  result
}


#' Draw Plausible Values from a Multidimensional DIRE Model
#'
#' Draws plausible values from an object fitted by [MDIRE_mml()]. The function
#' identifies the PV columns for every subtest without assuming five dimensions
#' or ten plausible values. Optional binary grouping variables can be supplied
#' to reproduce the subgroup summaries from the original `highDimDire()` code.
#'
#' @param object An object returned by [MDIRE_mml()].
#' @param npv Positive integer. Number of plausible values per respondent.
#' @param pvVariableNameSuffix Suffix passed to [Dire::drawPVs()].
#' @param group_vars Optional character names or numeric indices of binary
#'   columns in `object$X` used to form subgroup combinations.
#' @param ... Additional arguments passed to [Dire::drawPVs()].
#'
#' @return A list containing raw PV data, dimension-specific PV matrices,
#'   optional subgroup PV matrices and summaries, and the original DIRE PV
#'   result.
#'
#' @export
MDIRE_drawPVs <- function(
    object,
    npv = 10L,
    pvVariableNameSuffix = "_dire",
    group_vars = NULL,
    ...
) {
  if (!inherits(object, "MDIREfit")) {
    stop("object must be returned by MDIRE_mml().", call. = FALSE)
  }
  if (!requireNamespace("Dire", quietly = TRUE)) {
    stop("Package 'Dire' is required.", call. = FALSE)
  }

  npv <- as.integer(npv)
  if (length(npv) != 1L || is.na(npv) || npv < 1L) {
    stop("npv must be a positive integer.", call. = FALSE)
  }

  pv_result <- Dire::drawPVs(
    x = object$mmlcomp,
    npv = npv,
    pvVariableNameSuffix = pvVariableNameSuffix,
    ...
  )
  datPVs <- pv_result$data
  if (!is.data.frame(datPVs)) {
    datPVs <- as.data.frame(datPVs, check.names = FALSE)
  }
  if (!"id" %in% names(datPVs)) {
    stop("Dire::drawPVs() did not return an id column.", call. = FALSE)
  }

  D_latent <- object$D
  subtest_names <- object$subtest_names
  dimension_pvs <- vector("list", D_latent)
  names(dimension_pvs) <- subtest_names

  for (dimension in seq_len(D_latent)) {
    expected_names <- paste0(
      subtest_names[dimension],
      pvVariableNameSuffix,
      seq_len(npv)
    )
    missing_names <- setdiff(expected_names, names(datPVs))
    if (length(missing_names) > 0L) {
      stop(
        "Missing PV columns for subtest ", subtest_names[dimension], ": ",
        paste(missing_names, collapse = ", "),
        call. = FALSE
      )
    }
    dimension_pvs[[dimension]] <- as.matrix(
      datPVs[, expected_names, drop = FALSE]
    )
    colnames(dimension_pvs[[dimension]]) <- paste0("PV", seq_len(npv))
  }

  # Helpers equivalent to variance.sample() and variance.cal() in the
  # original implementation.
  sampling_variance <- function(pv_matrix) {
    if (nrow(pv_matrix) < 2L) {
      return(NA_real_)
    }
    mean(apply(pv_matrix, 2L, stats::var))
  }

  total_variance <- function(pv_matrix) {
    if (nrow(pv_matrix) < 2L) {
      return(NA_real_)
    }
    within <- sampling_variance(pv_matrix)
    between <- if (ncol(pv_matrix) > 1L) {
      stats::var(colMeans(pv_matrix))
    } else {
      0
    }
    within + (1 + 1 / ncol(pv_matrix)) * between
  }

  X <- object$X
  subject_ids <- as.character(object$stuDat$subject)
  pv_subject_ids <- as.character(datPVs$id)
  pv_to_X <- match(pv_subject_ids, subject_ids)
  if (anyNA(pv_to_X)) {
    stop("Could not align PV ids with the fitted subject data.", call. = FALSE)
  }

  if (is.null(group_vars)) {
    combinations <- data.frame(group = "overall")
    group_membership <- list(rep(TRUE, nrow(datPVs)))
    group_names <- "overall"
  } else {
    if (is.numeric(group_vars)) {
      group_indices <- unique(as.integer(group_vars))
      if (anyNA(group_indices) ||
          any(group_indices < 1L | group_indices > ncol(X))) {
        stop("Numeric group_vars must index columns of X.", call. = FALSE)
      }
      group_vars <- colnames(X)[group_indices]
    } else {
      group_vars <- unique(as.character(group_vars))
      if (any(!group_vars %in% colnames(X))) {
        stop("group_vars must name columns of object$X.", call. = FALSE)
      }
    }

    group_data <- X[pv_to_X, group_vars, drop = FALSE]
    if (any(!group_data %in% c(0, 1))) {
      stop("Every group_vars column must be binary and coded 0/1.", call. = FALSE)
    }

    combinations <- expand.grid(
      replicate(length(group_vars), 0:1, simplify = FALSE)
    )
    names(combinations) <- group_vars
    group_names <- apply(
      combinations,
      1L,
      function(values) paste0("result_", paste(values, collapse = ""))
    )

    group_membership <- lapply(
      seq_len(nrow(combinations)),
      function(group_index) {
        Reduce(
          `&`,
          lapply(
            seq_along(group_vars),
            function(variable_index) {
              group_data[, variable_index] ==
                combinations[group_index, variable_index]
            }
          )
        )
      }
    )
  }

  group_pvs <- lapply(
    dimension_pvs,
    function(dimension_matrix) {
      result <- lapply(
        group_membership,
        function(membership) dimension_matrix[membership, , drop = FALSE]
      )
      names(result) <- group_names
      result
    }
  )

  group_summary <- do.call(
    rbind,
    lapply(seq_len(D_latent), function(dimension) {
      do.call(
        rbind,
        lapply(seq_along(group_names), function(group_index) {
          values <- group_pvs[[dimension]][[group_index]]
          data.frame(
            dimension = subtest_names[dimension],
            group = group_names[group_index],
            sample_size = nrow(values),
            mean = if (length(values)) mean(values) else NA_real_,
            var = sampling_variance(values),
            mean_v = total_variance(values),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  rownames(group_summary) <- NULL

  list(
    data = datPVs,
    rawPV = datPVs,
    dimension_pvs = dimension_pvs,
    pv = group_pvs,
    combinations = combinations,
    group_summary = group_summary,
    mean = group_summary$mean,
    var = group_summary$var,
    mean_v = group_summary$mean_v,
    sample_size = group_summary$sample_size,
    mmlcomp = object$mmlcomp,
    pv_result = pv_result,
    call = match.call()
  )
}
