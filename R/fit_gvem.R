library(paran)
library(lavaan)
library(mirt)
library(glmnet)
library(torch)
#' FARLR Marginal Maximum Likelihood Estimation
#'
#' Fits a Factor-Adjusted Regularized Latent Regression (FARLR) model using a
#' marginal maximum likelihood (MML) framework. This function serves as a unified
#' interface that supports multiple estimation strategies, including
#' \code{"FARLR_EMM"} and \code{"FARLR_Debias"}. Latent traits are integrated out
#' using Monte Carlo approximation, and regression parameters are estimated under
#' regularization. Depending on the specified method, the algorithm either
#' employs an EM–M–type iterative scheme or a post-selection debiasing procedure.
#'
#' @param X Matrix. Covariate design matrix for the latent regression model,
#'   typically of dimension \code{n x p}.
#' @param Y Matrix. Observed item response matrix of dimension \code{n x J}.
#' @param parTab Data frame. Item parameter table containing at least the columns
#'   \code{slope}, \code{difficulty}, and \code{guessin}, used to define the item
#'   response model.
#' @param n_sam Integer. Number of Monte Carlo samples per individual used to
#'   approximate integrals over latent traits. Defaults to \code{5}.
#' @param method Character string. Estimation method to be used. Supported values
#'   include:
#'   \describe{
#'     \item{\code{"FARLR_EMM"}}{Iterative FARLR estimation based on an EM–M–type
#'     updating scheme with regularization.}
#'     \item{\code{"FARLR_Debias"}}{FARLR estimation with regularized variable
#'     selection followed by a post-selection debiased refit.}
#'   }
#'   Defaults to \code{"FARLR_EMM"}.
#' @param lambda Numeric vector. Candidate regularization parameters used for
#'   penalized regression. Defaults to \code{seq(0.1, 0.5, by = 0.1)}.
#' @param delta.criteria Numeric. Convergence tolerance for iterative updates.
#'   Defaults to \code{1e-3}.
#' @param iter.max Integer. Maximum number of iterations allowed for each value of
#'   \code{lambda}. Defaults to \code{200}.
#' @param window.size Integer. Window size for sliding-window averaging of
#'   coefficient updates used to stabilize iterative estimation. Defaults to
#'   \code{50}.
#' @param verbose Logical. If \code{TRUE}, progress messages and iteration status
#'   are displayed during model fitting. Defaults to \code{TRUE}.
#'
#' @return A list containing estimation results. The exact contents depend on the
#'   selected \code{method}, but typically include:
#' \describe{
#'   \item{\code{coefficients}}{Estimated regression coefficients.}
#'   \item{\code{sigma}}{Estimated residual standard deviation.}
#'   \item{\code{LogLik}}{Value of the objective function evaluated at the selected
#'   regularization parameter.}
#'   \item{\code{lambda}}{Selected regularization parameter.}
#'   \item{\code{Convergence}}{Indicator of convergence status.}
#' }
#'
#' @details
#' The function marginalizes over latent variables using Monte Carlo integration
#' and estimates regression parameters under regularization. When
#' \code{method = "FARLR_EMM"}, parameters are updated iteratively using an
#' EM–M–style procedure. When \code{method = "FARLR_Debias"}, a penalized estimator
#' is first used for variable selection, followed by a debiased refit on the
#' selected active set. The regularization parameter is selected by minimizing a
#' BIC-type criterion over the supplied \code{lambda} grid.
#'
#' @seealso
#' \code{\link{farlr_debias}}, \code{\link{glmnet}}, \code{\link[mirt]{simdata}}
#'
#' @examples
#' \dontrun{
#' fit_emm <- farlr_mml(
#'   X = X,
#'   Y = Y,
#'   parTab = parTab,
#'   method = "FARLR_EMM",
#'   verbose = TRUE
#' )
#'
#' fit_debias <- farlr_mml(
#'   X = X,
#'   Y = Y,
#'   parTab = parTab,
#'   method = "FARLR_Debias",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
Farlr_mml <- function(
    X,
    Y,
    parTab,
    n_sam = 25L,
    method = base::c("FARLR_EMM", "FARLR_Debias"),
    lambda = seq(0.05, 0.5, length.out = 10),
    main = integer(0),
    delta.criteria = 1e-3,
    iter.max = 500L,
    window.size = 200L,
    verbose = FALSE,
    progress = TRUE,
    seed = NULL,
    K_hat = NULL,
    proposal_inflation = 0.2,
    n_sam_final = 60L,
    paran_iterations = 500L,
    paran_centile = 0
) {
  method <- match.arg(method)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required.", call. = FALSE)
  }
  if (is.null(K_hat) && !requireNamespace("paran", quietly = TRUE)) {
    stop(
      "Package 'paran' is required when K_hat is not supplied.",
      call. = FALSE
    )
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.", call. = FALSE)
  }

  suppressWarnings(storage.mode(Y) <- "numeric")
  if (any(!is.finite(Y[!is.na(Y)]))) {
    stop("Observed values in Y must be numeric and finite.", call. = FALSE)
  }

  N <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)

  if (nrow(Y) != N) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }
  if (N < 2L || p < 1L || J < 1L) {
    stop("X and Y have invalid dimensions.", call. = FALSE)
  }

  n_sam <- as.integer(n_sam)
  iter.max <- as.integer(iter.max)
  window.size <- as.integer(window.size)
  if (length(n_sam) != 1L || is.na(n_sam) || n_sam < 1L) {
    stop("n_sam must be a positive integer.", call. = FALSE)
  }
  if (length(iter.max) != 1L || is.na(iter.max) || iter.max < 1L) {
    stop("iter.max must be a positive integer.", call. = FALSE)
  }
  if (length(window.size) != 1L || is.na(window.size) || window.size < 1L) {
    stop("window.size must be a positive integer.", call. = FALSE)
  }

  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L ||
      any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must contain finite nonnegative values.", call. = FALSE)
  }

  main <- unique(as.integer(main))
  if (length(main) > 0L &&
      (anyNA(main) || any(main < 1L | main > p))) {
    stop("main must contain column indices between 1 and ncol(X).", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Names and item-parameter extraction
  # ---------------------------------------------------------------------------
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("item", seq_len(J))
  }
  item_names <- colnames(Y)

  # Standardize X exactly once and use this same scale everywhere below.
  X_original <- X
  X_center <- colMeans(X_original)
  X_scale <- apply(X_original, 2L, stats::sd)

  if (any(!is.finite(X_scale)) || any(X_scale <= 0)) {
    bad_columns <- colnames(X_original)[
      !is.finite(X_scale) | X_scale <= 0
    ]
    stop(
      "X contains zero-variance or invalid columns: ",
      paste(bad_columns, collapse = ", "),
      call. = FALSE
    )
  }

  X <- sweep(X_original, 2L, X_center, "-")
  X <- sweep(X, 2L, X_scale, "/")
  colnames(X) <- colnames(X_original)

  get_parameter <- function(name, default = NULL, required = FALSE) {
    value <- NULL

    if (is.data.frame(parTab) || is.list(parTab)) {
      if (name %in% names(parTab)) {
        value <- parTab[[name]]
      }
    }

    if (is.null(value)) {
      if (required) {
        stop("parTab must contain '", name, "'.", call. = FALSE)
      }
      return(default)
    }

    value
  }

  a_item <- as.numeric(get_parameter("a", required = TRUE))

  # Accept either the mirt intercept d directly or an IRT difficulty b.
  d_item <- get_parameter("d", default = NULL)
  if (is.null(d_item)) {
    b_item <- get_parameter("b", default = NULL)
    if (is.null(b_item)) {
      # d is not used by GPCM items. Binary items are checked after type is built.
      d_item <- rep(NA_real_, J)
    } else {
      b_item <- as.numeric(b_item)
      if (length(b_item) != J) {
        stop("parTab$b must contain one value per item.", call. = FALSE)
      }
      d_item <- -a_item * b_item
    }
  }
  d_item <- as.numeric(d_item)
  c_item <- as.numeric(get_parameter("c", default = rep(0, J)))
  b1_item <- get_parameter("b1", default = rep(NA_real_, J))
  b2_item <- get_parameter("b2", default = rep(NA_real_, J))

  c_item <- rep_len(c_item, J)
  b1_item <- rep_len(as.numeric(b1_item), J)
  b2_item <- rep_len(as.numeric(b2_item), J)

  if (length(a_item) != J || length(d_item) != J) {
    stop(
      "The item parameters a and d/b must each have length ncol(Y).",
      call. = FALSE
    )
  }

  # Construct item type item by item from parTab:
  #   finite b1 and b2 -> GPCM
  #   otherwise c != 0 -> 3PL
  #   otherwise -> 2PL
  has_two_steps <- is.finite(b1_item) & is.finite(b2_item)

  # Missing c is interpreted as c = 0 for non-GPCM items.
  c_item[!has_two_steps & is.na(c_item)] <- 0

  type <- ifelse(
    has_two_steps,
    "gpcm",
    ifelse(c_item != 0, "3pl", "2pl")
  )

  type <- tolower(trimws(as.character(type)))
  if (length(type) != J || anyNA(type)) {
    stop("type must contain one nonmissing value per item.", call. = FALSE)
  }

  allowed_types <- base::c("1pl", "2pl", "3pl", "rasch", "gpcm")
  unknown_types <- setdiff(unique(type), allowed_types)
  if (length(unknown_types) > 0L) {
    stop(
      "Unknown item type(s): ", paste(unknown_types, collapse = ", "),
      call. = FALSE
    )
  }

  idx_gpcm <- type == "gpcm"
  idx_binary <- !idx_gpcm

  if (any(idx_binary) &&
      (any(!is.finite(a_item[idx_binary])) ||
       any(!is.finite(d_item[idx_binary])))) {
    stop("Binary items require finite a and d parameters.", call. = FALSE)
  }
  if (any(idx_binary) &&
      (any(!is.finite(c_item[idx_binary])) ||
       any(c_item[idx_binary] < 0 | c_item[idx_binary] >= 1))) {
    stop("Binary-item c must satisfy 0 <= c < 1.", call. = FALSE)
  }
  if (any(idx_gpcm) &&
      (any(!is.finite(a_item[idx_gpcm])) ||
       any(!is.finite(b1_item[idx_gpcm])) ||
       any(!is.finite(b2_item[idx_gpcm])))) {
    stop("GPCM items require finite a, b1, and b2 parameters.", call. = FALSE)
  }

  if (any(idx_binary)) {
    observed_binary <- Y[, idx_binary, drop = FALSE]
    observed_binary <- observed_binary[!is.na(observed_binary)]
    if (any(!observed_binary %in% base::c(0, 1))) {
      stop("Binary responses must be 0, 1, or NA.", call. = FALSE)
    }
  }
  if (any(idx_gpcm)) {
    observed_gpcm <- Y[, idx_gpcm, drop = FALSE]
    observed_gpcm <- observed_gpcm[!is.na(observed_gpcm)]
    if (any(!observed_gpcm %in% 0:2)) {
      stop("GPCM responses must be 0, 1, 2, or NA.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Factor-adjusted design
  # ---------------------------------------------------------------------------
  if (is.null(K_hat)) {
    parallel_result <- paran::paran(
      X,
      iterations = as.integer(paran_iterations),
      centile = paran_centile,
      quiet = TRUE
    )
    K_hat <- as.integer(parallel_result$Retained)
  } else {
    K_hat <- as.integer(K_hat)
    parallel_result <- NULL
  }

  if (length(K_hat) != 1L || is.na(K_hat) ||
      K_hat < 1L || K_hat > p) {
    stop("K_hat must be an integer between 1 and ncol(X).", call. = FALSE)
  }

  if (!exists("factor.analysis", mode = "function")) {
    stop("Function factor.analysis() was not found.", call. = FALSE)
  }

  fa <- factor.analysis(X, K_hat, method = "ml")
  Wupdate_t <- as.matrix(fa$Gamma)

  if (nrow(Wupdate_t) != p || ncol(Wupdate_t) != K_hat) {
    stop(
      "factor.analysis() returned Gamma with unexpected dimensions.",
      call. = FALSE
    )
  }

  sigma_values <- if (is.matrix(fa$Sigma)) {
    diag(fa$Sigma)
  } else {
    as.numeric(fa$Sigma)
  }
  if (length(sigma_values) != p ||
      any(!is.finite(sigma_values)) || any(sigma_values <= 0)) {
    stop(
      "factor.analysis() must return p positive uniqueness variances in Sigma.",
      call. = FALSE
    )
  }

  Sgm_inv <- diag(1 / sigma_values, nrow = p, ncol = p)
  orthg <- crossprod(Wupdate_t, Sgm_inv %*% Wupdate_t) / p
  orthg <- (orthg + t(orthg)) / 2
  V <- eigen(orthg, symmetric = TRUE)$vectors

  Wupdate <- t(Wupdate_t %*% V)
  score_system <- Wupdate %*% Sgm_inv %*% t(Wupdate)

  Uupdate <- t(
    solve(
      score_system,
      Wupdate %*% Sgm_inv %*% t(X)
    )
  )

  factor_score_names <- paste0("FS", seq_len(K_hat))
  colnames(Uupdate) <- factor_score_names
  Z <- cbind(Uupdate, X)

  # ---------------------------------------------------------------------------
  # IRT proposal distribution for importance sampling
  # ---------------------------------------------------------------------------
  mirt_itemtype <- vapply(
    type,
    function(item_type) {
      switch(
        item_type,
        gpcm = "gpcm",
        `3pl` = "3PL",
        rasch = "Rasch",
        `1pl` = "Rasch",
        "2PL"
      )
    },
    character(1)
  )

  est_mirt <- suppressMessages(
    suppressWarnings(
      mirt::mirt(
        data = Y,
        model = 1,
        itemtype = mirt_itemtype,
        verbose = FALSE
      )
    )
  )

  theta_est_irt <- mirt::fscores(
    est_mirt,
    full.scores = TRUE,
    full.scores.SE = TRUE
  )

  theta_est_irt.mean <- as.numeric(theta_est_irt[, 1L])
  theta_est_irt.se <- as.numeric(theta_est_irt[, 2L])

  if (length(theta_est_irt.mean) != N ||
      length(theta_est_irt.se) != N ||
      any(!is.finite(theta_est_irt.mean)) ||
      any(!is.finite(theta_est_irt.se))) {
    stop("mirt::fscores() returned invalid scores or standard errors.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Method-specific design matrix
  # ---------------------------------------------------------------------------
  if (method == "FARLR_EMM") {
    hatB <- NA
    hatU <- NA
    Fan <- NA
    design_base <- Z
    target_function_name <- "farlr_emm"
  } else {
    hatB <- crossprod(Uupdate, X) / N
    hatU <- X - Uupdate %*% hatB
    colnames(hatU) <- colnames(X)
    Fan <- cbind(Uupdate, hatU)
    colnames(Fan) <- base::c(factor_score_names, colnames(X))
    design_base <- Fan
    target_function_name <- "farlr_debias"
  }

  if (!exists(target_function_name, mode = "function")) {
    stop(
      "Function ", target_function_name, "() was not found. ",
      "Source its R file before calling Farlr_mml().",
      call. = FALSE
    )
  }
  target_function <- get(target_function_name, mode = "function")

  repeated_rows <- rep(seq_len(N), times = n_sam)
  resp_rep <- Y[repeated_rows, , drop = FALSE]
  design_em <- design_base[repeated_rows, , drop = FALSE]

  fit_arguments <- list(
    N = N,
    resp = Y,
    resp.em = resp_rep,
    a = a_item,
    d = d_item,
    c = c_item,
    b1 = b1_item,
    b2 = b2_item,
    type = type,
    K_hat = K_hat,
    p = p,
    lambda_all = lambda,
    delta.criteria = delta.criteria,
    iter.max = iter.max,
    n_sam = n_sam,
    window.size = window.size,
    theta_est_irt.mean = theta_est_irt.mean,
    theta_est_irt.se = theta_est_irt.se,
    Uupdate = Uupdate,
    hatU = hatU,
    Fan = Fan,
    main = main,
    Fan.em = design_em,
    verbose = verbose,
    progress = progress,
    proposal_inflation = proposal_inflation,
    n_sam_final = n_sam_final,
    seed = seed
  )

  # Pass only arguments supported by the selected method. This preserves
  # compatibility with an older farlr_debias() interface that lacks progress,
  # seed, type, or GPCM-specific arguments.
  target_formals <- names(formals(target_function))
  if (!"..." %in% target_formals) {
    fit_arguments <- fit_arguments[
      intersect(names(fit_arguments), target_formals)
    ]
  }

  result <- tryCatch(
    do.call(target_function, fit_arguments),
    error = function(e) {
      stop(
        target_function_name, "() failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # ---------------------------------------------------------------------------
  # Standardize the coefficient output for both fitting methods
  # ---------------------------------------------------------------------------
  coefficient_names <- base::c(
    factor_score_names,
    colnames(X)
  )

  coefficients <- as.numeric(result$coefficients)

  if (length(coefficients) != length(coefficient_names)) {
    stop(
      target_function_name, "() returned ", length(coefficients),
      " coefficients, but ", length(coefficient_names),
      " were expected.",
      call. = FALSE
    )
  }

  names(coefficients) <- coefficient_names
  result$coefficients <- coefficients

  # farlr_debias() historically returned both $coefficients and $coef.
  # Keep one consistently named coefficient component for both methods.
  result$coef <- NULL

  if (is.list(result$all_results)) {
    result$all_results <- lapply(
      result$all_results,
      function(candidate) {
        if (!is.null(candidate$coefficients)) {
          candidate_coefficients <- as.numeric(candidate$coefficients)
          if (length(candidate_coefficients) == length(coefficient_names)) {
            names(candidate_coefficients) <- coefficient_names
            candidate$coefficients <- candidate_coefficients
          }
        }
        candidate$coef <- NULL
        candidate
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Long-format and subject-level outputs
  # ---------------------------------------------------------------------------
  subject <- factor(seq_len(N))
  stuItems <- data.frame(
    subject = rep(subject, times = J),
    key = factor(rep(item_names, each = N), levels = item_names),
    score = as.vector(Y),
    row.names = NULL
  )

  subject_design <- data.frame(
    subject = subject,
    Z,
    check.names = FALSE
  )

  result$stuDat <- subject_design
  result$stuItems <- stuItems
  result$X_original <- X_original
  result$X <- X
  result$X_center <- X_center
  result$X_scale <- X_scale
  result$item_params <- parTab
  result$item_type <- type
  result$factor_scores <- Uupdate
  result$hatB <- hatB
  result$hatU <- hatU
  result$Fan <- Fan
  result$K_hat <- K_hat
  result$parallel_analysis <- parallel_result
  result$method <- method
  result$call <- match.call()

  invisible(result)
}

