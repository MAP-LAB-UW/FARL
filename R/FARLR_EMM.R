#' FARLR EM-M Algorithm for Latent Regression with Regularization
#'
#' Fits a FARLR-style latent regression model using an EM-M algorithm.
#' In the E-step, latent traits are sampled from a normal approximation based on
#' IRT estimates (\code{theta_est_irt.mean}, \code{theta_est_irt.se}). In the M-step,
#' regression coefficients are updated using weighted \code{glmnet} (LASSO) with a
#' debiasing refit via weighted least squares. A sliding window averaging scheme
#' is used to stabilize coefficient updates across iterations. The tuning
#' parameter \code{lambda} is selected by minimizing a BIC-like criterion over
#' \code{lambda_all}.
#'
#' @param N Integer. Sample size (number of persons).
#' @param resp Matrix. Observed item responses of dimension \code{n x J}.
#' @param parTab Data frame. Item parameter table containing at least
#' \code{slope}, \code{difficulty}, and \code{guessin}.
#' @param K_hat Integer. Number of latent factors/components included in \code{Z.em}.
#' @param p Integer. Number of covariates/predictors included in \code{Z.em}.
#' @param lambda_all Numeric vector. Candidate regularization parameters passed to
#' \code{glmnet}.
#' @param delta.criteria Numeric. Convergence tolerance for parameter updates.
#' Default is \code{1e-3}.
#' @param iter.max Integer. Maximum number of EM iterations for each \code{lambda}.
#' Default is \code{500}.
#' @param n_sam Integer. Number of Monte Carlo samples per subject used in the
#' E-step. Default is \code{50}.
#' @param window.size Integer. Sliding window size used to average coefficient
#' updates across iterations. Default is \code{50}.
#' @param theta_est_irt.mean Numeric vector of length \code{n}. IRT-based posterior
#' mean estimates of latent trait \eqn{\theta}.
#' @param theta_est_irt.se Numeric vector of length \code{n}. IRT-based posterior
#' standard error estimates of \eqn{\theta}.
#' @param resp.em Matrix. Replicated/expanded responses used for Monte Carlo
#' computations in the E-step (see \code{q_num_NA}).
#' @param Z.em Matrix. Design matrix used in the regression step, typically of
#' dimension \code{n x (K_hat + p)}.
#' @param main Integer index (or indices). Predictor(s) to be treated as
#' unpenalized via \code{penalty.factor} (set to 0). All other predictors are
#' penalized unless already unpenalized in the first \code{K_hat} columns.
#' @param verbose Logical. If \code{TRUE}, prints progress messages and a progress bar.
#' Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{coefficients}}{Estimated regression coefficients at the selected
#'   \code{lambda} (length \code{K_hat + p}).}
#'   \item{\code{sigma}}{Estimated residual standard deviation.}
#'   \item{\code{LogLik}}{BIC-like objective value corresponding to the selected \code{lambda}
#'   (named \code{LogLik} for compatibility).}
#'   \item{\code{minBIC}}{Index of \code{lambda_all} achieving the minimum BIC criterion.}
#'   \item{\code{Convergence}}{A character flag indicating convergence status.}
#' }
#'
#' @details
#' For each \code{lambda} in \code{lambda_all}, the algorithm iterates until
#' \code{delta} falls below \code{delta.criteria} or \code{iter.max} is reached.
#' The maximum change in coefficient and residual scale estimates is monitored:
#' \eqn{\delta = \max(|\sigma^{(t)}-\sigma^{(t-1)}|, \max_j |\beta_j^{(t)}-\beta_j^{(t-1)}|)}.
#'
#' This function relies on helper functions (not shown here), including
#' \code{q_num_NA()} and \code{add_to_window()}.
#'
#' @seealso \code{\link[glmnet]{glmnet}}, \code{\link[mirt]{simdata}}
#'
#' @examples
#' \dontrun{
#' fit <- farlr_emm(
#'   n = nrow(resp),
#'   resp = resp,
#'   parTab = parTab,
#'   K_hat = K_hat,
#'   p = p,
#'   lambda_all = seq(0.001, 0.1, length.out = 10),
#'   theta_est_irt.mean = theta_mean,
#'   theta_est_irt.se = theta_se,
#'   resp_rep = resp_rep,
#'   Z.em = Z.em,
#'   main = 1,
#'   verbose = TRUE
#' )
#' }
#'
#' @keywords internal
#' @noRd
#'
farlr_emm <- function(
    N,
    resp,
    resp.em,
    a,
    d,
    c,
    b1,
    b2,
    type,
    K_hat,
    p,
    lambda_all,
    delta.criteria = 1e-3,
    iter.max = 500L,
    n_sam = 50L,
    window.size = 50L,
    theta_est_irt.mean,
    theta_est_irt.se,
    Uupdate = NA,
    hatU = NA,
    Fan = NA,
    main = integer(0),
    Fan.em,
    verbose = FALSE,
    progress = TRUE,
    proposal_inflation = 0.2,
    seed = NULL
) {
  # Uupdate, hatU, and Fan are retained for compatibility with existing calls.
  # This function uses Fan.em as the repeated latent-regression design matrix.

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---------------------------------------------------------------------------
  # Input preparation and validation
  # ---------------------------------------------------------------------------
  resp <- as.matrix(resp)
  resp.em <- as.matrix(resp.em)
  Z.em <- as.matrix(Fan.em)

  N <- as.integer(N)
  K_hat <- as.integer(K_hat)
  p <- as.integer(p)
  n_sam <- as.integer(n_sam)
  iter.max <- as.integer(iter.max)
  window.size <- as.integer(window.size)

  if (length(N) != 1L || is.na(N) || N < 1L) {
    stop("N must be a positive integer.", call. = FALSE)
  }
  if (length(K_hat) != 1L || is.na(K_hat) || K_hat < 0L) {
    stop("K_hat must be a nonnegative integer.", call. = FALSE)
  }
  if (length(p) != 1L || is.na(p) || p < 1L) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  if (length(n_sam) != 1L || is.na(n_sam) || n_sam < 1L) {
    stop("n_sam must be a positive integer.", call. = FALSE)
  }
  if (length(iter.max) != 1L || is.na(iter.max) || iter.max < 1L) {
    stop("iter.max must be a positive integer.", call. = FALSE)
  }
  if (length(window.size) != 1L || is.na(window.size) || window.size < 1L) {
    stop("window.size must be a positive integer.", call. = FALSE)
  }
  if (length(delta.criteria) != 1L ||
      !is.finite(delta.criteria) || delta.criteria <= 0) {
    stop("delta.criteria must be a positive finite number.", call. = FALSE)
  }
  if (length(proposal_inflation) != 1L ||
      !is.finite(proposal_inflation) || proposal_inflation < 0) {
    stop("proposal_inflation must be a finite nonnegative number.", call. = FALSE)
  }

  J <- ncol(resp)
  Q <- N * n_sam
  n_coef <- K_hat + p

  if (nrow(resp) != N) {
    stop("nrow(resp) must equal N.", call. = FALSE)
  }
  if (nrow(resp.em) != Q || ncol(resp.em) != J) {
    stop(
      "resp.em must have N * n_sam rows and ncol(resp) columns.",
      call. = FALSE
    )
  }
  if (nrow(Z.em) != Q || ncol(Z.em) != n_coef) {
    stop(
      "Fan.em must have N * n_sam rows and K_hat + p columns.",
      call. = FALSE
    )
  }
  if (!is.numeric(Z.em) || any(!is.finite(Z.em))) {
    stop("Fan.em must be a finite numeric matrix.", call. = FALSE)
  }

  # The following ordering is required for both resp.em and Fan.em:
  # subjects 1,...,N for draw 1, then subjects 1,...,N for draw 2, etc.
  expected_resp_em <- resp[rep(seq_len(N), times = n_sam), , drop = FALSE]
  same_response <- isTRUE(all.equal(
    unname(resp.em),
    unname(expected_resp_em),
    check.attributes = FALSE
  ))
  if (!same_response) {
    stop(
      "resp.em has the wrong row order. Construct it with ",
      "resp[rep(seq_len(N), times = n_sam), , drop = FALSE].",
      call. = FALSE
    )
  }

  if (length(theta_est_irt.mean) != N ||
      any(!is.finite(theta_est_irt.mean))) {
    stop(
      "theta_est_irt.mean must contain N finite values.",
      call. = FALSE
    )
  }
  if (length(theta_est_irt.se) != N ||
      any(!is.finite(theta_est_irt.se))) {
    stop(
      "theta_est_irt.se must contain N finite values.",
      call. = FALSE
    )
  }

  proposal_sd <- theta_est_irt.se + proposal_inflation
  if (any(proposal_sd <= 0)) {
    stop(
      "theta_est_irt.se + proposal_inflation must be positive.",
      call. = FALSE
    )
  }

  lambda_all <- as.numeric(lambda_all)
  if (length(lambda_all) == 0L ||
      any(!is.finite(lambda_all)) || any(lambda_all < 0)) {
    stop(
      "lambda_all must contain finite nonnegative values.",
      call. = FALSE
    )
  }

  main <- as.integer(main)
  if (length(main) > 0L &&
      (anyNA(main) || any(main < 1L | main > p))) {
    stop("main must contain indices between 1 and p.", call. = FALSE)
  }
  main <- unique(main)

  type <- tolower(trimws(as.character(type)))
  if (length(type) != J || anyNA(type)) {
    stop("type must contain one nonmissing value per item.", call. = FALSE)
  }

  accepted_types <- base::c("1pl", "2pl", "3pl", "rasch", "gpcm")
  unknown_types <- setdiff(unique(type), accepted_types)
  if (length(unknown_types) > 0L) {
    stop(
      "Unknown item type(s): ",
      paste(unknown_types, collapse = ", "),
      call. = FALSE
    )
  }

  idx_gpcm <- type == "gpcm"
  idx_binary <- !idx_gpcm

  a <- as.numeric(a)
  d <- as.numeric(d)
  guessing <- rep_len(as.numeric(c), J)

  if (length(a) != J || any(!is.finite(a))) {
    stop("a must contain one finite value per item.", call. = FALSE)
  }
  if (length(d) != J) {
    stop("d must contain one value per item.", call. = FALSE)
  }
  if (any(idx_binary) && any(!is.finite(d[idx_binary]))) {
    stop("Binary items require finite d parameters.", call. = FALSE)
  }
  if (any(idx_binary) &&
      (any(!is.finite(guessing[idx_binary])) ||
       any(guessing[idx_binary] < 0 | guessing[idx_binary] >= 1))) {
    stop(
      "Binary-item c parameters must satisfy 0 <= c < 1.",
      call. = FALSE
    )
  }

  if (is.null(b1)) b1 <- rep(NA_real_, J)
  if (is.null(b2)) b2 <- rep(NA_real_, J)
  b1 <- as.numeric(b1)
  b2 <- as.numeric(b2)

  if (any(idx_gpcm)) {
    if (length(b1) != J || length(b2) != J) {
      stop(
        "When GPCM items are present, b1 and b2 must have length J.",
        call. = FALSE
      )
    }
    if (any(!is.finite(b1[idx_gpcm])) ||
        any(!is.finite(b2[idx_gpcm]))) {
      stop("GPCM items require finite b1 and b2 values.", call. = FALSE)
    }
  }

  # Validate observed response categories once rather than in every iteration.
  if (any(idx_binary)) {
    observed_binary <- resp.em[, idx_binary, drop = FALSE]
    observed_binary <- observed_binary[!is.na(observed_binary)]
    if (any(!observed_binary %in% base::c(0, 1))) {
      stop("Binary responses must be 0, 1, or NA.", call. = FALSE)
    }
  }
  if (any(idx_gpcm)) {
    observed_gpcm <- resp.em[, idx_gpcm, drop = FALSE]
    observed_gpcm <- observed_gpcm[!is.na(observed_gpcm)]
    if (any(!observed_gpcm %in% 0:2)) {
      stop("GPCM responses must be 0, 1, 2, or NA.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Mixed binary/GPCM log numerator
  # ---------------------------------------------------------------------------
  log_q_num <- function(theta, beta, sigma) {
    Q_local <- length(theta)
    log_item_likelihood <- numeric(Q_local)

    if (any(idx_binary)) {
      a_binary <- a[idx_binary]
      d_binary <- d[idx_binary]
      c_binary <- guessing[idx_binary]

      eta_binary <-
        outer(theta, a_binary, "*") +
        matrix(rep(d_binary, each = Q_local), nrow = Q_local)

      c_matrix <- matrix(
        rep(c_binary, each = Q_local),
        nrow = Q_local
      )
      probability_binary <-
        c_matrix + (1 - c_matrix) * plogis(eta_binary)
      probability_binary <- pmin(
        pmax(probability_binary, 1e-12),
        1 - 1e-12
      )

      y_binary <- resp.em[, idx_binary, drop = FALSE]
      ll_binary <-
        y_binary * log(probability_binary) +
        (1 - y_binary) * log1p(-probability_binary)
      ll_binary[is.na(y_binary)] <- 0

      log_item_likelihood <-
        log_item_likelihood + rowSums(ll_binary)
    }

    if (any(idx_gpcm)) {
      a_gpcm <- a[idx_gpcm]
      b1_gpcm <- b1[idx_gpcm]
      b2_gpcm <- b2[idx_gpcm]
      y_gpcm <- resp.em[, idx_gpcm, drop = FALSE]

      theta_a <- outer(theta, a_gpcm, "*")

      # P1/P0 = exp{a(theta-b1)}
      eta1 <-
        theta_a -
        matrix(rep(a_gpcm * b1_gpcm, each = Q_local), nrow = Q_local)

      # P2/P0 = exp{a(theta-b1) + a(theta-b2)}
      eta2 <-
        2 * theta_a -
        matrix(
          rep(a_gpcm * (b1_gpcm + b2_gpcm), each = Q_local),
          nrow = Q_local
        )

      zero_matrix <- matrix(0, nrow = nrow(eta1), ncol = ncol(eta1))
      max_eta <- pmax(zero_matrix, eta1, eta2)
      log_denominator <-
        max_eta +
        log(
          exp(-max_eta) +
            exp(eta1 - max_eta) +
            exp(eta2 - max_eta)
        )

      log_P0 <- -log_denominator
      log_P1 <- eta1 - log_denominator
      log_P2 <- eta2 - log_denominator

      ll_gpcm <- matrix(0, nrow = Q_local, ncol = sum(idx_gpcm))
      obs0 <- !is.na(y_gpcm) & y_gpcm == 0
      obs1 <- !is.na(y_gpcm) & y_gpcm == 1
      obs2 <- !is.na(y_gpcm) & y_gpcm == 2
      ll_gpcm[obs0] <- log_P0[obs0]
      ll_gpcm[obs1] <- log_P1[obs1]
      ll_gpcm[obs2] <- log_P2[obs2]

      log_item_likelihood <-
        log_item_likelihood + rowSums(ll_gpcm)
    }

    mu <- drop(Z.em %*% beta)
    log_regression_density <- dnorm(
      theta,
      mean = mu,
      sd = sigma,
      log = TRUE
    )

    log_regression_density + log_item_likelihood
  }

  # ---------------------------------------------------------------------------
  # Progress bar. It is created once and closed only after every lambda finishes.
  # ---------------------------------------------------------------------------
  progress_enabled <- isTRUE(progress)
  total_steps <- length(lambda_all) * iter.max
  progress_bar <- NULL

  if (progress_enabled) {
    progress_bar <- utils::txtProgressBar(
      min = 0,
      max = total_steps,
      initial = 0,
      style = 3
    )
  }

  close_progress_bar <- function() {
    if (inherits(progress_bar, "txtProgressBar")) {
      close(progress_bar)
      progress_bar <<- NULL
    }
    invisible(NULL)
  }

  update_progress_bar <- function(value) {
    if (inherits(progress_bar, "txtProgressBar")) {
      utils::setTxtProgressBar(
        progress_bar,
        min(as.numeric(value), total_steps)
      )
    }
    invisible(NULL)
  }

  on.exit(close_progress_bar(), add = TRUE)

  # ---------------------------------------------------------------------------
  # Penalization setup
  # ---------------------------------------------------------------------------
  penalty_factors_x <- rep(1, p)
  if (length(main) > 0L) {
    penalty_factors_x[main] <- 0
  }
  penalty_factors <- base::c(rep(0, K_hat), penalty_factors_x)

  results <- vector("list", length(lambda_all))

  # ---------------------------------------------------------------------------
  # Monte Carlo EM for each lambda
  # ---------------------------------------------------------------------------
  for (ll in seq_along(lambda_all)) {
    lambda <- lambda_all[ll]
    beta_old <- rep(0, n_coef)
    beta_window <- list()
    sigma_old <- 1
    delta <- Inf
    iter <- 0L

    covariance_current <- matrix(0, nrow = n_coef, ncol = n_coef)
    bic_current <- Inf
    logLik_current <- -Inf
    ess_summary_current <- base::c(
      minimum = NA_real_,
      median = NA_real_,
      mean = NA_real_
    )

    if (isTRUE(verbose) && !progress_enabled) {
      message(
        "Starting lambda ", ll, " of ", length(lambda_all),
        " | lambda = ", signif(lambda, 6)
      )
    }

    while (delta > delta.criteria && iter < iter.max) {
      iter <- iter + 1L

      # E-step: subjects 1,...,N are repeated for each Monte Carlo draw.
      theta_sample <- rnorm(
        Q,
        mean = theta_est_irt.mean,
        sd = proposal_sd
      )

      log_q <- log_q_num(
        theta = theta_sample,
        beta = beta_old,
        sigma = sigma_old
      )
      log_h <- dnorm(
        theta_sample,
        mean = theta_est_irt.mean,
        sd = proposal_sd,
        log = TRUE
      )

      # Normalize importance weights in log space, separately for each subject.
      log_ratio_matrix <- matrix(
        log_q - log_h,
        nrow = N,
        ncol = n_sam
      )
      row_max <- apply(log_ratio_matrix, 1L, max)

      if (any(!is.finite(row_max))) {
        stop(
          "Non-finite log importance ratios at lambda index ", ll,
          ", iteration ", iter, ".",
          call. = FALSE
        )
      }

      scaled_ratio <- exp(sweep(log_ratio_matrix, 1L, row_max, "-"))
      row_mean <- rowMeans(scaled_ratio)
      if (any(!is.finite(row_mean)) || any(row_mean <= 0)) {
        stop(
          "Invalid importance weights at lambda index ", ll,
          ", iteration ", iter, ".",
          call. = FALSE
        )
      }

      weight_matrix <- sweep(scaled_ratio, 1L, row_mean, "/")
      w_ik <- as.vector(weight_matrix)

      # Effective sample size, using weights normalized within each subject.
      normalized_weight_matrix <-
        weight_matrix / rowSums(weight_matrix)
      ess <- 1 / rowSums(normalized_weight_matrix^2)
      ess_summary_current <- base::c(
        minimum = min(ess),
        median = stats::median(ess),
        mean = mean(ess)
      )

      # Penalized weighted regression.
      glmnet_fit <- glmnet::glmnet(
        x = Z.em,
        y = theta_sample,
        family = "gaussian",
        weights = w_ik,
        penalty.factor = penalty_factors,
        intercept = FALSE,
        standardize = FALSE,
        lambda = lambda
      )

      coef_lasso <- as.numeric(
        stats::coef(glmnet_fit, s = lambda)
      )[-1L]
      active <- which(coef_lasso != 0)

      # Debias on the active set using weighted least squares.
      beta_debias <- rep(0, n_coef)
      covariance_current <- matrix(0, nrow = n_coef, ncol = n_coef)

      if (length(active) > 0L) {
        Z_active <- Z.em[, active, drop = FALSE]
        debias_fit <- stats::lm(
          theta_sample ~ 0 + Z_active,
          weights = w_ik
        )

        active_coef <- stats::coef(debias_fit)
        active_coef[is.na(active_coef)] <- 0
        beta_debias[active] <- active_coef

        active_covariance <- tryCatch(
          stats::vcov(debias_fit),
          error = function(e) {
            matrix(
              NA_real_,
              nrow = length(active),
              ncol = length(active)
            )
          }
        )
        covariance_current[active, active] <- active_covariance
      }

      # Moving-window average of debiased coefficients.
      beta_window[[length(beta_window) + 1L]] <- beta_debias
      if (length(beta_window) > window.size) {
        beta_window <- tail(beta_window, window.size)
      }
      beta_mean <- rowMeans(do.call(base::cbind, beta_window))

      # Sigma update. Use the active model degrees of freedom rather than all p.
      residuals <- theta_sample - drop(Z.em %*% beta_mean)
      weighted_ss <- sum(w_ik * residuals^2)
      model_df <- sum(beta_debias != 0)
      sigma_denominator <- sum(w_ik) - model_df

      if (!is.finite(sigma_denominator) || sigma_denominator <= 0) {
        stop("Invalid denominator in sigma update.", call. = FALSE)
      }

      sigma_new <- sqrt(weighted_ss / sigma_denominator)
      if (!is.finite(sigma_new) || sigma_new <= 0) {
        stop(
          "Invalid sigma estimate at lambda index ", ll,
          ", iteration ", iter, ".",
          call. = FALSE
        )
      }

      # Preserve the original BIC definition based on the current sparse model.
      current_residuals <- theta_sample - drop(Z.em %*% beta_debias)
      logLik_current <-
        -N / 2 * log(2 * pi * sigma_new^2) -
        sum(w_ik * current_residuals^2) /
        (2 * n_sam * sigma_new^2)
      bic_current <-
        -2 * logLik_current + model_df * log(N)

      delta <- max(
        abs(sigma_old - sigma_new),
        max(abs(beta_old - beta_mean))
      )
      beta_old <- beta_mean
      sigma_old <- sigma_new

      update_progress_bar((ll - 1L) * iter.max + iter)

      if (isTRUE(verbose) && !progress_enabled) {
        message(
          "  iter = ", iter,
          " | delta = ", signif(delta, 4),
          " | sigma = ", signif(sigma_old, 4),
          " | BIC = ", signif(bic_current, 6),
          " | median ESS = ", signif(ess_summary_current[["median"]], 4)
        )
      }
    }

    converged <- is.finite(delta) && delta <= delta.criteria
    results[[ll]] <- list(
      coefficients = beta_old,
      sigma = sigma_old,
      BIC = bic_current,
      LogLik = logLik_current,
      covariance = covariance_current,
      lambda = lambda,
      lambda_index = ll,
      iterations = iter,
      delta = delta,
      converged = converged,
      effective_sample_size = ess_summary_current
    )

    # Mark the whole lambda block complete if it converged before iter.max.
    update_progress_bar(ll * iter.max)

    if (isTRUE(verbose) && !progress_enabled) {
      message(
        "Finished lambda ", ll,
        " | iterations = ", iter,
        " | final delta = ", signif(delta, 4),
        " | BIC = ", signif(bic_current, 6),
        " | ", if (converged) "converged" else "did not converge"
      )
    }
  }

  update_progress_bar(total_steps)
  close_progress_bar()

  # ---------------------------------------------------------------------------
  # Select lambda by minimum finite BIC
  # ---------------------------------------------------------------------------
  BIC_values <- vapply(results, function(x) x$BIC, numeric(1))
  if (all(!is.finite(BIC_values))) {
    stop("No finite BIC value was obtained.", call. = FALSE)
  }
  BIC_for_selection <- BIC_values
  BIC_for_selection[!is.finite(BIC_for_selection)] <- Inf
  minBIC <- which.min(BIC_for_selection)
  best <- results[[minBIC]]

  if (isTRUE(verbose)) {
    message(
      "Selected lambda index: ", minBIC,
      " | lambda = ", signif(best$lambda, 6),
      " | BIC = ", signif(best$BIC, 6)
    )
  }

  list(
    coefficients = best$coefficients,
    sigma = best$sigma,
    BIC = best$BIC,
    LogLik = best$LogLik,
    covariance = best$covariance,
    minBIC = minBIC,
    lambda = best$lambda,
    iterations = best$iterations,
    delta = best$delta,
    effective_sample_size = best$effective_sample_size,
    Convergence = if (best$converged) "Converged" else "Did not converge",
    all_results = results
  )
}


