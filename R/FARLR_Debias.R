library(torch)
#' FARLR Debiased Estimation for Regularized Latent Regression
#'
#' Fits a Factor-Augmented Regularized Latent Regression (FARLR) model with a
#' post-selection debiasing step. The method uses Monte Carlo samples from a
#' normal approximation to the IRT posterior distribution of the latent trait,
#' based on \code{theta_est_irt.mean} and \code{theta_est_irt.se}. For each
#' candidate value in \code{lambda_all}, regression coefficients are estimated
#' by weighted LASSO using \code{glmnet}. The selected coefficients are then
#' debiased using a weighted correction step. Coefficient updates are smoothed
#' across iterations using a sliding-window average, and the final tuning
#' parameter is chosen by minimizing a BIC-type criterion.
#'
#' This implementation avoids forming large projection and diagonal weight
#' matrices directly. Instead, it uses equivalent matrix products, which improves
#' computational efficiency while preserving the intended calculations.
#'
#' @param N Integer. Number of individuals.
#' @param resp Matrix. Observed item response matrix of dimension \code{N x J}.
#' @param resp.em Matrix. Expanded response matrix used in the Monte Carlo
#'   integration step.
#' @param a Numeric vector. Item slope parameters.
#' @param d Numeric vector. Item difficulty or intercept parameters.
#' @param c Numeric vector. Item guessing parameters.
#' @param b1 Numeric vector. Additional item parameter used by \code{q_num_NA()}.
#' @param b2 Numeric vector. Additional item parameter used by \code{q_num_NA()}.
#' @param type Object specifying the item model type used by \code{q_num_NA()}.
#' @param K_hat Integer. Number of estimated latent factor components.
#' @param p Integer. Number of observed covariates.
#' @param lambda_all Numeric vector. Candidate regularization parameters passed
#'   to \code{glmnet}.
#' @param delta.criteria Numeric. Convergence tolerance for the iterative updates.
#'   Defaults to \code{1e-3}.
#' @param iter.max Integer. Maximum number of iterations for each value of
#'   \code{lambda}. Defaults to \code{500}.
#' @param n_sam Integer. Number of Monte Carlo samples per individual used during
#'   tuning. Defaults to \code{5}.
#' @param window.size Integer. Window size used for sliding-window averaging of
#'   coefficient updates. Defaults to \code{50}.
#' @param theta_est_irt.mean Numeric vector of length \code{n}. IRT posterior mean
#'   estimates of the latent trait.
#' @param theta_est_irt.se Numeric vector of length \code{n}. IRT posterior
#'   standard error estimates of the latent trait.
#' @param Uupdate Matrix. Factor/design matrix used in the projection and
#'   debiasing updates.
#' @param hatU Matrix. Design matrix used in the penalized regression step.
#' @param Fan Matrix. Final regression design matrix used when recomputing the
#'   residual scale after selecting \code{lambda}.
#' @param main Integer vector. Indices of covariates that are left unpenalized in
#'   the \code{glmnet} fit through \code{penalty.factor = 0}.
#' @param Fan.em Matrix. Expanded regression design matrix used during the main
#'   iterative estimation step.
#' @param verbose boolean. Output the intermediate steps or not.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{coef}}{Estimated FARLR debiased regression coefficients for the
#'   selected value of \code{lambda}.}
#'   \item{\code{sigma}}{Estimated residual standard deviation after selecting
#'   \code{lambda}.}
#'   \item{\code{minBIC}}{Index of the value in \code{lambda_all} that minimizes
#'   the BIC-type criterion.}
#' }
#'
#' @details
#' For each value of \code{lambda_all}, the algorithm repeatedly samples latent
#' trait values, computes importance weights, estimates penalized regression
#' coefficients, applies a debiasing correction, and updates the residual
#' standard deviation. Iteration stops when the maximum change in the coefficient
#' vector or residual standard deviation is below \code{delta.criteria}, or when
#' \code{iter.max} is reached.
#'
#' The tuning parameter is selected by minimizing a BIC-type objective. After
#' selection, the residual standard deviation is recomputed using a larger Monte
#' Carlo sample size.
#'
#' This function depends on auxiliary routines, including \code{q_num_NA()} for
#' Monte Carlo integration and \code{add_to_window()} for sliding-window
#' averaging.
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
#' @examples
#' \dontrun{
#' fit <- farlr_debias(
#'   N = nrow(resp),
#'   resp = resp,
#'   resp.em = resp.em,
#'   a = a,
#'   d = d,
#'   c = c,
#'   b1 = b1,
#'   b2 = b2,
#'   type = type,
#'   K_hat = K_hat,
#'   p = p,
#'   lambda_all = seq(0.001, 0.1, length.out = 10),
#'   theta_est_irt.mean = theta_mean,
#'   theta_est_irt.se = theta_se,
#'   Uupdate = Uupdate,
#'   hatU = hatU,
#'   Fan = Fan,
#'   main = 1,
#'   Fan.em = Fan.em,
#'   verbose = TRUE
#' )
#' }
#'
#' @keywords internal
#' @noRd
#'
farlr_debias <- function(
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
    Uupdate,
    hatU,
    Fan,
    main = integer(0),
    Fan.em,
    verbose = FALSE,
    progress = TRUE,
    proposal_inflation = 0.2,
    n_sam_final = 60L,
    seed = NULL
) {
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
  Uupdate <- as.matrix(Uupdate)
  hatU <- as.matrix(hatU)
  Fan <- as.matrix(Fan)
  Z.em <- as.matrix(Fan.em)

  N <- as.integer(N)
  K_hat <- as.integer(K_hat)
  p <- as.integer(p)
  n_sam <- as.integer(n_sam)
  n_sam_final <- as.integer(n_sam_final)
  iter.max <- as.integer(iter.max)
  window.size <- as.integer(window.size)

  if (length(N) != 1L || is.na(N) || N < 1L) {
    stop("N must be a positive integer.", call. = FALSE)
  }
  if (length(K_hat) != 1L || is.na(K_hat) || K_hat < 1L) {
    stop("K_hat must be a positive integer.", call. = FALSE)
  }
  if (length(p) != 1L || is.na(p) || p < 1L) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  if (length(n_sam) != 1L || is.na(n_sam) || n_sam < 1L) {
    stop("n_sam must be a positive integer.", call. = FALSE)
  }
  if (length(n_sam_final) != 1L || is.na(n_sam_final) || n_sam_final < 1L) {
    stop("n_sam_final must be a positive integer.", call. = FALSE)
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
    stop("proposal_inflation must be finite and nonnegative.", call. = FALSE)
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
  if (!identical(dim(Uupdate), base::c(N, K_hat))) {
    stop("Uupdate must have dimensions N by K_hat.", call. = FALSE)
  }
  if (!identical(dim(hatU), base::c(N, p))) {
    stop("hatU must have dimensions N by p.", call. = FALSE)
  }
  if (!identical(dim(Fan), base::c(N, n_coef))) {
    stop("Fan must have dimensions N by K_hat + p.", call. = FALSE)
  }
  if (!identical(dim(Z.em), base::c(Q, n_coef))) {
    stop("Fan.em must have dimensions N * n_sam by K_hat + p.", call. = FALSE)
  }
  if (any(!is.finite(Uupdate)) || any(!is.finite(hatU)) ||
      any(!is.finite(Fan)) || any(!is.finite(Z.em))) {
    stop("Uupdate, hatU, Fan, and Fan.em must be finite.", call. = FALSE)
  }

  repeated_rows <- rep(seq_len(N), times = n_sam)
  expected_resp_em <- resp[repeated_rows, , drop = FALSE]
  expected_Fan_em <- Fan[repeated_rows, , drop = FALSE]

  if (!isTRUE(all.equal(
    unname(resp.em), unname(expected_resp_em), check.attributes = FALSE
  ))) {
    stop(
      "resp.em has the wrong row order. Repeat subjects 1,...,N for each draw.",
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(
    unname(Z.em), unname(expected_Fan_em),
    tolerance = 1e-10, check.attributes = FALSE
  ))) {
    stop(
      "Fan.em must be Fan[rep(seq_len(N), times = n_sam), ].",
      call. = FALSE
    )
  }

  if (length(theta_est_irt.mean) != N ||
      any(!is.finite(theta_est_irt.mean))) {
    stop("theta_est_irt.mean must contain N finite values.", call. = FALSE)
  }
  if (length(theta_est_irt.se) != N ||
      any(!is.finite(theta_est_irt.se))) {
    stop("theta_est_irt.se must contain N finite values.", call. = FALSE)
  }
  proposal_sd <- theta_est_irt.se + proposal_inflation
  if (any(proposal_sd <= 0)) {
    stop("The proposal standard deviations must be positive.", call. = FALSE)
  }

  lambda_all <- as.numeric(lambda_all)
  if (length(lambda_all) == 0L ||
      any(!is.finite(lambda_all)) || any(lambda_all < 0)) {
    stop("lambda_all must contain finite nonnegative values.", call. = FALSE)
  }

  main <- unique(as.integer(main))
  if (length(main) > 0L &&
      (anyNA(main) || any(main < 1L | main > p))) {
    stop("main must contain indices between 1 and p.", call. = FALSE)
  }

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
    stop("Binary-item c must satisfy 0 <= c < 1.", call. = FALSE)
  }

  if (is.null(b1)) b1 <- rep(NA_real_, J)
  if (is.null(b2)) b2 <- rep(NA_real_, J)
  b1 <- as.numeric(b1)
  b2 <- as.numeric(b2)
  if (any(idx_gpcm)) {
    if (length(b1) != J || length(b2) != J ||
        any(!is.finite(b1[idx_gpcm])) ||
        any(!is.finite(b2[idx_gpcm]))) {
      stop("GPCM items require finite b1 and b2 values.", call. = FALSE)
    }
  }

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
  # Stable mixed binary/GPCM log numerator
  # ---------------------------------------------------------------------------
  log_q_num <- function(theta, resp_q, Z_q, beta, sigma) {
    Q_local <- length(theta)
    log_item_likelihood <- numeric(Q_local)

    if (any(idx_binary)) {
      a_binary <- a[idx_binary]
      d_binary <- d[idx_binary]
      c_binary <- guessing[idx_binary]

      eta <-
        outer(theta, a_binary, "*") +
        matrix(rep(d_binary, each = Q_local), nrow = Q_local)
      c_matrix <- matrix(rep(c_binary, each = Q_local), nrow = Q_local)
      probability <- c_matrix + (1 - c_matrix) * plogis(eta)
      probability <- pmin(pmax(probability, 1e-12), 1 - 1e-12)

      y <- resp_q[, idx_binary, drop = FALSE]
      ll <- y * log(probability) + (1 - y) * log1p(-probability)
      ll[is.na(y)] <- 0
      log_item_likelihood <- log_item_likelihood + rowSums(ll)
    }

    if (any(idx_gpcm)) {
      a_gpcm <- a[idx_gpcm]
      b1_gpcm <- b1[idx_gpcm]
      b2_gpcm <- b2[idx_gpcm]
      y <- resp_q[, idx_gpcm, drop = FALSE]

      theta_a <- outer(theta, a_gpcm, "*")
      eta1 <-
        theta_a -
        matrix(rep(a_gpcm * b1_gpcm, each = Q_local), nrow = Q_local)
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
      ll <- matrix(0, nrow = Q_local, ncol = sum(idx_gpcm))
      obs0 <- !is.na(y) & y == 0
      obs1 <- !is.na(y) & y == 1
      obs2 <- !is.na(y) & y == 2
      ll[obs0] <- log_P0[obs0]
      ll[obs1] <- log_P1[obs1]
      ll[obs2] <- log_P2[obs2]
      log_item_likelihood <- log_item_likelihood + rowSums(ll)
    }

    mu <- drop(Z_q %*% beta)
    log_regression_density <- dnorm(
      theta,
      mean = mu,
      sd = sigma,
      log = TRUE
    )

    log_regression_density + log_item_likelihood
  }

  importance_weights <- function(
    theta,
    resp_q,
    Z_q,
    beta,
    sigma,
    n_draws
  ) {
    log_q <- log_q_num(theta, resp_q, Z_q, beta, sigma)
    log_h <- dnorm(
      theta,
      mean = theta_est_irt.mean,
      sd = proposal_sd,
      log = TRUE
    )

    log_ratio_matrix <- matrix(
      log_q - log_h,
      nrow = N,
      ncol = n_draws
    )
    row_max <- apply(log_ratio_matrix, 1L, max)
    if (any(!is.finite(row_max))) {
      stop("Non-finite log importance ratios.", call. = FALSE)
    }

    scaled_ratio <- exp(sweep(log_ratio_matrix, 1L, row_max, "-"))
    row_mean <- rowMeans(scaled_ratio)
    if (any(!is.finite(row_mean)) || any(row_mean <= 0)) {
      stop("Invalid importance weights.", call. = FALSE)
    }

    weight_matrix <- sweep(scaled_ratio, 1L, row_mean, "/")
    normalized <- weight_matrix / rowSums(weight_matrix)
    ess <- 1 / rowSums(normalized^2)

    list(
      weights = as.vector(weight_matrix),
      ess = ess
    )
  }

  # ---------------------------------------------------------------------------
  # Progress bar: created once, closed after tuning and final sigma update
  # ---------------------------------------------------------------------------
  progress_enabled <- isTRUE(progress)
  total_steps <- length(lambda_all) * iter.max + 1L
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
  # Fixed expanded matrices and penalty factors
  # ---------------------------------------------------------------------------
  Uupdate.em <- Uupdate[repeated_rows, , drop = FALSE]
  hatU.em <- hatU[repeated_rows, , drop = FALSE]
  inv_nnsam <- 1 / Q

  penalty_factors <- rep(1, p)
  if (length(main) > 0L) {
    penalty_factors[main] <- 0
  }

  results <- vector("list", length(lambda_all))

  # ---------------------------------------------------------------------------
  # Monte Carlo EM for every lambda
  # ---------------------------------------------------------------------------
  for (ll in seq_along(lambda_all)) {
    lambda <- lambda_all[ll]
    beta_old <- rep(0, n_coef)
    beta_window <- list()
    sigma_old <- 1
    delta <- Inf
    iter <- 0L
    bic_current <- Inf
    logLik_current <- -Inf
    ess_summary_current <- base::c(
      minimum = NA_real_, median = NA_real_, mean = NA_real_
    )

    if (isTRUE(verbose) && !progress_enabled) {
      message(
        "Starting lambda ", ll, " of ", length(lambda_all),
        " | lambda = ", signif(lambda, 6)
      )
    }

    while (delta > delta.criteria && iter < iter.max) {
      iter <- iter + 1L

      theta_sample <- rnorm(
        Q,
        mean = theta_est_irt.mean,
        sd = proposal_sd
      )

      weight_result <- importance_weights(
        theta = theta_sample,
        resp_q = resp.em,
        Z_q = Z.em,
        beta = beta_old,
        sigma = sigma_old,
        n_draws = n_sam
      )
      w_ik <- weight_result$weights
      ess <- weight_result$ess
      ess_summary_current <- base::c(
        minimum = min(ess),
        median = stats::median(ess),
        mean = mean(ess)
      )

      # Original projection, computed without forming a Q by Q matrix.
      projection <- drop(
        Uupdate.em %*% (crossprod(Uupdate.em, theta_sample) * inv_nnsam)
      )
      theta_1 <- theta_sample - projection

      lasso_fit <- glmnet::glmnet(
        x = hatU.em,
        y = theta_1,
        weights = w_ik,
        penalty.factor = penalty_factors,
        intercept = FALSE,
        # Preserve the original farlr_debias() behavior. glmnet's original
        # default is standardize = TRUE; changing this changes the selected
        # variables and the debiased coefficients.
        standardize = TRUE,
        family = "gaussian",
        lambda = lambda
      )
      coef_lasso <- as.numeric(
        stats::coef(lasso_fit, s = lambda)
      )[-1L]

      # Diagonal debiasing correction from the original implementation.
      T_diagonal <- colSums(w_ik * hatU.em^2) * inv_nnsam
      if (any(!is.finite(T_diagonal)) || any(T_diagonal <= 0)) {
        stop(
          "Invalid diagonal debiasing matrix at lambda index ", ll,
          ", iteration ", iter, ".",
          call. = FALSE
        )
      }

      lasso_residuals <- theta_1 - drop(hatU.em %*% coef_lasso)
      correction <-
        colSums(hatU.em * (w_ik * lasso_residuals)) *
        inv_nnsam / T_diagonal
      coef_debias <- coef_lasso + correction

      # Preserve sparsity selected by the penalized fit.
      coef_debias[coef_lasso == 0] <- 0
      coef_debias[!is.finite(coef_debias)] <- 0

      # Weighted factor coefficients.
      Uw <- Uupdate.em * w_ik
      lhs <- crossprod(Uupdate.em, Uw)
      rhs <- crossprod(Uupdate.em, theta_sample * w_ik)
      phi_hat <- tryCatch(
        as.numeric(solve(lhs, rhs)),
        error = function(e) {
          stop(
            "Factor coefficient solve failed at lambda index ", ll,
            ", iteration ", iter, ": ", conditionMessage(e),
            call. = FALSE
          )
        }
      )

      beta_current <- base::c(phi_hat, coef_debias)
      beta_window[[length(beta_window) + 1L]] <- beta_current
      if (length(beta_window) > window.size) {
        beta_window <- tail(beta_window, window.size)
      }
      beta_mean <- rowMeans(do.call(base::cbind, beta_window))

      residuals <- theta_sample - drop(Z.em %*% beta_mean)
      weighted_ss <- sum(w_ik * residuals^2)
      sigma_denominator <- sum(w_ik) - p
      if (!is.finite(sigma_denominator) || sigma_denominator <= 0) {
        stop("Invalid denominator in sigma update.", call. = FALSE)
      }
      sigma_new <- sqrt(weighted_ss / sigma_denominator)
      if (!is.finite(sigma_new) || sigma_new <= 0) {
        stop("Invalid sigma estimate.", call. = FALSE)
      }

      model_df <- sum(beta_mean != 0)
      logLik_current <-
        -N / 2 * log(2 * pi * sigma_new^2) -
        weighted_ss / (2 * n_sam * sigma_new^2)
      bic_current <- -2 * logLik_current + model_df * log(N)

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
      coef = beta_old,
      sigma = sigma_old,
      BIC = bic_current,
      LogLik = logLik_current,
      lambda = lambda,
      lambda_index = ll,
      iterations = iter,
      delta = delta,
      converged = converged,
      effective_sample_size = ess_summary_current
    )

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

  # ---------------------------------------------------------------------------
  # Final sigma update with n_sam_final draws
  # ---------------------------------------------------------------------------
  final_rows <- rep(seq_len(N), times = n_sam_final)
  resp_final <- resp[final_rows, , drop = FALSE]
  Z_final <- Fan[final_rows, , drop = FALSE]
  Q_final <- N * n_sam_final

  theta_final <- rnorm(
    Q_final,
    mean = theta_est_irt.mean,
    sd = proposal_sd
  )
  final_weights <- importance_weights(
    theta = theta_final,
    resp_q = resp_final,
    Z_q = Z_final,
    beta = best$coefficients,
    sigma = best$sigma,
    n_draws = n_sam_final
  )
  final_residuals <-
    theta_final - drop(Z_final %*% best$coefficients)
  final_denominator <- sum(final_weights$weights) - p
  if (!is.finite(final_denominator) || final_denominator <= 0) {
    stop("Invalid denominator in final sigma update.", call. = FALSE)
  }
  sigma_final <- sqrt(
    sum(final_weights$weights * final_residuals^2) / final_denominator
  )

  final_ess <- base::c(
    minimum = min(final_weights$ess),
    median = stats::median(final_weights$ess),
    mean = mean(final_weights$ess)
  )

  update_progress_bar(total_steps)
  close_progress_bar()

  if (isTRUE(verbose)) {
    message(
      "Selected lambda index: ", minBIC,
      " | lambda = ", signif(best$lambda, 6),
      " | BIC = ", signif(best$BIC, 6),
      " | final sigma = ", signif(sigma_final, 6)
    )
  }

  list(
    coefficients = best$coefficients,
    coef = best$coefficients,
    sigma = sigma_final,
    tuning_sigma = best$sigma,
    BIC = best$BIC,
    LogLik = best$LogLik,
    minBIC = minBIC,
    lambda = best$lambda,
    iterations = best$iterations,
    delta = best$delta,
    effective_sample_size = best$effective_sample_size,
    final_effective_sample_size = final_ess,
    Convergence = if (best$converged) "Converged" else "Did not converge",
    factor_scores = Uupdate,
    fit_score = Uupdate,
    hatU = hatU,
    Fan = Fan,
    all_results = results
  )
}


