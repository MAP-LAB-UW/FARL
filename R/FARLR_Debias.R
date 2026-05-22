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
#' @export
farlr_debias <- function(N, resp, resp.em, a, d, c, b1, b2, type,
                              K_hat, p, lambda_all,
                              delta.criteria = 1e-3,
                              iter.max = 500,
                              n_sam = 5,
                              window.size = 50,
                              theta_est_irt.mean,
                              theta_est_irt.se,
                              Uupdate,
                              hatU,
                              Fan,
                              main,
                              Fan.em,
                              verbose = TRUE) {

  results <- vector("list", length(lambda_all))

  # Expanded matrices used during tuning
  Uupdate.em <- rep(1, n_sam) %x% Uupdate
  hatU.em <- rep(1, n_sam) %x% hatU
  Z.em <- Fan.em

  Uupdate.em.torch <- torch_tensor(Uupdate.em)
  hatU.em.torch <- torch_tensor(hatU.em)

  Ut_update <- Uupdate.em.torch$transpose(1, 2)
  inv_nnsam <- 1 / (N * n_sam)

  for (ll in seq_along(lambda_all)) {

    beta_gamma_old <- rep(0, K_hat + p)
    sigma_old <- 1
    beta_gamma_t <- list()

    delta <- 1
    iter <- 1
    lambda <- lambda_all[ll]

    if (verbose) {
      message("Starting lambda ", ll, " of ", length(lambda_all),
              " | lambda = ", lambda)
    }

    while (delta > delta.criteria && iter < iter.max) {

      # E-step: sample latent proficiency values
      theta_sample <- rnorm(
        N * n_sam,
        mean = theta_est_irt.mean,
        sd = theta_est_irt.se + 0.2
      )

      theta_sample.torch <- torch_tensor(theta_sample)

      q_num_sample <- q_num_NA(
        a, d, c, b1, b2, type,
        theta_sample, resp.em, Z.em,
        beta_gamma_old, sigma_old
      )

      h_sample <- dnorm(
        theta_sample,
        mean = theta_est_irt.mean,
        sd = theta_est_irt.se + 0.2
      )

      den_all <- q_num_sample / h_sample

      # Compute importance weights by individual
      den_mat <- matrix(den_all, nrow = N, ncol = n_sam)
      den_i <- rowMeans(den_mat)
      w_mat <- den_mat / den_i
      w_ik <- as.vector(w_mat)

      w_ik.torch <- torch_tensor(w_ik)

      # Compute theta_1 without forming the large projection matrix
      proj <- inv_nnsam *
        Uupdate.em.torch$matmul(
          Ut_update$matmul(theta_sample.torch$unsqueeze(2))
        )$squeeze(2)

      theta_1 <- theta_sample.torch - proj

      # Weighted LASSO fit
      penalty.factor <- rep(1, p)
      penalty.factor[main] <- 0

      fit <- glmnet(
        x = hatU.em,
        y = as.matrix(theta_1),
        weights = w_ik,
        penalty.factor = penalty.factor,
        intercept = FALSE,
        lambda = lambda
      )

      coef_hat_em <- coef(fit, complete = TRUE)[-1]
      coef_hat_em.torch <- torch_tensor(coef_hat_em)

      # Debiasing correction
      weighted_sq <- w_ik.torch$unsqueeze(2) * hatU.em.torch$pow(2)
      T <- inv_nnsam * weighted_sq$sum(dim = 1)

      Theta1 <- torch_diag(1 / T)
      temp <- inv_nnsam * torch_matmul(Theta1, hatU.em.torch$t())

      residuals <- theta_1 - hatU.em.torch$matmul(coef_hat_em.torch)
      w_resid <- w_ik.torch * residuals

      coef_hat_debias <- coef_hat_em.torch +
        temp$matmul(w_resid$unsqueeze(2))$squeeze(2)

      coef_hat_debias <- as.array(coef_hat_debias)

      # Estimate factor coefficients using weighted least squares
      Uw <- Uupdate.em * w_ik
      lhs <- crossprod(Uupdate.em, Uw)
      rhs <- crossprod(Uupdate.em, theta_sample * w_ik)
      phi_hat <- as.numeric(solve(lhs, rhs))

      # Keep only coefficients selected by the penalized fit
      coef_hat_em_debias <- ifelse(coef_hat_em != 0, coef_hat_debias, 0)
      coef_hat_em_debias <- append(coef_hat_em_debias, phi_hat, after = 0)

      # Smooth coefficient updates using a sliding window
      beta_gamma_t <- add_to_window(
        coef_hat_em_debias,
        beta_gamma_t,
        window.size
      )

      beta_gamma_means <- rowMeans(do.call(cbind, beta_gamma_t))

      fitted_values.em <- matrix(beta_gamma_means, nrow = 1) %*% t(Z.em)
      residuals_2 <- theta_sample - fitted_values.em

      WSSR <- sum(w_ik * residuals_2^2)
      sigma_means <- sqrt(WSSR / (nrow(Z.em) - p))

      bic <- -2 * (
        -N / 2 * log(2 * pi * sigma_means^2) -
          1 / (2 * n_sam * sigma_means^2) *
          sum(w_ik * (theta_sample - Z.em %*% beta_gamma_means)^2)
      ) + sum(beta_gamma_means != 0) * log(N)

      delta_s <- abs(sigma_old - sigma_means)
      delta_b <- max(abs(beta_gamma_old - beta_gamma_means))
      delta <- max(delta_s, delta_b)

      sigma_old <- sigma_means
      beta_gamma_old <- beta_gamma_means
      iter <- iter + 1

      if (verbose) {
        message("  iter = ", iter,
                " | delta = ", signif(delta, 4),
                " | sigma = ", signif(sigma_old, 4),
                " | BIC = ", signif(bic, 4))
      }
    }

    if (verbose) {
      message("Finished lambda ", ll,
              " | iterations = ", iter,
              " | final delta = ", signif(delta, 4),
              " | BIC = ", signif(bic, 4))
    }

    results[[ll]] <- list(
      coef = beta_gamma_old,
      sigma = sigma_old,
      bic = bic
    )
  }

  # Select lambda by minimum BIC
  minBIC <- which.min(sapply(results, function(x) x$bic))
  beta_gamma_old <- results[[minBIC]]$coef
  sigma_old <- results[[minBIC]]$sigma

  if (verbose) {
    message("Selected lambda index: ", minBIC,
            " | lambda = ", lambda_all[minBIC])
  }

  # Recompute sigma using a larger Monte Carlo sample
  n_sam_final <- 60

  resp_rep <- rep(1, n_sam_final) %x% resp
  Z.em <- rep(1, n_sam_final) %x% Fan

  theta_sample <- rnorm(
    N * n_sam_final,
    mean = theta_est_irt.mean,
    sd = theta_est_irt.se + 0.2
  )

  q_num_sample <- q_num_NA(
    a, d, c, b1, b2, type,
    theta_sample, resp_rep, Z.em,
    beta_gamma_old, sigma_old
  )

  h_sample <- dnorm(
    theta_sample,
    mean = theta_est_irt.mean,
    sd = theta_est_irt.se + 0.2
  )

  den_all <- q_num_sample / h_sample

  den_i <- (1 / n_sam_final) * as.numeric(
    tapply(den_all, (seq_along(den_all) - 1) %% N + 1, sum)
  )

  w_ik <- den_all / den_i

  fitted_values.em <- matrix(beta_gamma_old, nrow = 1) %*% t(Z.em)
  residuals_2 <- theta_sample - fitted_values.em

  WSSR <- sum(w_ik * residuals_2^2)
  sigma_final <- sqrt(WSSR / (nrow(Z.em) - p))

  if (verbose) {
    message("Final sigma = ", signif(sigma_final, 4))
  }

  return(list(
    coef = beta_gamma_old,
    sigma = sigma_final,
    minBIC = minBIC
  ))
}


