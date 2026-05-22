q_num_NA <- function(a, d, c, theta, resp, Z, beta_gamma_t, sigma_t) {
  # pij
  pij <- c + (1 - c) * (1 / (1 + exp(-(as.matrix(theta) %*% a + rep(d, each = length(theta))))))

  # mask non NA
  mask <- !is.na(resp)

  lpij_irt <- matrix(0, nrow = nrow(pij), ncol = ncol(pij))  # 0
  lpij_irt[mask] <- resp[mask] * log(pij[mask]) + (1 - resp[mask]) * log(1 - pij[mask])

  pi_irt <- exp(rowSums(lpij_irt))  # summing over items

  p_reg <- as.vector(dnorm(theta, beta_gamma_t %*% t(Z), sigma_t))
  num <- p_reg * pi_irt
  return(num)
}
# Function to add new results to the sliding window
add_to_window <- function(new_result, results, window_size) {
  # Add the new result
  results <- c(results, list(new_result))

  # If the window exceeds the size, remove the oldest result
  if (length(results) > window_size) {
    results <- results[-1]  # Remove the first element
  }

  return(results)
}
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
#' @param n Integer. Sample size (number of persons).
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
#' @param resp_rep Matrix. Replicated/expanded responses used for Monte Carlo
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
#' @export
farlr_emm <- function(n,resp, parTab, K_hat, p, lambda_all, delta.criteria = 1e-3, iter.max = 500, n_sam = 50,  window.size = 50, theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, Uupdate = NA,hatU = NA, Fan = NA, main, verbose = TRUE){
  set.seed(234)
  results <- list()
  a <- parTab$slope
  d <- -parTab$slope * parTab$difficulty
  c <- parTab$guessing
  output <- if (verbose) stderr() else nullfile()
  cat(file = output, 'Fitting the latent regression model...\n')
  if (verbose) pb <- txtProgressBar(file = output, 0, length(lambda_all), style = 3)
  for (ll in 1:length(lambda_all)) {
    # initial value
    beta_gamma_old <- rep(0,K_hat+p)
    beta_gamma_t <- c()
    sigma_old <- 1

    delta <- 1
    iter <- 1
    theta_t <- 1
    cov_t <- c()
    cv.fit.em <- list()
    cv.fit.em.debias <- list()
    lambda <- lambda_all[ll]
    window_size <- window.size

    while ((delta > delta.criteria) && (iter < iter.max)) {
      # E-step : sample the theta

      theta_sample <- rnorm(n*n_sam, theta_est_irt.mean, (theta_est_irt.se + 0.2))
      q_num_sample <- q_num_NA(a,d,c, theta_sample,resp_rep,Z.em,beta_gamma_old,sigma_old)

      h_sample <- dnorm(theta_sample, theta_est_irt.mean, (theta_est_irt.se + 0.2))

      den_all <- q_num_sample/h_sample

      den_i <- (1/n_sam)*as.numeric(tapply(den_all, (seq_along(den_all) - 1) %% nrow(resp) + 1, sum))
      w_ik <- (1/den_i)*q_num_sample/h_sample

      factors <- rep(1,p)
      factors[main] <- 0
      cv.fit.em[[iter]] <- glmnet(Z.em, theta_sample, standardize=FALSE, family="gaussian",
                                  penalty.factor=c(rep(0,K_hat),factors), weights = w_ik, intercept = FALSE, lambda = lambda)
      #plot(cv.fit.em)
      coef_hat_em <- coef(cv.fit.em[[iter]],complete=TRUE)[-1]

      #lm ----
      cv.fit.em.debias[[iter]] <- lm(theta_sample~0+Z.em[,which(coef_hat_em!=0)], weights = w_ik)
      coef_hat_em_debias <- rep(0,p+K_hat)
      coef_hat_em_debias[which(coef_hat_em!=0)] <- coef(cv.fit.em.debias[[iter]])
      coef_hat_em_debias <- ifelse(is.na(coef_hat_em_debias), 0 ,coef_hat_em_debias )
      # cross validation
      # cv.fit
      beta_gamma_t <- add_to_window(coef_hat_em_debias, beta_gamma_t, window_size)

      #calculate mean
      matrix_beta_gamma <- do.call(cbind, beta_gamma_t)

      beta_gamma_means <- rowMeans(matrix_beta_gamma)

      cov <- vcov(cv.fit.em.debias[[iter]])
      fitted_values.em <- matrix(beta_gamma_means, nrow = 1, ncol = p+K_hat)%*%t(Z.em)
      residuals <- theta_sample - fitted_values.em
      WSSR <- sum(w_ik*(residuals^2))#
      #var <- posterior_var
      sigma_temp <- sqrt((WSSR)/(sum(w_ik)-p))

      #joint response ???hg
      bic <- -2*(-n/2*log(2*pi*sigma_temp^2)- 1/(2*n_sam*sigma_temp^2) * sum(w_ik*(theta_sample - Z.em%*%coef_hat_em_debias)^2)) + sum(coef_hat_em_debias!=0)*log(n)

      delta_s <- abs(sigma_old-sigma_temp)
      delta_b <- max(abs(beta_gamma_old - beta_gamma_means))
      #delta_t <- abs(theta_t-theta_debiased)
      delta <- max(delta_s, delta_b)
      sigma_old <- sigma_temp
      cov_t <- cov
      beta_gamma_old <- beta_gamma_means
      #theta_t <- theta_debiased
      iter <- iter+1
      #print(c(delta))
    }
    #print(ll)
    #print(iter)
    if (verbose) {setTxtProgressBar(pb, pb$getVal() + 1)}
    results[[ll]] <- list(beta_gamma_old,sigma_old,bic,cov_t)
  }
  minBIC <- which.min(unlist(lapply(results, function(x) x[[3]])))
  return(list(
    coefficients  = results[[minBIC]][[1]],
    sigma = results[[minBIC]][[2]],
    LogLik = results[[minBIC]][[3]],
    minBIC = minBIC,
    Convergence = ifelse(iter<iter.max,"Converged","Did not converge")
  ))
}


