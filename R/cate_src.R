## This file contains all the factor analysis methods.
## R comes with a factor analysis function factanal()
## There is also fa() in the package psych.
## Here we implement some recent algorithms in the case of n < p and build a high level wrapper of them.
## For p >> n, empirical results show that using all three methods give similar results

#' Factor analysis
#'
#' The main function for factor analysis with potentially high dimensional variables.
#' Here we implement some recent algorithms that is optimized for the high dimensional problem where
#' the number of samples n is less than the number of variables p.
#'
#' @param Y data matrix, a n*p matrix
#' @param r number of factors
#' @param method algorithm to be used
#'
#' @details The three methods are quasi-maximum likelihood (ml),
#' principal component analysis (pc),
#' and factor analysis using an early stopping criterion (esa).
#'
#' The ml is iteratively solved the Expectation-Maximization algorithm
#' using the PCA solution as the initial value.
#' See Bai and Li (2012) and for more details. For the esa method, see
#' Owen and Wang (2015) for more details.
#'
#' @return a list of objects
#' \describe{
#' \item{Gamma}{estimated factor loadings}
#' \item{Z}{estimated latent factors}
#' \item{Sigma}{estimated noise variance matrix}
#' }
#'
#' @references {
#' Bai, J. and Li, K. (2012). Statistical analysis of factor models of high dimension. \emph{The Annals of Statistics 40}, 436-465.
#' Owen, A. B. and Wang, J. (2015). Bi-cross-validation for factor analysis. \emph{arXiv:1503.03515}.
#' }
#'
#' @examples
#' ## a factor model
#' n <- 100
#' p <- 1000
#' r <- 5
#' Z <- matrix(rnorm(n * r), n, r)
#' Gamma <- matrix(rnorm(p * r), p, r)
#' Y <- Z %*% t(Gamma) + rnorm(n * p)
#'
#' ## to check the results, verify the true factors are in the linear span of the estimated factors.
#' pc.results <- factor.analysis(Y, r = 5, "pc")
#' sapply(summary(lm(Z ~ pc.results$Z)), function(x) x$r.squared)
#'
#' ml.results <- factor.analysis(Y, r = 5, "ml")
#' sapply(summary(lm(Z ~ ml.results$Z)), function(x) x$r.squared)
#'
#' esa.results <- factor.analysis(Y, r = 5, "esa")
#' sapply(summary(lm(Z ~ esa.results$Z)), function(x) x$r.squared)
#'
#' @seealso \code{\link{fa.pc}}, \code{\link{fa.em}}, \code{\link[esaBcv]{ESA}}
#'
#' @import esaBcv
#' @keywords internal
#' @noRd
#'
factor.analysis <- function(Y,
                            r,
                            method = c("ml", "pc", "esa")) {

  # match the arguments
  if (r == 0) {
    return(list(Gamma = NULL,
                Z = NULL,
                Sigma = apply(Y, 2, function(v) mean(v^2))))
  }
  method <- match.arg(method, c("ml", "pc", "esa"))

  if (method == "pc") {
    fa.pc(Y, r)
  } else if (method == "ml") {
    fa.em(Y, r)
  } else if(method == "esa") {
    result <- ESA(Y, r)
    return(list(Gamma = result$estV %*% diag(result$estD, r,r),
                Z = sqrt(nrow(Y)) * result$estU,
                Sigma = result$estSigma))
  }
}
#' @keywords internal
fa.pc <- function(Y, r) {

  svd.Y <- svd(Y)

  Gamma <- svd.Y$v[, 1:r] %*% diag(svd.Y$d[1:r], r, r) / sqrt(nrow(Y))
  Z <- sqrt(nrow(Y)) * svd.Y$u[, 1:r]

  Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) mean(x^2))

  return(list(Gamma = Gamma,
              Z = Z,
              Sigma = Sigma))

}
#' @keywords internal
fa.em <- function(Y, r, tol = 1e-6, maxiter = 1000) {

  ## A Matlab version of this EM algorithm was in http://www.mathworks.com/matlabcentral/fileexchange/28906-factor-analysis/content/fa.m

  ## The EM algorithm:
  ##
  ## Y = Z Gamma' + E Sigma^{1/2}                        (n * p)
  ## mle to estimate Gamma (p * r) and Sigma             (p * p)

  ## E step:
  ## EZ = Y (Gamma Gamma' + Sigma)^{-1} Gamma            (n * r)
  ## VarZ = I - Gamma' (Gamma Gamma' + Sigma)^{-1} Gamma (r * r)
  ## EZ'Z = n VarZ + EZ' * EZ                            (r * r)

  ## M step:
  ## update Gamma: Gamma = (Y' EZ)(EZ'Z)^{-1}
  ## update Sigma:
  ## Sigma = 1/n diag(Y'Y - Gamma EZ' Y - Y' EZ' Gamma' + Gamma EZ'Z Gamma')

  ## The log-likelihood (simplified)
  ## llh <- -log det(Gamma Gamma' + Sigma) - tr((Gamma Gamma' + Sigma)^{-1} S)
  ##

  ## For details, see http://cs229.stanford.edu/notes/cs229-notes9.pdf

  p <- ncol(Y)
  n <- nrow(Y)

  ## initialize parameters
  init <- fa.pc(Y, r)
  Gamma <- init$Gamma
  #Gamma <- matrix(runif(p * r), nrow = p)
  #invSigma <- 1/colMeans(Y^2)
  invSigma <- 1/init$Sigma
  #invSigma <- 1/colMeans((Y - sqrt(n) * start.svd$u %*% t(Gamma))^2)
  llh <- -Inf # log-likelihood

  ## precompute quantitites
  I <- diag(rep(1, r))
  ## diagonal of the sample Covariance S
  #sample.var <- apply(Y, 2, var)
  sample.var <- colMeans(Y^2)

  ## compute quantities needed
  tilde.Gamma <- sqrt(invSigma) * Gamma
  M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
  eigenM <- eigen(M, symmetric = TRUE)
  YSG <- Y %*% (invSigma * Gamma)

  logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
  B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
  logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
  llh <- -logdetY - logtrY

  converged <- FALSE
  for (iter in 1:maxiter) {

    ## E step:
    ## Using Woodbury matrix identity:
    ## tilde.Gamma = Sigma^{-1/2} Gamma
    ## VarZ = (I + tilde.Gamma' tilde.Gamma)^{-1}
    ## EZ = Y Sigma^{-1} Gamma VarZ
    varZ <- eigenM$vectors %*% (1/eigenM$values * t(eigenM$vectors))
    EZ <- YSG %*% varZ
    EZZ <- n * varZ + t(EZ) %*% EZ

    ## M step:
    eigenEZZ <- eigen(EZZ, symmetric = TRUE)
    YEZ <- t(Y) %*% EZ
    ## EZ'Z = G'G
    G <- sqrt(eigenEZZ$values) * t(eigenEZZ$vectors)
    ## updating invSigma
    invSigma <- 1/(sample.var - 2/n * rowSums(YEZ * Gamma) +
                     1/n * rowSums((Gamma %*% t(G))^2))
    ## updating Gamma
    Gamma <- YEZ %*% eigenEZZ$vectors %*%
      (1/eigenEZZ$values * t(eigenEZZ$vectors))

    ## compute quantities needed
    tilde.Gamma <- sqrt(invSigma) * Gamma
    M <- diag(r) + t(tilde.Gamma) %*% tilde.Gamma
    eigenM <- eigen(M, T)
    YSG <- Y %*% (invSigma * Gamma)

    ## compute likelihood and check for convergence
    old.llh <- llh
    ## Ussing Woodbury matrix identity
    ## log det(Gamma Gamma' + Sigma) =
    ## log [det(invSigma^{-1})det(M)]
    ## tr((Gamma Gamma' + Sigma)^{-1} S) =
    ## tr(invSigma S) - tr(B'B)/n
    ## where B = (M')^{-1/2} YSG'
    logdetY <- -sum(log(invSigma)) + sum(log(eigenM$values))
    B <- 1/sqrt(eigenM$values) * t(eigenM$vectors) %*% t(YSG)
    logtrY <- sum(invSigma * sample.var) - sum(B^2)/n
    llh <- -logdetY - logtrY

    if (abs(llh - old.llh) < tol * abs(llh)) {
      converged <- TRUE
      break
    }
  }

  # GLS to estimate factor loadings
  svd.H <- svd(t(Gamma) %*% (invSigma * Gamma))
  Z <- Y %*% (invSigma * Gamma) %*% (svd.H$u %*% (1/svd.H$d * t(svd.H$v)))

  return(list(Gamma = Gamma,
              Sigma = 1/invSigma,
              Z = Z,
              niter = iter,
              converged = converged))

}
#' @keywords internal
est.factor.num <- function(Y,
                           method = c("bcv", "ed"),
                           rmax = 20,
                           nRepeat = 12,
                           bcv.plot = TRUE, log = "") {

  method <- match.arg(method, c("bcv", "ed"))

  n <- nrow(Y)
  p <- ncol(Y)

  if (identical(method, "bcv")) {
    result <- EsaBcv(Y, nRepeat = nRepeat, r.limit = rmax)
    r <- result$best.r
    avg.sigma2 <- mean(result$estSigma)
    errors <- colMeans(result$result.list)/avg.sigma2
    if (bcv.plot) {
      plot(0:(length(errors)-1), errors, log = log, xlab = "r",
           ylab = "bcv MSE relative to the noise", type = "o")
    }
    return(list(r = r, errors = errors))
  } else if (identical(method, "ed")) {
    return(EigenDiff(Y, rmax))
  }

}
#' @keywords internal
EigenDiff <- function(Y, rmax = 20, niter = 10) {
  n <- nrow(Y)
  p <- ncol(Y)
  ev <- svd(Y)$d^2 / n
  n <- length(ev)
  if (is.null(rmax))
    rmax <- 3 * sqrt(n)
  j <- rmax + 1

  diffs <- ev - c(ev[-1], 0)

  for (i in 1:niter) {
    y <- ev[j:(j+4)]
    x <- ((j-1):(j+3))^(2/3)
    lm.coef <- lm(y ~ x)
    delta <- 2 * abs(lm.coef$coef[2])
    idx <- which(diffs[1:rmax] > delta)
    if (length(idx) == 0)
      hatr <- 0
    else hatr <- max(idx)

    newj = hatr + 1
    if (newj == j) break
    j = newj
  }

  return(hatr)

}

# q_num_NA <- function(a, d, c, theta, resp, Z, beta_gamma_t, sigma_t) {
#   # pij
#   pij <- c + (1 - c) * (1 / (1 + exp(-(as.matrix(theta) %*% a + rep(d, each = length(theta))))))
#
#   # mask non NA
#   mask <- !is.na(resp)
#
#   lpij_irt <- matrix(0, nrow = nrow(pij), ncol = ncol(pij))  # 0
#   lpij_irt[mask] <- resp[mask] * log(pij[mask]) + (1 - resp[mask]) * log(1 - pij[mask])
#
#   pi_irt <- exp(rowSums(lpij_irt))  # summing over items
#
#   p_reg <- as.vector(dnorm(theta, beta_gamma_t %*% t(Z), sigma_t))
#   num <- p_reg * pi_irt
#   return(num)
# }
#' @keywords internal
q_num_NA <- function(a, d, c, b1, b2, type, theta, resp, Z, beta_gamma_t, sigma_t) {

  Q <- length(theta)
  J <- length(type)

  # resp: Q x J (theta-repeated); NA allowed
  mask <- !is.na(resp)

  idx_g <- (type == "gpcm")
  idx_d <- !idx_g

  pij <- matrix(NA_real_, nrow = Q, ncol = J)

  # ---- (1) non-GPCM (2PL/3PL): pij = c + (1-c)*logit^{-1}(a*theta + d) ----
  if (any(idx_d)) {
    eta <- as.matrix(theta) %*% a[idx_d] + rep(d[idx_d], each = Q) # Q x Jd
    pij[, idx_d] <- c[idx_d] + (1 - c[idx_d]) * plogis(eta)
  }

  # ---- (2) GPCM (0/1/2): pij = P(Y = observed category | theta) ----
  if (any(idx_g)) {
    ag <- a[idx_g]
    b1g <- b1[idx_g]
    b2g <- b2[idx_g]
    yg <- resp[, idx_g] # Q x Jg

    # Q x Jg
    e1 <- exp( outer(theta, ag, "*") - matrix(rep(ag * b1g, each = Q), nrow = Q) )
    e2 <- exp( 2*outer(theta, ag, "*") - matrix(rep(ag * (b1g + b2g), each = Q), nrow = Q) )
    D <- 1 + e1 + e2

    P0 <- 1 / D
    P1 <- e1 / D
    P2 <- e2 / D

    pij[, idx_g] <- P0 * (yg == 0) + P1 * (yg == 1) + P2 * (yg == 2)
  }

  # numeric safety (won't affect NA)
  pij <- pmin(pmax(pij, 1e-12), 1 - 1e-12)

  # ---- log-likelihood (mask-based, your style) ----
  lpij_irt <- matrix(0, nrow = Q, ncol = J)

  # non-gpcm: Bernoulli log-likelihood on observed cells
  mask_d <- mask & matrix(rep(idx_d, each = Q), nrow = Q)
  lpij_irt[mask_d] <- resp[mask_d] * log(pij[mask_d]) +
    (1 - resp[mask_d]) * log(1 - pij[mask_d])

  # gpcm: log prob of observed category on observed cells
  mask_g <- mask & matrix(rep(idx_g, each = Q), nrow = Q)
  lpij_irt[mask_g] <- log(pij[mask_g])

  pi_irt <- exp(rowSums(lpij_irt)) # length Q

  # prior p(theta | Z)
  mu <- as.numeric(beta_gamma_t %*% t(Z))
  p_reg <- dnorm(theta, mean = mu, sd = sigma_t)

  num <- p_reg * pi_irt
  return(num)
}
# Function to add new results to the sliding window
#' @keywords internal
add_to_window <- function(new_result, results, window_size) {
  # Add the new result
  results <- c(results, list(new_result))

  # If the window exceeds the size, remove the oldest result
  if (length(results) > window_size) {
    results <- results[-1]  # Remove the first element
  }

  return(results)
}
# drawPVs(mmlcomp, num_of_PVs) <- function(){
#   return(draw_PVs_manual_merge(
#     nrow(mmlcomp$X),
#     mmlcomp$stuItems$score,
#     mmlcomp$coef,
#     mmlcomp$sigma,
#     fit_score,
#     scale(mmlcomp$X),
#     mmlcomp$X,
#     mmlcomp$item_params$a,
#     mmlcomp$item_params$d,
#     mmlcomp$item_params$c,
#     main,
#     step = 10,
#     theta,
#     X[,main]
#   ))
# }
# draw_PVs_manual_merge <- function(n, resp, beta_hat, sigma, U, X_p, X, a, d, c,
#                                   main, step, theta, X_discrete) {
#   # --------------------------------------------------
#   # Draw plausible values manually based on EAP
#   # and summarize them within binary subgroups
#   # --------------------------------------------------
#   # Quadrature grid
#   Q <- 30
#   t_q <- seq(-4, 4, length.out = Q)
#
#   # Store EAP posterior mean and variance for each examinee
#   EAP_estimates <- numeric(n)
#   EAP_variance  <- numeric(n)
#   Likelihood <- numeric(n)
#
#   # --------------------------------------------------
#   # Likelihood function for one examinee
#   # --------------------------------------------------
#   compute_likelihood <- function(response_vector, a, b, t_q) {
#     likelihood <- sapply(t_q, function(t) {
#       P <- c + (1 - c) * (1 / (1 + exp(-a * (t - b))))
#       prod((P^response_vector) * ((1 - P)^(1 - response_vector)))
#     })
#     return(likelihood)
#   }
#
#   # --------------------------------------------------
#   # Compute EAP mean and variance for each examinee
#   # --------------------------------------------------
#   for (i in 1:n) {
#     prior_mean <- c(as.numeric(U[i, ]), as.numeric(X_p[i, ])) %*% beta_hat
#     prior_f2 <- dnorm(t_q, mean = prior_mean, sd = sigma)
#
#     likelihood2 <- compute_likelihood(resp[i, ], a, -(d / a), t_q)
#
#     posterior_weight <- likelihood2 * prior_f2
#
#     EAP_estimates[i] <- sum(t_q * posterior_weight) / sum(posterior_weight)
#     EAP_variance[i]  <- sum((t_q - EAP_estimates[i])^2 * posterior_weight) /
#       sum(posterior_weight)
#     Likelihood[i] <- log(sum(likelihood2 * prior_f2 * (t_q[2] - t_q[1] )))
#   }
#   sum_ll <- sum(Likelihood)
#
#   # Compare estimated mean with true mean
#   true_mean <- mean(theta)
#   eap_mean  <- mean(EAP_estimates)
#   print(paste("Difference using estimate:", round(eap_mean - true_mean, 3)))
#
#   # --------------------------------------------------
#   # Draw plausible values from N(EAP, posterior variance)
#   # --------------------------------------------------
#   theta_pv <- numeric(n * step)
#
#   for (i in 1:n) {
#     for (j in 1:step) {
#       #set.seed(234)
#       theta_pv[(i - 1) * step + j] <- rnorm(
#         n = 1,
#         mean = EAP_estimates[i],
#         sd   = sqrt(EAP_variance[i])
#       )
#     }
#   }
#
#   theta_pv_matrix <- matrix(theta_pv, ncol = step, byrow = TRUE)
#
#   # --------------------------------------------------
#   # Repeat discrete covariates to match PV expansion
#   # --------------------------------------------------
#   X.rep <- vector("list", length(main))
#   for (i in 1:length(main)) {
#     X.rep[[i]] <- unlist(lapply(X_discrete[, i], function(x) rep(x, step)))
#   }
#
#   # --------------------------------------------------
#   # Create all binary combinations of the main variables
#   # --------------------------------------------------
#   n_main <- length(main)
#   combinations <- expand.grid(replicate(n_main, 0:1, simplify = FALSE))
#   colnames(combinations) <- main
#
#   # Store subgroup PV matrices
#   result_list <- vector("list", nrow(combinations))
#
#   # --------------------------------------------------
#   # Split PVs into subgroups based on binary patterns
#   # --------------------------------------------------
#   for (i in 1:nrow(combinations)) {
#     condition <- rep(TRUE, length(theta_pv))
#
#     for (j in 1:n_main) {
#       condition <- condition & (X.rep[[j]] == combinations[i, j])
#     }
#
#     subgroup_theta <- theta_pv[condition]
#     n_subjects <- sum(condition) / step
#
#     result_list[[i]] <- matrix(
#       subgroup_theta,
#       nrow = n_subjects,
#       ncol = step,
#       byrow = TRUE
#     )
#   }
#
#   # --------------------------------------------------
#   # Compute subgroup summaries
#   # --------------------------------------------------
#   mean_vector   <- numeric(length(result_list))
#   mean_v_vector <- numeric(length(result_list))
#   var_vector    <- numeric(length(result_list))
#
#   for (i in 1:length(result_list)) {
#     result <- result_list[[i]]
#     mean_vector[i]   <- mean(colMeans(result))
#     var_vector[i]    <- variance.sample(result)
#     mean_v_vector[i] <- variance.cal(result)
#   }
#
#   # --------------------------------------------------
#   # Return results
#   # --------------------------------------------------
#   PVs <- matrix(theta_pv, ncol = step, byrow = TRUE)
#
#   return(list(
#     mean   = mean_vector,
#     var    = var_vector,
#     mean_v = mean_v_vector,
#     PVs    = PVs,
#     group_PVs = result_list,
#     sum_ll = sum_ll
#   ))
# }


