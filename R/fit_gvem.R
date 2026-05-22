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
Farlr_mml <- function(X, Y, parTab, n_sam = 15, method = "FARLR_EMM", lambda = seq(0.1, 0.5, by = 0.1), main, delta.criteria = 1e-3,iter.max = 500, window.size = 50, verbose = TRUE) {
  PA <- paran(X, iterations = 500, centile = 0, quiet = TRUE)
  K_hat <- PA$Retained
  fa <- factor.analysis(X, K_hat, method = "ml")
  Wupdate.t <- fa$Gamma
  Sgm_inv <- solve(diag(fa$Sigma))
  p = ncol(X)
  N = nrow(X)
  orthg <- t(Wupdate.t) %*% Sgm_inv %*% Wupdate.t/p
  V <- eigen(orthg)$vectors
  Wupdate <- t(Wupdate.t %*% V)
  Uupdate <- t(solve(Wupdate %*% Sgm_inv %*% t(Wupdate)) %*% Wupdate %*% Sgm_inv %*% t(X))
  Z <- cbind(Uupdate,X)
  colnames(Y) <- paste0("item", c(1:ncol(Y)), sep = "")
  storage.mode(Y) <- "integer"
  est_mirt <- suppressMessages(
    suppressWarnings(
      mirt(Y, 1, verbose = FALSE)
    )
  )
  theta_est_irt <- fscores(est_mirt,full.scores.SE = TRUE)
  theta_est_irt.mean <- theta_est_irt[,1]
  theta_est_irt.se <- theta_est_irt[,2]
  main <- c(1, 2, 29, 15, 45)
  itemNames <- if (is.null(colnames(Y))) {
    paste0("i", sprintf("%03d", seq_len(ncol(Y))))
  } else {
    colnames(Y)
  }
  colnames(Y) <- itemNames
  subject <- factor(c(1:nrow(X)))
  stuItems <- reshape(data=data.frame(cbind(Y,subject)), varying=itemNames, idvar="subject",
                      direction="long", v.names="score",
                      times=itemNames, timevar="key")
  new_itemNames <- (paste0("item",1:ncol(Y)))
  stuItems$key <- rep(new_itemNames,each=n)
  stuDat <- X
  subject <- factor(c(1:n))
  stuDat <- data.frame(cbind(subject, X))
  colnames(stuDat) <-c("subject",paste("X", c(1:p), sep = ""))

  Y_back <- reshape(
    stuItems,
    idvar = "subject",
    timevar = "key",
    direction = "wide"
  )
  rownames(Y_back) <- subject
  Y_back <- Y_back[,-1]
  if (method == "FARLR_Debias") {
    hatB <- t( (1/n) * t(Uupdate) %*% X )
    hatU <- X - Uupdate %*% t(hatB)
    Fan  <- cbind(Uupdate, hatU)
    Z.em <- NA
    resp_rep <- NA
    fn <- farlr_debias

  } else if (method == "FARLR_EMM") {
    hatB <- NA
    hatU <- NA
    Fan  <- NA
    resp_rep <- rep(1, n_sam) %x% Y
    Z.em <- rep(1, n_sam) %x% Z
    fn <- farlr_emm
  } else {
    stop("Unknown method: ", method)
  }

  result <-fn(
    N = nrow(resp),
    resp = Y,
    resp.em = resp_rep,
    a = a,
    d = d,
    c = c,
    b1 = b1,
    b2 = b2,
    type = type,
    K_hat = K_hat,
    p = ncol(X),
    lambda_all = seq(0.001, 0.1, length.out = 10),
    delta.criteria = 1e-3,
    iter.max = 500,
    n_sam = n_sam,
    window.size = window.size,
    theta_est_irt.mean = theta_mean,
    theta_est_irt.se = theta_se,
    Uupdate = Uupdate,
    hatU = hatU,
    Fan = Fan,
    main = main,
    Fan.em = Fan.em,
    verbose = TRUE
  )
  #resultII <- fn(nrow(Y), Y_back, parTab, K_hat, ncol(X), lambda, delta.criteria = 1e-3,iter.max = 200, n_sam = 30, window.size = 50,theta_est_irt.mean, theta_est_irt.se, resp_rep, Z.em, bin, verbose = TRUE)
  result$stuDat <- cbind(subject,Z)
  result$stuItems <- stuItems
  colnames(X) <- paste0("X", c(1:ncol(sim_a1$X)))
  result$X <- X
  result$item_params <- parTab
  invisible(result)
}
mml_test <- function(){
  mmlcomp <-  with(sim_a1, Farlr_mml(X, Y, parTab, method = "FARLR_Debias")) #mml(X = sim_a1$X, Y = sim_a1$Y, parTab = sim_a1$parTab,  method = "FARLR_Debias")
  invisible(mmlcomp)
  PVs <- drawPVs(mmlcomp, 10L)
  return(PVs)
}
mml_test2 <- function(){
  with(sim_a1, Farlr_mml(X, Y, parTab, method = "FARLR_EMM"))
  }
