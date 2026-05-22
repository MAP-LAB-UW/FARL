library(Dire)
#' DIRE Marginal Maximum Likelihood Estimation
#'
#' Fits a Direct Item Regression Effects (DIRE) model using a marginal maximum
#' likelihood (MML) framework. The DIRE model extends standard item response
#' theory by allowing item parameters to depend directly on person-level
#' covariates, enabling the assessment of covariate-induced differential item
#' functioning (DIF). This function provides a unified estimation interface for
#' dichotomous and polytomous item responses, supports complex survey designs,
#' and allows flexible numerical integration and optimization options.
#'
#' @param formula Formula. A model formula specifying person-level covariates
#'   entering the DIRE model. The left-hand side is ignored; the right-hand side
#'   defines covariates with potential direct effects on item parameters.
#' @param stuItems Data frame. Long-format item response data containing item
#'   identifiers and response values for each individual.
#' @param stuDat Data frame. Person-level data frame containing covariates
#'   referenced in \code{formula}.
#' @param idVar Character. Name of the variable in \code{stuItems} identifying
#'   individuals.
#' @param dichotParamTab Data frame or \code{NULL}. Item parameter table for
#'   dichotomous items. If \code{NULL}, dichotomous items are not modeled.
#' @param polyParamTab Data frame or \code{NULL}. Item parameter table for
#'   polytomous items. If \code{NULL}, polytomous items are not modeled.
#' @param testScale Character vector or \code{NULL}. Optional specification of
#'   test or scale membership for items.
#' @param Q Integer. Number of quadrature nodes used for numerical integration.
#'   Defaults to \code{30}.
#' @param minNode Numeric. Lower bound of the quadrature nodes. Defaults to
#'   \code{-4}.
#' @param maxNode Numeric. Upper bound of the quadrature nodes. Defaults to
#'   \code{4}.
#' @param polyModel Character. Polytomous item response model to be used.
#'   Supported options include \code{"GPCM"} and \code{"GRM"}.
#' @param weightVar Character or \code{NULL}. Optional sampling weight variable
#'   in \code{stuDat} used for weighted likelihood estimation.
#' @param multiCore Logical. If \code{TRUE}, enables parallel computation for
#'   likelihood evaluation where supported. Defaults to \code{FALSE}.
#' @param bobyqaControl List or \code{NULL}. Optional control parameters passed
#'   to the \code{bobyqa} optimizer for numerical optimization.
#' @param composite Logical. If \code{TRUE}, uses a composite likelihood
#'   approximation for estimation. Defaults to \code{TRUE}.
#' @param strataVar Character or \code{NULL}. Optional stratification variable
#'   for complex survey designs.
#' @param PSUVar Character or \code{NULL}. Optional primary sampling unit (PSU)
#'   variable for complex survey designs.
#' @param fast Logical. If \code{TRUE}, uses computational shortcuts to accelerate
#'   estimation. Defaults to \code{TRUE}.
#' @param calcCor Logical. If \code{TRUE}, computes correlation matrices for
#'   estimated item effects. Defaults to \code{TRUE}.
#' @param verbose Integer. Verbosity level controlling diagnostic output.
#'   \code{0} suppresses output; larger values produce more detailed messages.
#'
#' @return A list containing estimation results from the DIRE model, typically
#' including:
#' \describe{
#'   \item{\code{item.par}}{Estimated baseline item parameters.}
#'   \item{\code{dire.coef}}{Estimated direct item regression effect coefficients
#'   associated with person-level covariates.}
#'   \item{\code{vcov}}{Estimated variance–covariance matrix of parameter
#'   estimates.}
#'   \item{\code{LogLik}}{Maximized marginal (or composite) log-likelihood value.}
#'   \item{\code{convergence}}{Indicator of convergence status of the numerical
#'   optimization.}
#' }
#'
#' @details
#' The DIRE model allows person-level covariates to enter item parameter models
#' directly, providing a flexible framework for assessing covariate-related DIF
#' without requiring anchor items. Latent variables are integrated out using
#' Gaussian quadrature with \code{Q} nodes over the interval
#' [\code{minNode}, \code{maxNode}]. For large-scale assessments or complex survey
#' data, composite likelihood and survey design adjustments can be employed to
#' improve computational feasibility.
#'
#' @seealso
#' \code{\link[mirt]{mirt}}, \code{\link[lme4]{lmer}}, \code{\link[stats]{optim}}
#'
#' @examples
#' \dontrun{
#' fit <- Dire_mml(
#'   formula = ~ gender + ses,
#'   stuItems = stuItems,
#'   stuDat = stuDat,
#'   idVar = "student_id",
#'   dichotParamTab = dichotTab,
#'   polyParamTab = polyTab,
#'   polyModel = "GPCM",
#'   Q = 30,
#'   verbose = 1
#' )
#' }
#'
#' @export
Dire_mml <- function(formula, stuItems, stuDat, idVar, dichotParamTab = NULL,
                     polyParamTab = NULL, testScale = NULL, Q = 30, minNode = -4,
                     maxNode = 4, polyModel = c("GPCM", "GRM"), weightVar = NULL,
                     multiCore = FALSE, bobyqaControl = NULL, composite = TRUE,
                     strataVar = NULL, PSUVar = NULL, fast = TRUE, calcCor = TRUE,
                     verbose = 0){
  Dire::mml(formula, stuItems, stuDat, idVar, dichotParamTab,
            polyParamTab, testScale, Q, minNode,
            maxNode, polyModel , weightVar,
            multiCore, bobyqaControl, composite,
            strataVar, PSUVar, fast, calcCor,
            verbose)

}
#' Draw Plausible Values from a Fitted DIRE Model
#'
#' Generates plausible values (PVs) of latent traits from a fitted
#' Direct Item Regression Effects (DIRE) model. Plausible values are
#' random draws from the posterior distribution of the latent variable
#' given observed responses and estimated model parameters, and are
#' commonly used for secondary analyses to properly account for
#' measurement uncertainty.
#'
#' @param x Object. A fitted DIRE model object returned by
#'   \code{\link{Dire_mml}}.
#' @param npv Integer. Number of plausible values to draw for each
#'   individual.
#' @param pvVariableNameSuffix Character. Suffix appended to the names of
#'   generated plausible value variables. Defaults to \code{"_dire"}.
#' @param ... Additional arguments passed to internal sampling routines.
#'
#' @return A data frame containing the generated plausible values. Each
#' row corresponds to an individual, and each column corresponds to a
#' plausible value draw. Column names follow the pattern
#' \code{PV1<pvVariableNameSuffix>}, \code{PV2<pvVariableNameSuffix>}, …,
#' \code{PVnpv<pvVariableNameSuffix>}.
#'
#' @details
#' Plausible values are drawn by sampling from the posterior distribution
#' of the latent trait implied by the fitted DIRE model, conditional on
#' observed item responses and estimated item and regression parameters.
#' The resulting plausible values can be used in downstream analyses
#' (e.g., regression, group comparisons, or secondary modeling) by
#' combining results across draws using Rubin’s rules or other multiple
#' imputation techniques.
#'
#' This function is intended for post-estimation use and does not refit
#' the DIRE model.
#'
#' @seealso
#' \code{\link{Dire_mml}}, \code{\link[mirt]{fscores}}
#'
#' @examples
#' \dontrun{
#' fit <- Dire_mml(
#'   formula = ~ gender + ses,
#'   stuItems = stuItems,
#'   stuDat = stuDat,
#'   idVar = "student_id",
#'   dichotParamTab = dichotTab
#' )
#'
#' pv <- Dire_drawPVs(
#'   x = fit,
#'   npv = 5
#' )
#' }
#'
#' @export
Dire_drawPVs <- function(x, npv, pvVariableNameSuffix = "_dire", ...){
  Dire::drawPVs(x, npv, pvVariableNameSuffix = "_dire", ...)
}

