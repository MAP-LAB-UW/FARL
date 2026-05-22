#' FARL: Factor-Augmented Regularized Latent Regression for Large-Scale Assessment
#'
#' FARL is developed to support large-scale assessment (LSA) analyses, in which
#' latent regression models are used to integrate students’ background information
#' and to generate plausible values (PVs) for secondary analysis. In LSA settings,
#' background variables are often high-dimensional and highly correlated, posing
#' substantial challenges for traditional latent regression approaches.
#'
#' FARL implements factor-augmented regularized latent regression (FARLR), an
#' innovative framework that jointly models common factors and idiosyncratic
#' components. Regularization is employed to select relevant idiosyncratic
#' predictors, yielding interpretable regression results while maintaining
#' congeniality and estimation stability.
#'
#' @section Latent regression models:
#' \itemize{
#'   \item \code{\link{Farlr_mml}} fits the factor-augmented regularized latent regression model via marginal maximum likelihood
#'     \itemize{
#'       \item Method 1: \code{FARLR_EMM}
#'       \item Method 2: \code{FARLR_Debias}
#'     }
#' }
#' @section Latent regression models:
#' \itemize{
#'  \item \code{\link{Dire_mml}} fits the DIRE latent regression model
#' }
#' @section Plausible value generation:
#' \itemize{
#'   \item \code{\link{Farlr_drawPVs}} generates plausible values under the FARLR framework
#' }
#'
#' @section Simulation and example data:
#' \itemize{
#'   \item Built-in example datasets illustrating FARLR estimation and PV generation
#'   \item Utility functions for simulation studies in large-scale assessment settings
#' }
#'
#' @section Plausible value generation:
#' \itemize{
#'   \item \code{\link{Dire_drawPVs}} generates plausible values under the Dire framework
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib FARL
#' @import Rcpp
#' @import Matrix
#' @import MASS
#' @import mirt
#' @import torch
#' @import mvQuad
#' @importFrom tibble tibble lst
#' @importFrom abind abind
## usethis namespace: end
NULL
