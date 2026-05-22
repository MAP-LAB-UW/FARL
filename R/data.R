#' Simulated Dataset: 1D FARLR-Style Dichotomous Item Responses
#'
#' A simulated one-dimensional dichotomous response dataset generated under a
#' FARLR-style latent regression formulation. The latent trait \eqn{\theta} is
#' constructed from a linear component plus scaled noise to achieve a target
#' signal-to-noise ratio (SNR), and is then centered and standardized.
#' Item responses are subsequently generated using \code{mirt::simdata()} with
#' dichotomous (2PL) items.
#'
#' @author Yijun Cheng <chengxb@uw.edu>
#'
#' @format A list with the following components:
#' \tabular{ll}{
#' \code{X}        \tab Covariate/design matrix used for generating \eqn{\theta}.\cr
#' \code{Y}        \tab Simulated dichotomous item response matrix (\code{N} by \code{J}).\cr
#' \code{a}        \tab True item discrimination parameters (slopes), length \code{J}.\cr
#' \code{b}        \tab True item difficulty parameters, length \code{J}.\cr
#' \code{d}        \tab True item intercept parameters used by \code{mirt::simdata()},
#'                  computed as \code{d = -a*b}.\cr
#' \code{parTab}   \tab Parameter table used internally for model fitting.\cr
#' }
#'
#' @details
#' The latent trait is generated from a linear predictor and additive noise:
#' \itemize{
#'   \item \eqn{\theta = F \beta + E \nu + e}, where \eqn{e} is scaled to match a target SNR.
#' }
#'
#' Item parameters are generated as:
#' \itemize{
#'   \item \eqn{a_j \sim \mathrm{Lognormal}(0, 0.25)}
#'   \item \eqn{b_j \sim \mathrm{Uniform}(-2, 2)}
#'   \item \eqn{d_j = -a_j b_j}
#' }
#' Responses are simulated via \code{mirt::simdata(itemtype = "dich")}.
#'
#' @usage data(sim_a1)
#'
"sim_a1"
