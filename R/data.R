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

#' Simulated Dataset: 1D FARLR-Style Mixed-Format Item Responses
#'
#' A simulated one-dimensional mixed-format response dataset generated under a
#' FARLR-style latent regression formulation. The latent trait \eqn{\theta} is
#' constructed from a linear component plus scaled noise to achieve a target
#' signal-to-noise ratio (SNR), and is then centered and standardized. Item
#' responses are subsequently generated using \code{mirt::simdata()} with a
#' mixture of three-parameter logistic (3PL) and generalized partial credit
#' model (GPCM) items.
#'
#' @author Yijun Cheng <chengxb@uw.edu>
#'
#' @format A list with the following components:
#' \tabular{ll}{
#' \code{X}        \tab Covariate/design matrix used for generating
#'                       \eqn{\theta}.\cr
#' \code{Y}        \tab Simulated mixed-format item response matrix
#'                       (\code{N} by \code{J}). The 3PL items have scores 0/1,
#'                       while the GPCM items have scores 0/1/2.\cr
#' \code{a}        \tab True item discrimination parameters (slopes), length
#'                       \code{J}.\cr
#' \code{b}        \tab True difficulty parameters for 2PL/3PL items. Entries
#'                       corresponding to GPCM items are \code{NA}.\cr
#' \code{c}        \tab True lower-asymptote (guessing) parameters for 3PL
#'                       items. Entries for non-3PL items are zero.\cr
#' \code{b1}       \tab First GPCM step-difficulty parameter. Entries for
#'                       dichotomous items are \code{NA}.\cr
#' \code{b2}       \tab Second GPCM step-difficulty parameter. Entries for
#'                       dichotomous items are \code{NA}.\cr
#' \code{d}        \tab Item intercept matrix used by
#'                       \code{mirt::simdata()}, with columns \code{d0},
#'                       \code{d1}, and \code{d2}. Unused entries are
#'                       \code{NA}.\cr
#' \code{itemtype} \tab Character vector identifying each item as
#'                       \code{"3PL"}, \code{"2PL"}, or \code{"gpcm"}.\cr
#' \code{parTab}   \tab Parameter table containing item types, slopes,
#'                       difficulties, guessing parameters, GPCM step
#'                       difficulties, and mirt intercepts.\cr
#' }
#'
#' @details
#' The latent trait is generated from a linear predictor and additive noise:
#' \itemize{
#'   \item \eqn{\theta = F \beta + E \nu + e}, where \eqn{e} is scaled to
#'   match a target SNR.
#' }
#'
#' The default dataset contains 15 items: eight 3PL items and seven
#' three-category GPCM items. Item slopes are generated as
#' \eqn{a_j \sim \mathrm{Lognormal}(0, 0.25)}.
#'
#' For a 3PL item, the probability of a correct response is
#' \deqn{P(Y_{ij}=1 \mid \theta_i) = c_j + (1-c_j)
#'       \mathrm{logit}^{-1}\{a_j(\theta_i-b_j)\},}
#' and the mirt intercept is computed as \eqn{d_j=-a_jb_j}.
#'
#' Each GPCM item has three response categories, 0, 1, and 2, with two step
#' difficulties, \eqn{b_{1j}} and \eqn{b_{2j}}. These are converted to the
#' mirt intercept parameterization as
#' \deqn{d_{0j}=0,\quad d_{1j}=-a_jb_{1j},\quad
#'       d_{2j}=-a_j(b_{1j}+b_{2j}).}
#'
#' Responses are simulated using the item-specific values in
#' \code{itemtype}, with the 3PL lower-asymptote parameters supplied through
#' the \code{guess} argument to \code{mirt::simdata()}.
#'
#' @usage data(sim_a2)
#'
"sim_a2"

