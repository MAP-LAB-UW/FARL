#' FARL: Use Factor-Augmented Regularized Latent Regression for Large-Scale Assessment
#'
#' FARL supports latent regression analyses for large-scale assessments (LSAs),
#' where student background information is incorporated into the population
#' model used to generate plausible values (PVs). These background variables
#' are often high-dimensional and strongly dependent, which can make
#' conventional latent regression unstable and complicate variable selection.
#'
#' FARL implements factor-augmented regularized latent regression (FARLR).
#' FARLR decomposes covariate variation into common latent factors and
#' idiosyncratic components, includes both components in the latent regression,
#' and uses sparse regularization to select a parsimonious set of relevant
#' covariates. This structure is designed to stabilize estimation while
#' retaining interpretable background-variable effects.
#'
#' @section Model estimation:
#' \itemize{
#'   \item \code{\link{Farlr_mml}} is the unified interface for fitting FARLR
#'   models. It supports two estimation methods:
#'   \itemize{
#'     \item \code{method = "FARLR_EMM"}: importance-sampling-based
#'     expectation-maximization-maximization estimation, with a regularized
#'     variable-selection step followed by an unpenalized coefficient update.
#'     \item \code{method = "FARLR_Debias"}: regularized estimation using the
#'     factor and idiosyncratic components, followed by correction of LASSO
#'     shrinkage bias.
#'   }
#'   \item \code{\link{Dire_mml}} fits a PCA-based DIRE latent regression model
#'   for comparison with the FARLR approaches.
#'   \item \code{\link{MFARLR_mml}} fits multidimensional FARLR by combining
#'   dimension-specific sparse covariate screening with joint Gaussian
#'   variational expectation-maximization estimation.
#'   \item \code{\link{MDIRE_mml}} fits the multidimensional PCA-based DIRE
#'   comparison model.
#' }
#'
#' @section Plausible-value generation:
#' \itemize{
#'   \item \code{\link{Farlr_drawPVs}} generates plausible values from the
#'   posterior latent-trait distributions implied by a fitted FARLR model.
#'   \item \code{\link{Dire_drawPVs}} generates plausible values from a fitted
#'   DIRE model.
#'   \item \code{\link{MFARLR_drawPVs}} generates correlated multidimensional
#'   plausible values from a fitted MFARLR model.
#'   \item \code{\link{MDIRE_drawPVs}} generates plausible values from a fitted
#'   multidimensional DIRE model.
#' }
#'
#' @section Simulation and example data:
#' \itemize{
#'   \item \code{\link{sim_a1}} contains simulated dichotomous 2PL responses,
#'   item parameters, and high-dimensional background covariates.
#'   \item \code{\link{sim_a2}} contains simulated mixed 3PL and GPCM responses,
#'   item parameters, and high-dimensional background covariates.
#'   \item \code{\link{sim_m1}} contains simulated multidimensional item
#'   responses, item parameters, latent traits, and background covariates.
#' }
#'
#' @section Methodological scope:
#' FARLR is intended to improve compatibility between plausible-value
#' generation and common secondary analyses, but it does not guarantee
#' congeniality in the strict statistical sense. The package supports both
#' unidimensional and multidimensional latent proficiency models; currently,
#' \code{MFARLR_mml} supports multidimensional binary 2PL items.
#'
#' @keywords internal
"_PACKAGE"
