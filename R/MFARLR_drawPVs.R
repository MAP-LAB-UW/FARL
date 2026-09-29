#' Draw Multidimensional Plausible Values from an MFARLR-GVEM Model
#'
#' Draws Monte Carlo samples from each respondent's fitted multidimensional
#' latent-regression prior, weights the samples by the multidimensional 2PL
#' item-response likelihood, and computes respondent-specific posterior EAP
#' means and covariance matrices. Plausible values are then drawn from the
#' multivariate distribution determined by these posterior moments.
#'
#' @param mmlcomp A fitted object returned by [MFARLR_GVEM()]. The object
#'   retains `X`, `Y`, `parTab`, the multidimensional loading matrix, and the
#'   item intercepts, so these inputs do not need to be supplied again.
#' @param npv Positive integer. Number of plausible values per respondent.
#' @param n_mc Positive integer. Number of prior Monte Carlo samples used for
#'   each respondent. The original implementation used `20000`.
#' @param theta Optional true `N` by `D` latent-trait matrix used only to report
#'   the difference between the mean EAP estimate and the mean true score.
#' @param seed Optional random-number seed.
#' @param covariance_floor Small positive eigenvalue used to stabilize posterior
#'   covariance matrices before drawing plausible values.
#' @param progress Logical; display a respondent-level progress bar.
#' @param verbose Logical; print a completion message and, when `theta` is
#'   supplied, the EAP mean differences.
#' @param return_details Logical. If `TRUE`, return posterior summaries along
#'   with the plausible values. If `FALSE`, return only the wide PV data frame.
#'
#' @return When `return_details = TRUE`, a list containing:
#'   \describe{
#'     \item{data}{An `N` by `npv` by `D` plausible-value array.}
#'     \item{datPVs}{A wide data frame with one row per respondent.}
#'     \item{rawPV}{An `(N * npv)` by `D` stacked plausible-value matrix.}
#'     \item{EAP_estimates}{An `N` by `D` matrix of posterior EAP means.}
#'     \item{EAP_covariance}{An `N` by `D` by `D` array of posterior
#'       covariance matrices.}
#'     \item{effective_sample_size}{The importance-sampling effective sample
#'       size for every respondent.}
#'     \item{prior_mean}{The fitted latent-regression prior means.}
#'     \item{sigma}{The fitted latent residual covariance matrix.}
#'     \item{mean_difference}{Dimension-specific EAP minus true-score mean
#'       differences when `theta` is supplied.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- MFARLR_GVEM(
#'   X = sim_m1$X,
#'   Y = sim_m1$Y,
#'   parTab = sim_m1$parTab,
#'   farlr_args = list(main = sim_m1$main)
#' )
#'
#' PVs <- MFARLR_drawPVs(
#'   mmlcomp = fit,
#'   npv = 10L
#' )
#' }
#'
#' @export
MFARLR_drawPVs <- function(
    mmlcomp,
    npv = 10L,
    n_mc = 20000L,
    theta = NULL,
    seed = NULL,
    covariance_floor = 1e-8,
    progress = TRUE,
    verbose = TRUE,
    return_details = TRUE
) {
  if (!is.list(mmlcomp)) {
    stop("mmlcomp must be an object returned by MFARLR_GVEM().", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  npv <- as.integer(npv)
  n_mc <- as.integer(n_mc)

  if (length(npv) != 1L || is.na(npv) || npv < 1L) {
    stop("npv must be a positive integer.", call. = FALSE)
  }
  if (length(n_mc) != 1L || is.na(n_mc) || n_mc < 2L) {
    stop("n_mc must be an integer of at least 2.", call. = FALSE)
  }
  if (length(covariance_floor) != 1L ||
      !is.finite(covariance_floor) || covariance_floor <= 0) {
    stop("covariance_floor must be one positive finite number.", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.", call. = FALSE)
  }

  get_component <- function(candidates) {
    for (name in candidates) {
      value <- mmlcomp[[name]]
      if (!is.null(value)) {
        return(value)
      }
    }
    stop(
      "mmlcomp does not contain any of: ",
      paste(candidates, collapse = ", "),
      ". Refit the model with the current MFARLR_GVEM().",
      call. = FALSE
    )
  }

  response <- as.matrix(get_component(c("Y", "response")))
  a <- as.matrix(get_component("a"))
  d <- as.numeric(get_component("d"))
  Z <- as.matrix(get_component("Z"))
  beta_hat <- as.matrix(get_component("coefficients"))
  sigma_hat <- as.matrix(get_component("sigma"))

  suppressWarnings(storage.mode(response) <- "numeric")
  storage.mode(Z) <- "numeric"
  storage.mode(beta_hat) <- "numeric"
  storage.mode(a) <- "numeric"
  storage.mode(sigma_hat) <- "numeric"

  N <- nrow(response)
  J <- ncol(response)
  D <- ncol(a)

  if (nrow(Z) != N) {
    stop("mmlcomp$Y and mmlcomp$Z must have the same number of rows.", call. = FALSE)
  }
  if (!identical(dim(a), c(J, D))) {
    stop("a must be a J by D matrix.", call. = FALSE)
  }
  if (length(d) != J || any(!is.finite(d))) {
    stop("d must contain one finite intercept per item.", call. = FALSE)
  }
  if (!identical(dim(sigma_hat), c(D, D))) {
    stop("mmlcomp$sigma must be a D by D matrix.", call. = FALSE)
  }
  if (ncol(beta_hat) != D || nrow(beta_hat) != ncol(Z)) {
    stop(
      "mmlcomp$coefficients must be a ncol(Z) by D matrix.",
      call. = FALSE
    )
  }
  if (any(!is.finite(Z)) || any(!is.finite(beta_hat)) ||
      any(!is.finite(a)) || any(!is.finite(sigma_hat))) {
    stop("Z, coefficients, a, and sigma must be finite.", call. = FALSE)
  }

  observed_response <- response[!is.na(response)]
  if (any(!observed_response %in% c(0, 1))) {
    stop("response must contain only 0, 1, or NA.", call. = FALSE)
  }

  # Align the coefficient rows with the design columns when names are present.
  if (!is.null(colnames(Z)) && !is.null(rownames(beta_hat))) {
    coefficient_order <- match(colnames(Z), rownames(beta_hat))
    if (anyNA(coefficient_order)) {
      stop("Could not align Z columns with coefficient rows.", call. = FALSE)
    }
    beta_hat <- beta_hat[coefficient_order, , drop = FALSE]
  }

  stabilize_covariance <- function(covariance_matrix) {
    covariance_matrix <-
      (as.matrix(covariance_matrix) + t(as.matrix(covariance_matrix))) / 2
    decomposition <- eigen(covariance_matrix, symmetric = TRUE)
    eigenvalues <- pmax(decomposition$values, covariance_floor)
    stabilized <-
      decomposition$vectors %*%
      (eigenvalues * t(decomposition$vectors))
    (stabilized + t(stabilized)) / 2
  }

  sigma_hat <- stabilize_covariance(sigma_hat)
  prior_mean <- Z %*% beta_hat
  dimension_names <- colnames(beta_hat)
  if (is.null(dimension_names) || length(dimension_names) != D) {
    dimension_names <- paste0("Dimension", seq_len(D))
  }
  colnames(prior_mean) <- dimension_names

  EAP_estimates <- matrix(
    0,
    nrow = N,
    ncol = D,
    dimnames = list(rownames(response), colnames(prior_mean))
  )
  EAP_covariance <- array(
    0,
    dim = c(N, D, D),
    dimnames = list(
      rownames(response),
      colnames(prior_mean),
      colnames(prior_mean)
    )
  )
  effective_sample_size <- numeric(N)
  plausible_values <- array(
    NA_real_,
    dim = c(N, npv, D),
    dimnames = list(
      rownames(response),
      paste0("PV", seq_len(npv)),
      colnames(prior_mean)
    )
  )

  progress_bar <- NULL
  if (isTRUE(progress)) {
    progress_bar <- utils::txtProgressBar(
      min = 0,
      max = N,
      initial = 0,
      style = 3
    )
  } else if (isTRUE(verbose)) {
    message("Drawing multidimensional plausible values...")
  }

  close_progress_bar <- function() {
    if (inherits(progress_bar, "txtProgressBar")) {
      close(progress_bar)
      progress_bar <<- NULL
    }
    invisible(NULL)
  }
  on.exit(close_progress_bar(), add = TRUE)

  for (i in seq_len(N)) {
    theta_samples <- MASS::mvrnorm(
      n = n_mc,
      mu = prior_mean[i, ],
      Sigma = sigma_hat
    )
    theta_samples <- matrix(theta_samples, ncol = D)

    observed_items <- which(!is.na(response[i, ]))

    if (length(observed_items) == 0L) {
      log_likelihood <- numeric(n_mc)
    } else {
      a_observed <- a[observed_items, , drop = FALSE]
      d_observed <- d[observed_items]
      y_observed <- as.numeric(response[i, observed_items])

      linear_predictor <- theta_samples %*% t(a_observed)
      linear_predictor <- sweep(linear_predictor, 2L, d_observed, "+")
      probability <- stats::plogis(linear_predictor)
      probability <- pmin(pmax(probability, 1e-12), 1 - 1e-12)

      log_likelihood <- drop(
        log(probability) %*% y_observed +
          log1p(-probability) %*% (1 - y_observed)
      )
    }

    log_likelihood <- log_likelihood - max(log_likelihood)
    weights <- exp(log_likelihood)
    weight_sum <- sum(weights)

    if (!is.finite(weight_sum) || weight_sum <= 0) {
      stop(
        "Non-finite posterior weights for respondent ", i, ".",
        call. = FALSE
      )
    }

    weights <- weights / weight_sum
    effective_sample_size[i] <- 1 / sum(weights^2)

    EAP_estimates[i, ] <- colSums(theta_samples * weights)
    centered_samples <- sweep(
      theta_samples,
      2L,
      EAP_estimates[i, ],
      "-"
    )
    posterior_covariance <-
      crossprod(centered_samples, centered_samples * weights)
    posterior_covariance <- stabilize_covariance(posterior_covariance)
    EAP_covariance[i, , ] <- posterior_covariance

    plausible_values[i, , ] <- matrix(
      MASS::mvrnorm(
        n = npv,
        mu = EAP_estimates[i, ],
        Sigma = posterior_covariance
      ),
      nrow = npv,
      ncol = D
    )

    if (inherits(progress_bar, "txtProgressBar")) {
      utils::setTxtProgressBar(progress_bar, i)
    }
  }

  close_progress_bar()

  rawPV <- matrix(
    NA_real_,
    nrow = N * npv,
    ncol = D,
    dimnames = list(NULL, colnames(prior_mean))
  )
  for (i in seq_len(N)) {
    rows <- ((i - 1L) * npv + 1L):(i * npv)
    rawPV[rows, ] <- plausible_values[i, , ]
  }

  subject_ids <- rownames(response)
  if (is.null(subject_ids)) {
    subject_ids <- as.character(seq_len(N))
  }

  datPVs <- data.frame(id = subject_ids, check.names = FALSE)
  for (dimension in seq_len(D)) {
    for (pv_index in seq_len(npv)) {
      variable_name <- paste0(
        "_mfarlr_D", dimension,
        "_PV", pv_index
      )
      datPVs[[variable_name]] <- plausible_values[, pv_index, dimension]
    }
  }

  mean_difference <- rep(NA_real_, D)
  names(mean_difference) <- colnames(prior_mean)
  if (!is.null(theta)) {
    theta <- as.matrix(theta)
    if (!identical(dim(theta), c(N, D)) || any(!is.finite(theta))) {
      stop("theta must be a finite N by D matrix.", call. = FALSE)
    }
    mean_difference <- colMeans(EAP_estimates) - colMeans(theta)
    if (isTRUE(verbose)) {
      message(
        "EAP minus theta mean differences: ",
        paste(round(mean_difference, 3), collapse = ", ")
      )
    }
  }

  if (isTRUE(verbose)) {
    message(
      "Median importance-sampling ESS: ",
      round(stats::median(effective_sample_size), 1),
      " of ", n_mc, "."
    )
  }

  if (!isTRUE(return_details)) {
    return(datPVs)
  }

  list(
    data = plausible_values,
    datPVs = datPVs,
    rawPV = rawPV,
    EAP_estimates = EAP_estimates,
    EAP_covariance = EAP_covariance,
    effective_sample_size = effective_sample_size,
    prior_mean = prior_mean,
    sigma = sigma_hat,
    mean_difference = mean_difference,
    npv = npv,
    n_mc = n_mc
  )
}
