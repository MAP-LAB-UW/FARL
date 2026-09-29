#' Fit a high-dimensional multidimensional FARLR model
#'
#' Each latent dimension is first screened separately using the
#' `FARLR_Debias` method in `Farlr_mml()`. The union of the dimension-specific
#' nonzero patterns is stored in a `P` by `D` index matrix. Conditional on this
#' fixed sparsity pattern, multidimensional GVEM iterations jointly update the
#' regression coefficients and latent residual covariance matrix.
#'
#' @param X `N` by `p` covariate matrix.
#' @param Y `N` by `J` binary response matrix.
#' @param parTab Item parameter data frame with one row per item. For the
#'   simple-structure format used by `sim_m1`, it must contain `a` (or `slope`),
#'   `d` (or `b`/`difficulty`), and `dimension` (or `subtest`). A full loading
#'   matrix can alternatively be supplied through `a1`, ..., `aD` columns.
#'   The current GVEM implementation supports multidimensional 2PL items, so
#'   `c`/`guessing`, when present, must be zero.
#' @param farlr_args Named list of additional arguments passed to each
#'   `Farlr_mml()` call. Do not include `X`, `Y`, `parTab`, or `method`.
#' @param selection_tol Coefficients with absolute value greater than this
#'   value are selected in `index`.
#' @param sigma_start Optional positive-definite `D` by `D` starting covariance.
#' @param xi_start Optional positive `N` by `J` starting matrix for `xi`.
#' @param max_iter Maximum number of joint GVEM iterations.
#' @param threshold Convergence threshold for the maximum absolute parameter
#'   change in `sigma_hat` and `beta_hat`.
#' @param ridge Nonnegative ridge constant used in the selected regressions.
#' @param progress Logical; display a progress bar for the joint GVEM
#'   iterations when `TRUE`.
#' @param verbose Logical; report dimension-specific screening results and the
#'   final convergence status when `TRUE`. Iteration differences are printed
#'   when `verbose = TRUE` and `progress = FALSE`.
#'
#' @return A list containing the dimension-specific FARLR fits, selection
#'   `index`, joint coefficient matrix, residual covariance, variational
#'   moments, convergence information, the design matrix `Z`, and the original
#'   `X`, `Y`, and `parTab` inputs required for plausible-value generation.
#'
#' @examples
#' \dontrun{
#' fit <- MFARLR_mml(
#'   X = sim_m1$X,
#'   Y = sim_m1$Y,
#'   parTab = sim_m1$parTab,
#'   farlr_args = list(main = sim_m1$main)
#' )
#' }
#'
#' @export
MFARLR_mml<- function(
    X,
    Y,
    parTab,
    farlr_args = list(),
    selection_tol = 0,
    sigma_start = NULL,
    xi_start = NULL,
    max_iter = 1000L,
    threshold = 1e-4,
    ridge = 1e-6,
    progress = TRUE,
    verbose = TRUE
) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.", call. = FALSE)
  }
  suppressWarnings(storage.mode(Y) <- "numeric")
  if (nrow(Y) != nrow(X)) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }

  N <- nrow(Y)
  J <- ncol(Y)
  item_parameters <- multidim_item_parameters(
    parTab = parTab,
    J = J,
    allow_guessing = FALSE
  )
  a <- item_parameters$a
  d <- item_parameters$d
  item_index <- item_parameters$item_index
  D <- ncol(a)
  dimension_names <- item_parameters$dimension_names

  if (nrow(a) != J || D < 1L || any(!is.finite(a))) {
    stop("a must be a finite J by D matrix.", call. = FALSE)
  }
  if (length(d) != J || any(!is.finite(d))) {
    stop("d must contain one finite intercept per item.", call. = FALSE)
  }
  observed_response <- Y[!is.na(Y)]
  if (any(!observed_response %in% c(0, 1))) {
    stop(
      "The current multidimensional GVEM implementation supports only ",
      "binary responses coded 0/1 (with optional NA).",
      call. = FALSE
    )
  }
  if (!is.list(farlr_args) || is.null(names(farlr_args)) && length(farlr_args)) {
    stop("farlr_args must be a named list.", call. = FALSE)
  }

  protected_arguments <- c("X", "Y", "parTab", "method")
  duplicated_arguments <- intersect(names(farlr_args), protected_arguments)
  if (length(duplicated_arguments) > 0L) {
    stop(
      "Do not include these arguments in farlr_args: ",
      paste(duplicated_arguments, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  max_iter <- as.integer(max_iter)
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    stop("max_iter must be a positive integer.", call. = FALSE)
  }
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("threshold must be one positive finite number.", call. = FALSE)
  }
  if (length(selection_tol) != 1L ||
      !is.finite(selection_tol) || selection_tol < 0) {
    stop("selection_tol must be one nonnegative finite number.", call. = FALSE)
  }
  if (length(progress) != 1L || is.na(progress) || !is.logical(progress)) {
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(verbose) != 1L || is.na(verbose) || !is.logical(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  if (!exists("Farlr_mml", mode = "function")) {
    stop("Farlr_mml() was not found.", call. = FALSE)
  }

  dimension_fits <- vector("list", D)
  names(dimension_fits) <- dimension_names
  Z <- NULL
  coefficient_names <- NULL
  coefficient_matrix <- NULL
  index <- NULL

  gvem_progress <- NULL
  on.exit({
    if (!is.null(gvem_progress)) {
      close(gvem_progress)
    }
  }, add = TRUE)

  # Preliminary dimension-by-dimension FARLR_Debias screening.
  for (dimension in seq_len(D)) {
    items <- item_index[[dimension]]
    parTab_dimension <- data.frame(
      a = a[items, dimension],
      d = d[items],
      c = 0
    )

    fit_arguments <- c(
      list(
        X = X,
        Y = Y[, items, drop = FALSE],
        parTab = parTab_dimension,
        method = "FARLR_Debias"
      ),
      farlr_args
    )

    fit_dimension <- do.call(Farlr_mml, fit_arguments)
    dimension_fits[[dimension]] <- fit_dimension

    Z_dimension <- fit_dimension$Fan
    if (!is.matrix(Z_dimension)) {
      factor_scores <- as.matrix(fit_dimension$factor_scores)
      residual_covariates <- as.matrix(fit_dimension$hatU)
      Z_dimension <- cbind(factor_scores, residual_covariates)
    }

    coefficients_dimension <- as.numeric(fit_dimension$coefficients)
    names_dimension <- names(fit_dimension$coefficients)
    if (is.null(names_dimension)) {
      names_dimension <- colnames(Z_dimension)
    }

    if (ncol(Z_dimension) != length(coefficients_dimension)) {
      stop(
        "Dimension ", dimension,
        " returned incompatible design and coefficient dimensions.",
        call. = FALSE
      )
    }
    if (is.null(names_dimension) || length(names_dimension) != ncol(Z_dimension)) {
      stop(
        "Dimension ", dimension,
        " does not provide valid coefficient names.",
        call. = FALSE
      )
    }
    colnames(Z_dimension) <- names_dimension

    if (dimension == 1L) {
      Z <- Z_dimension
      coefficient_names <- names_dimension
      coefficient_matrix <- matrix(
        0,
        nrow = ncol(Z),
        ncol = D,
        dimnames = list(coefficient_names, dimension_names)
      )
      index <- matrix(
        0L,
        nrow = ncol(Z),
        ncol = D,
        dimnames = dimnames(coefficient_matrix)
      )
    } else {
      name_order <- match(coefficient_names, names_dimension)
      if (anyNA(name_order)) {
        stop(
          "Coefficient names are inconsistent across dimensions.",
          call. = FALSE
        )
      }
      Z_dimension <- Z_dimension[, name_order, drop = FALSE]
      coefficients_dimension <- coefficients_dimension[name_order]

      if (!isTRUE(all.equal(Z_dimension, Z, tolerance = 1e-8))) {
        stop(
          "The FARLR design matrices differ across dimensions. ",
          "Use the same X, K_hat, and factor-analysis settings in every fit.",
          call. = FALSE
        )
      }
    }

    selected <- which(abs(coefficients_dimension) > selection_tol)
    index[selected, dimension] <- 1L
    coefficient_matrix[selected, dimension] <-
      coefficients_dimension[selected]

    if (isTRUE(verbose)) {
      message(
        "Dimension ", dimension, ": selected ", length(selected),
        " of ", length(coefficients_dimension), " coefficients."
      )
    }
  }

  beta_hat <- coefficient_matrix

  if (is.null(sigma_start)) {
    sigma_hat <- diag(D)
  } else {
    sigma_hat <- as.matrix(sigma_start)
    if (!identical(dim(sigma_hat), c(D, D)) ||
        any(!is.finite(sigma_hat))) {
      stop("sigma_start must be a finite D by D matrix.", call. = FALSE)
    }
    if (any(eigen(sigma_hat, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
      stop("sigma_start must be positive definite.", call. = FALSE)
    }
  }

  if (is.null(xi_start)) {
    xi_hat <- matrix(1, nrow = N, ncol = J)
  } else {
    xi_hat <- as.matrix(xi_start)
    if (!identical(dim(xi_hat), c(N, J)) ||
        any(!is.finite(xi_hat)) || any(xi_hat <= 0)) {
      stop("xi_start must be a positive finite N by J matrix.", call. = FALSE)
    }
  }

  old_parameter <- c(sigma_hat, beta_hat)
  converged <- FALSE
  difference <- Inf
  mu_hat <- matrix(NA_real_, nrow = N, ncol = D)
  omega_hat <- array(NA_real_, dim = c(N, D, D))

  if (isTRUE(progress)) {
    message("Running joint GVEM iterations")
    gvem_progress <- utils::txtProgressBar(
      min = 0L,
      max = max_iter,
      initial = 0L,
      style = 3L
    )
  }

  for (iteration in seq_len(max_iter)) {
    iteration_result <- gvem_ite_once(
      resp = Y,
      Z = Z,
      sigma_hat = sigma_hat,
      a = a,
      d0 = d,
      xi_hat = xi_hat,
      D = D,
      index = index,
      beta_hat = beta_hat,
      ridge = ridge
    )

    mu_hat <- iteration_result$mu
    sigma_hat <- iteration_result$sigma
    xi_hat <- iteration_result$xi
    omega_hat <- iteration_result$omega
    beta_hat <- iteration_result$coef

    new_parameter <- c(sigma_hat, beta_hat)
    difference <- max(abs(old_parameter - new_parameter))

    if (isTRUE(progress)) {
      utils::setTxtProgressBar(gvem_progress, iteration)
    } else if (isTRUE(verbose)) {
      message("Iteration ", iteration, ": diff = ", signif(difference, 6))
    }

    if (difference < threshold) {
      converged <- TRUE
      if (isTRUE(progress)) {
        utils::setTxtProgressBar(gvem_progress, max_iter)
      }
      break
    }

    old_parameter <- new_parameter
  }

  if (!is.null(gvem_progress)) {
    close(gvem_progress)
    gvem_progress <- NULL
  }

  if (isTRUE(verbose)) {
    if (converged) {
      message("Converged at iteration ", iteration, ".")
    } else {
      message(
        "Did not converge after ", max_iter,
        " iterations; final diff = ", signif(difference, 6), "."
      )
    }
  }

  result <- list(
    coefficients = beta_hat,
    sigma = sigma_hat,
    mu = mu_hat,
    omega = omega_hat,
    xi = xi_hat,
    index = index,
    Z = Z,
    X = X,
    Y = Y,
    parTab = parTab,
    response = Y,
    a = a,
    d = d,
    dimension_fits = dimension_fits,
    item_index = item_index,
    iterations = iteration,
    difference = difference,
    convergence = if (converged) "Converged" else "Did not converge",
    call = match.call()
  )

  class(result) <- c("MFARLR_GVEM", "list")
  result
}
