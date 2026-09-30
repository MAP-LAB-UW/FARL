#' Jaakkola--Jordan variational weight
#'
#' @param x Numeric scalar or vector.
#' @return A numeric vector with the same length as `x`.
#' @noRd
quaeta <- function(x) {
  x <- as.numeric(x)
  out <- numeric(length(x))
  near_zero <- abs(x) < sqrt(.Machine$double.eps)
  out[near_zero] <- 1 / 8
  out[!near_zero] <- tanh(x[!near_zero] / 2) / (4 * x[!near_zero])
  out
}


#' Extract multidimensional binary-item parameters
#'
#' Converts the public `parTab` representation into the loading matrix and
#' intercept vector used internally by MFARLR-GVEM.
#'
#' @param parTab Item parameter data frame.
#' @param J Number of items.
#' @param allow_guessing Whether nonzero 3PL guessing parameters are allowed.
#' @return A list containing `a`, `d`, `c`, `dimension`, `dimension_names`, and
#'   `item_index`.
#' @noRd
multidim_item_parameters <- function(parTab, J, allow_guessing = FALSE) {
  if (!is.data.frame(parTab)) {
    stop("parTab must be a data frame.", call. = FALSE)
  }
  if (nrow(parTab) != J) {
    stop("parTab must contain exactly one row per column of Y.", call. = FALSE)
  }

  get_column <- function(candidates, required = FALSE) {
    match_index <- match(candidates, names(parTab), nomatch = 0L)
    match_index <- match_index[match_index > 0L]
    if (length(match_index)) {
      return(parTab[[match_index[1L]]])
    }
    if (isTRUE(required)) {
      stop(
        "parTab must contain one of: ",
        paste(candidates, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    NULL
  }

  # A full loading matrix may be supplied as a matrix-valued `a` column or as
  # separate a1, ..., aD columns. Otherwise use one slope and one dimension
  # label per item, which is the representation used by sim_m1$parTab.
  a_value <- get_column(c("a"), required = FALSE)
  loading_columns <- grep("^(a|slope)[0-9]+$", names(parTab), value = TRUE)

  if (is.matrix(a_value)) {
    a_matrix <- as.matrix(a_value)
    if (nrow(a_matrix) != J) {
      stop("The matrix-valued parTab$a must have J rows.", call. = FALSE)
    }
    dimension_names <- colnames(a_matrix)
  } else if (length(loading_columns)) {
    a_matrix <- as.matrix(parTab[, loading_columns, drop = FALSE])
    dimension_names <- loading_columns
  } else {
    slope <- get_column(c("a", "slope"), required = TRUE)
    slope <- as.numeric(slope)
    if (length(slope) != J || any(!is.finite(slope))) {
      stop("parTab$a or parTab$slope must contain J finite values.", call. = FALSE)
    }

    dimension_value <- get_column(
      c("dimension", "Dimension", "dimension_index", "subtest"),
      required = TRUE
    )
    if (length(dimension_value) != J || anyNA(dimension_value)) {
      stop("The item dimension column must contain J nonmissing values.", call. = FALSE)
    }

    if (is.numeric(dimension_value) || is.integer(dimension_value)) {
      dimension <- as.integer(dimension_value)
      if (any(dimension < 1L)) {
        stop("Numeric item dimensions must be positive integers.", call. = FALSE)
      }
      D <- max(dimension)
      if (!identical(sort(unique(dimension)), seq_len(D))) {
        stop("Numeric item dimensions must be consecutive from 1 to D.", call. = FALSE)
      }
      dimension_names <- paste0("Dimension", seq_len(D))
    } else {
      dimension_character <- as.character(dimension_value)
      dimension_names <- unique(dimension_character)
      dimension <- match(dimension_character, dimension_names)
      D <- length(dimension_names)
    }

    a_matrix <- matrix(
      0,
      nrow = J,
      ncol = D,
      dimnames = list(NULL, dimension_names)
    )
    a_matrix[cbind(seq_len(J), dimension)] <- slope
  }

  storage.mode(a_matrix) <- "numeric"
  if (!ncol(a_matrix) || any(!is.finite(a_matrix))) {
    stop("The item loading matrix constructed from parTab is invalid.", call. = FALSE)
  }
  if (any(rowSums(a_matrix != 0) == 0L)) {
    stop("Every item must load on at least one dimension.", call. = FALSE)
  }

  D <- ncol(a_matrix)
  if (is.null(dimension_names) || length(dimension_names) != D ||
      anyNA(dimension_names) || any(dimension_names == "")) {
    dimension_names <- paste0("Dimension", seq_len(D))
  }
  colnames(a_matrix) <- make.unique(as.character(dimension_names))

  d_value <- get_column(c("d", "intercept"), required = FALSE)
  if (is.null(d_value)) {
    difficulty <- get_column(c("b", "difficulty"), required = TRUE)
    difficulty <- as.numeric(difficulty)
    if (length(difficulty) != J || any(!is.finite(difficulty))) {
      stop("parTab$b or parTab$difficulty must contain J finite values.", call. = FALSE)
    }
    if (any(rowSums(a_matrix != 0) != 1L)) {
      stop(
        "parTab must provide d/intercept for cross-loading items; difficulty ",
        "can be converted to d only for simple-structure items.",
        call. = FALSE
      )
    }
    active_slope <- rowSums(a_matrix)
    d_value <- -active_slope * difficulty
  }
  d_value <- as.numeric(d_value)
  if (length(d_value) != J || any(!is.finite(d_value))) {
    stop("parTab$d must contain one finite intercept per item.", call. = FALSE)
  }

  c_value <- get_column(c("c", "guessing"), required = FALSE)
  if (is.null(c_value)) {
    c_value <- rep(0, J)
  }
  c_value <- rep_len(as.numeric(c_value), J)
  if (any(!is.finite(c_value)) || any(c_value < 0 | c_value >= 1)) {
    stop("parTab$c/guessing must satisfy 0 <= c < 1.", call. = FALSE)
  }
  if (!isTRUE(allow_guessing) && any(c_value != 0)) {
    stop(
      "MFARLR_mml currently supports multidimensional 2PL items only; ",
      "parTab$c/guessing must be zero.",
      call. = FALSE
    )
  }

  item_index <- lapply(
    seq_len(D),
    function(dimension) which(a_matrix[, dimension] != 0)
  )
  if (any(lengths(item_index) == 0L)) {
    stop("Every latent dimension must have at least one item.", call. = FALSE)
  }
  names(item_index) <- colnames(a_matrix)

  primary_dimension <- max.col(abs(a_matrix), ties.method = "first")

  list(
    a = a_matrix,
    d = d_value,
    c = c_value,
    dimension = primary_dimension,
    dimension_names = colnames(a_matrix),
    item_index = item_index
  )
}


#' Perform one multidimensional GVEM iteration
#'
#' The item discrimination and intercept parameters are treated as fixed. The
#' function updates the subject-level variational moments, the selected latent
#' regression coefficients, the residual covariance matrix, and `xi`.
#'
#' @param resp `N` by `J` binary response matrix. Missing responses may be `NA`.
#' @param Z `N` by `P` latent-regression design matrix.
#' @param sigma_hat Current `D` by `D` residual covariance matrix.
#' @param a `J` by `D` item discrimination matrix.
#' @param d0 Length-`J` item intercept vector.
#' @param xi_hat Current `N` by `J` variational-parameter matrix.
#' @param D Number of latent dimensions.
#' @param index `P` by `D` binary matrix. A value of one indicates that the
#'   corresponding regression coefficient is estimated.
#' @param beta_hat Current `P` by `D` regression-coefficient matrix.
#' @param ridge Nonnegative ridge constant used in the selected regressions.
#'
#' @return A list containing `omega`, `mu`, `sigma`, `xi`, and `coef`.
#' @noRd
gvem_ite_once <- function(
    resp,
    Z,
    sigma_hat,
    a,
    d0,
    xi_hat,
    D,
    index,
    beta_hat,
    ridge = 1e-6
) {
  resp <- as.matrix(resp)
  Z <- as.matrix(Z)
  sigma_hat <- as.matrix(sigma_hat)
  a <- as.matrix(a)
  xi_hat <- as.matrix(xi_hat)
  index <- as.matrix(index)
  beta_hat <- as.matrix(beta_hat)
  d0 <- as.numeric(d0)

  N <- nrow(resp)
  J <- ncol(resp)
  P <- ncol(Z)

  if (!identical(dim(a), c(J, D))) {
    stop("a must be a J by D matrix.", call. = FALSE)
  }
  if (length(d0) != J || any(!is.finite(d0))) {
    stop("d0 must contain J finite item intercepts.", call. = FALSE)
  }
  if (!identical(dim(sigma_hat), c(D, D))) {
    stop("sigma_hat must be a D by D matrix.", call. = FALSE)
  }
  if (!identical(dim(xi_hat), c(N, J))) {
    stop("xi_hat must be an N by J matrix.", call. = FALSE)
  }
  if (!identical(dim(index), c(P, D))) {
    stop("index must be a ncol(Z) by D matrix.", call. = FALSE)
  }
  if (!identical(dim(beta_hat), c(P, D))) {
    stop("beta_hat must be a ncol(Z) by D matrix.", call. = FALSE)
  }
  if (any(!is.finite(Z)) || any(!is.finite(a)) ||
      any(!is.finite(sigma_hat)) || any(!is.finite(xi_hat))) {
    stop("Z, a, sigma_hat, and xi_hat must be finite.", call. = FALSE)
  }
  observed_response <- resp[!is.na(resp)]
  if (any(!observed_response %in% c(0, 1))) {
    stop("resp must contain only 0, 1, or NA.", call. = FALSE)
  }
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be one nonnegative finite number.", call. = FALSE)
  }

  sigma_inverse <- tryCatch(
    solve(sigma_hat),
    error = function(e) {
      stop("sigma_hat is singular: ", conditionMessage(e), call. = FALSE)
    }
  )

  prior_mean <- Z %*% beta_hat
  mu_hat <- matrix(0, nrow = N, ncol = D)
  omega_hat <- array(0, dim = c(N, D, D))
  hxi <- matrix(quaeta(xi_hat), nrow = N, ncol = J)

  # E-step: update the variational posterior moments.
  for (i in seq_len(N)) {
    posterior_precision <- sigma_inverse
    linear_term <- numeric(D)
    observed_items <- which(!is.na(resp[i, ]))

    for (j in observed_items) {
      loading <- a[j, ]
      weight <- hxi[i, j]
      posterior_precision <-
        posterior_precision + 2 * weight * tcrossprod(loading)
      linear_term <-
        linear_term +
        (resp[i, j] - 2 * weight * d0[j] - 0.5) * loading
    }

    omega_i <- tryCatch(
      solve(posterior_precision),
      error = function(e) {
        stop(
          "Posterior precision is singular for respondent ", i, ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    omega_hat[i, , ] <- omega_i
    mu_hat[i, ] <- drop(
      omega_i %*% (sigma_inverse %*% prior_mean[i, ] + linear_term)
    )
  }

  # Update the local variational parameters.
  for (i in seq_len(N)) {
    observed_items <- which(!is.na(resp[i, ]))
    for (j in observed_items) {
      loading <- a[j, ]
      posterior_variance <- drop(
        crossprod(loading, omega_hat[i, , ] %*% loading)
      )
      posterior_location <- sum(loading * mu_hat[i, ]) + d0[j]
      xi_hat[i, j] <- sqrt(max(
        .Machine$double.eps,
        posterior_location^2 + posterior_variance
      ))
    }
  }

  # M-step: update only coefficients selected by the dimension-specific
  # FARLR_Debias fits.
  beta_new <- matrix(
    0,
    nrow = P,
    ncol = D,
    dimnames = dimnames(beta_hat)
  )

  for (dimension in seq_len(D)) {
    selected <- which(index[, dimension] != 0)
    if (length(selected) == 0L) {
      next
    }

    Z_selected <- Z[, selected, drop = FALSE]
    cross_product <- crossprod(Z_selected)
    diag(cross_product) <- diag(cross_product) + ridge

    beta_new[selected, dimension] <- drop(
      solve(
        cross_product,
        crossprod(Z_selected, mu_hat[, dimension])
      )
    )
  }

  # Update the residual covariance using the posterior covariance and the
  # posterior residual outer products.
  fitted_mean <- Z %*% beta_new
  sigma_new <- matrix(0, nrow = D, ncol = D)

  for (i in seq_len(N)) {
    residual_i <- mu_hat[i, ] - fitted_mean[i, ]
    sigma_new <-
      sigma_new + omega_hat[i, , ] + tcrossprod(residual_i)
  }

  sigma_new <- sigma_new / N
  sigma_new <- (sigma_new + t(sigma_new)) / 2

  list(
    omega = omega_hat,
    mu = mu_hat,
    sigma = sigma_new,
    xi = xi_hat,
    coef = beta_new
  )
}

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
