

mml_test_1 <- function(){
  mmlcomp_1 <-  with(sim_a1, Farlr_mml(X, Y, parTab, method = "FARLR_Debias", main = c(1,2,15,29,45)))
  invisible(mmlcomp_1)
  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp_1,
    npv = 10
  )
  return(PVs)
}
mml_test_2 <- function(){
  mmlcomp_2 <- with(
    sim_a1,
    Farlr_mml(
      X = X,
      Y = Y,
      parTab = parTab,
      main = c(1,2,15,29,45),
      method = "FARLR_EMM"
    )
  )
  invisible(mmlcomp_2)
  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp_2,
    npv = 10)
  return(PVs)
}
mml_test_3 <- function(){
  mmlcomp_3 <-  with(sim_a2, Farlr_mml(X, Y, parTab, method = "FARLR_Debias", main = c(1,2,15,29,45)))
  #estimate_3 <- cbind(mmlcomp_3$factor_scores, mmlcomp_3$hatU) %*% mmlcomp_3$coef
  #sum((theta-estimate_3)^2)
  invisible(mmlcomp_3)
  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp_3,
    npv = 10
  )
  return(PVs)
}
mml_test_4 <- function(){
  mmlcomp_4 <- with(
    sim_a2,
    Farlr_mml(
      X = X,
      Y = Y,
      parTab = parTab,
      main = c(1,2,15,29,45),
      method = "FARLR_EMM"
    )
  )
  # estimate_4 <- cbind(mmlcomp_4$factor_scores, scale(X)) %*% mmlcomp_4$coef
  # sum((theta-estimate_4)^2)
  invisible(mmlcomp_4)
  PVs <- Farlr_drawPVs(
    mmlcomp = mmlcomp_4,
    npv = 10
  )
}
mml_test_Dire <- function() {

  X <- as.matrix(sim_a1$X)
  resp <- as.matrix(sim_a1$Y)
  n <- nrow(X)
  J <- ncol(resp)

  colnames(X) <- paste0("X", seq_len(ncol(X)))

  main <- c(1, 2, 15, 29, 45)
  main_vars <- paste0("X", main)

  subject <- factor(seq_len(n))
  stuDat <- data.frame(subject = subject, X)

  itemNames <- paste0("i", sprintf("%03d", seq_len(J)))
  colnames(resp) <- itemNames
  resp0 <- data.frame(resp, subject = subject)

  stuItems <- reshape(
    resp0,
    varying = itemNames,
    idvar = "subject",
    direction = "long",
    v.names = "score",
    times = itemNames,
    timevar = "key"
  )

  new_itemNames <- paste0("item", seq_len(J))
  stuItems$key <- rep(new_itemNames, each = n)

  parTab <- sim_a1$parTab[, c(
    "ItemID", "test", "subtest", "slope",
    "difficulty", "guessing", "D"
  )]
  parTab$ItemID <- new_itemNames

  testDat <- data.frame(
    test = "comp",
    subtest = "main",
    location = 0,
    scale = 1
  )

  mmlcomp <- Dire_mml(
    formula = reformulate(main_vars, response = "comp"),
    stuItems = stuItems,
    stuDat = stuDat,
    idVar = "subject",
    dichotParamTab = parTab,
    testScale = testDat,
    X = X,
    main_vars = main_vars,
    pca_type = "cov",
    use_residual = TRUE,
    var_threshold = 0.90
  )
  invisible(mmlcomp)
  PVs <- Dire_drawPVs(
    x = mmlcomp,
    npv = 10L,
    pvVariableNameSuffix = "_dire"
  )$data

  # Return PVs in the original sim_a1 subject order
  pv_order <- match(
    as.character(seq_len(n)),
    as.character(PVs$id)
  )

  PVs <- PVs[pv_order, , drop = FALSE]
  rownames(PVs) <- NULL

  return(PVs)
}
# Test FARLR-Debias, FARLR-EMM, and DIRE using sim_a1.
mml_test_all <- function(
    sim_data = sim_a1,
    theta = NULL,
    main = c(1L, 2L, 15L, 29L, 45L),
    npv = 10L,
    seed = 2026L,
    plot_file = NULL
) {
  X <- as.matrix(sim_data$X)
  response <- as.matrix(sim_data$Y)
  parTab <- sim_data$parTab

  n <- nrow(X)
  J <- ncol(response)

  if (nrow(response) != n) {
    stop("sim_data$X and sim_data$Y must have the same number of rows.")
  }

  if (is.null(theta)) {
    theta <- sim_data$theta
  }
  if (is.null(theta)) {
    stop(
      "Supply `theta`, or save the true latent trait as `sim_data$theta`."
    )
  }

  theta <- as.numeric(theta)
  if (length(theta) != n || any(!is.finite(theta))) {
    stop("theta must contain one finite value per respondent.")
  }

  npv <- as.integer(npv)
  if (length(npv) != 1L || is.na(npv) || npv < 1L) {
    stop("npv must be a positive integer.")
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }

  main <- as.integer(main)
  main_vars <- colnames(X)[main]

  # -------------------------------------------------------------------------
  # 1. FARLR-Debias
  # -------------------------------------------------------------------------
  mmlcomp_debias <- Farlr_mml(
    X = X,
    Y = response,
    parTab = parTab,
    method = "FARLR_Debias",
    main = main,
    seed = seed
  )

  PVs_debias <- Farlr_drawPVs(
    mmlcomp = mmlcomp_debias,
    npv = npv,
    theta = theta,
    seed = seed + 1L,
    progress = FALSE,
    verbose = FALSE
  )

  # -------------------------------------------------------------------------
  # 2. FARLR-EMM
  # -------------------------------------------------------------------------
  mmlcomp_emm <- Farlr_mml(
    X = X,
    Y = response,
    parTab = parTab,
    method = "FARLR_EMM",
    main = main,
    seed = seed + 2L
  )

  PVs_emm <- Farlr_drawPVs(
    mmlcomp = mmlcomp_emm,
    npv = npv,
    theta = theta,
    seed = seed + 3L,
    progress = FALSE,
    verbose = FALSE
  )

  # -------------------------------------------------------------------------
  # 3. DIRE
  # -------------------------------------------------------------------------
  subject <- factor(seq_len(n))
  stuDat <- data.frame(subject = subject, X, check.names = FALSE)

  item_names <- paste0("i", sprintf("%03d", seq_len(J)))
  colnames(response) <- item_names
  resp0 <- data.frame(response, subject = subject, check.names = FALSE)

  stuItems <- reshape(
    data = resp0,
    varying = item_names,
    idvar = "subject",
    direction = "long",
    v.names = "score",
    times = item_names,
    timevar = "key"
  )

  new_item_names <- paste0("item", seq_len(J))
  stuItems$key <- rep(new_item_names, each = n)

  dire_parTab <- parTab[, c(
    "ItemID", "test", "subtest", "slope",
    "difficulty", "guessing", "D"
  )]
  dire_parTab$ItemID <- new_item_names

  testDat <- data.frame(
    test = "comp",
    subtest = "main",
    location = 0,
    scale = 1
  )

  mmlcomp_dire <- Dire_mml(
    formula = stats::reformulate(main_vars, response = "comp"),
    stuItems = stuItems,
    stuDat = stuDat,
    idVar = "subject",
    dichotParamTab = dire_parTab,
    testScale = testDat,
    X = X,
    main_vars = main_vars,
    pca_type = "cov",
    use_residual = TRUE,
    var_threshold = 0.80
  )

  set.seed(seed + 4L)
  PVs_dire <- Dire_drawPVs(
    x = mmlcomp_dire,
    npv = npv,
    pvVariableNameSuffix = "_dire"
  )

  # -------------------------------------------------------------------------
  # Put every PV data frame into the original sim_a1 subject order.
  # -------------------------------------------------------------------------
  order_by_subject <- function(dat, n) {
    dat <- as.data.frame(dat)

    if (!"id" %in% names(dat)) {
      stop("The plausible-value data must contain an id column.")
    }

    row_order <- match(
      as.character(seq_len(n)),
      as.character(dat$id)
    )

    if (anyNA(row_order)) {
      stop("Some subject IDs are missing from the plausible-value data.")
    }

    dat <- dat[row_order, , drop = FALSE]
    rownames(dat) <- NULL
    dat
  }

  debias_datPVs <- order_by_subject(PVs_debias$datPVs, n)
  emm_datPVs <- order_by_subject(PVs_emm$datPVs, n)
  dire_datPVs <- order_by_subject(PVs_dire$data, n)

  debias_PV1 <- as.numeric(debias_datPVs[["_farl1"]])
  emm_PV1 <- as.numeric(emm_datPVs[["_farl1"]])
  dire_PV1 <- as.numeric(dire_datPVs[["_dire1"]])

  stopifnot(
    length(debias_PV1) == n,
    length(emm_PV1) == n,
    length(dire_PV1) == n,
    all(is.finite(debias_PV1)),
    all(is.finite(emm_PV1)),
    all(is.finite(dire_PV1))
  )

  # Numerical comparison based on the first plausible value.
  comparison <- data.frame(
    method = c("FARLR-Debias", "FARLR-EMM", "DIRE"),
    correlation = c(
      stats::cor(theta, debias_PV1),
      stats::cor(theta, emm_PV1),
      stats::cor(theta, dire_PV1)
    ),
    bias = c(
      mean(debias_PV1 - theta),
      mean(emm_PV1 - theta),
      mean(dire_PV1 - theta)
    ),
    RMSE = c(
      sqrt(mean((debias_PV1 - theta)^2)),
      sqrt(mean((emm_PV1 - theta)^2)),
      sqrt(mean((dire_PV1 - theta)^2))
    )
  )

  # -------------------------------------------------------------------------
  # Plot all three methods in one figure.
  # -------------------------------------------------------------------------
  if (!is.null(plot_file)) {
    grDevices::png(
      filename = plot_file,
      width = 3000,
      height = 2400,
      res = 300
    )
  }

  old_par <- graphics::par(no.readonly = TRUE)
  graphics::par(
    mar = c(5.5, 5.5, 5, 2) + 0.1,
    las = 1,
    cex.axis = 1.15,
    cex.lab = 1.25
  )

  plot_range <- range(
    theta,
    debias_PV1,
    emm_PV1,
    dire_PV1,
    finite = TRUE
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    xlim = plot_range,
    ylim = plot_range,
    xlab = "True latent trait",
    ylab = "First plausible value",
    main = paste0(
      "Plausible Value Comparison\n",
      "DIRE, FARLR-EMM, and FARLR-Debias"
    ),
    cex.main = 1.4
  )

  graphics::grid(
    col = "gray90",
    lty = 1
  )

  graphics::points(
    theta,
    dire_PV1,
    pch = 16,
    cex = 0.5,
    col = grDevices::adjustcolor("steelblue", alpha.f = 0.22)
  )

  graphics::points(
    theta,
    emm_PV1,
    pch = 16,
    cex = 0.5,
    col = grDevices::adjustcolor("forestgreen", alpha.f = 0.22)
  )

  graphics::points(
    theta,
    debias_PV1,
    pch = 16,
    cex = 0.5,
    col = grDevices::adjustcolor("darkorange", alpha.f = 0.22)
  )

  graphics::abline(
    stats::lm(dire_PV1 ~ theta),
    col = "steelblue",
    lwd = 3
  )

  graphics::abline(
    stats::lm(emm_PV1 ~ theta),
    col = "forestgreen",
    lwd = 3
  )

  graphics::abline(
    stats::lm(debias_PV1 ~ theta),
    col = "darkorange",
    lwd = 3
  )

  graphics::abline(
    0,
    1,
    col = "gray35",
    lwd = 2,
    lty = 2
  )

  graphics::legend(
    "bottomright",
    legend = c(
      "DIRE regression",
      "FARLR-EMM regression",
      "FARLR-Debias regression",
      "Identity"
    ),
    col = c("steelblue", "forestgreen", "darkorange", "gray35"),
    lty = c(1, 1, 1, 2),
    lwd = c(3, 3, 3, 2),
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.88),
    box.col = "gray70",
    inset = 0.02,
    cex = 1.05
  )

  graphics::par(old_par)

  if (!is.null(plot_file)) {
    grDevices::dev.off()
  }

  list(
    mmlcomp = list(
      FARLR_Debias = mmlcomp_debias,
      FARLR_EMM = mmlcomp_emm,
      DIRE = mmlcomp_dire
    ),
    PVs = list(
      FARLR_Debias = PVs_debias,
      FARLR_EMM = PVs_emm,
      DIRE = PVs_dire
    ),
    ordered_datPVs = list(
      FARLR_Debias = debias_datPVs,
      FARLR_EMM = emm_datPVs,
      DIRE = dire_datPVs
    ),
    comparison = comparison,
    plot_file = plot_file
  )
}
#
#
# # Run the complete test and display the plot in the active graphics device.
# # This also works with the original simulation script, where theta exists as
# # a separate object but was not included in sim_a1.
# theta_test <- if (!is.null(sim_a1$theta)) sim_a1$theta else theta
#
# test_result <- mml_test_all(
#   sim_data = sim_a1,
#   theta = theta_test,
#   main = c(1, 2, 15, 29, 45),
#   npv = 10L,
#   seed = 2026L
# )
#
# test_result$comparison
#
#
# # # To save the same plot, run:
test_result <- mml_test_all(
  sim_data = sim_a1,
  theta = theta_test,
  main = c(1, 2, 15, 29, 45),
  npv = 10L,
  seed = 2026L,
  plot_file = "DIRE_FARLR_EMM_Debias_PV_comparison.png"
)
