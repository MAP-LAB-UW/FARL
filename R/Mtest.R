# # ========================================================================
# # Compare MFARLR and MDIRE using sim_m1
# # ========================================================================
#
# # ------------------------------------------------------------------------
# # 1. Fit MFARLR-GVEM
# # ------------------------------------------------------------------------
#
# mfarlr_fit <- MFARLR_GVEM(
#   X = sim_m1$X,
#   Y = sim_m1$Y,
#   parTab = sim_m1$parTab,
#   farlr_args = list(
#     n_sam = 5L,
#     lambda = c(0.05, 0.10),
#     main = sim_m1$main,
#     delta.criteria = 1e-2,
#     iter.max = 20L,
#     window.size = 5L,
#     progress = FALSE,
#     verbose = FALSE,
#     seed = 123,
#     proposal_inflation = 0.2,
#     n_sam_final = 10L
#   ),
#   selection_tol = 1e-8,
#   max_iter = 500L,
#   threshold = 1e-3,
#   ridge = 1e-6,
#   verbose = TRUE
# )
#
# # ------------------------------------------------------------------------
# # 2. Draw MFARLR plausible values
# # ------------------------------------------------------------------------
#
# mfarlr_PVs <- MFARLR_drawPVs(
#   mmlcomp = mfarlr_fit,
#   npv = 10L,
#   n_mc = 20000L,
#   theta = sim_m1$theta,
#   seed = 123,
#   progress = TRUE
# )
#
#
# # ------------------------------------------------------------------------
# # 3. Fit multidimensional DIRE
# # ------------------------------------------------------------------------
#
# dire_fit <- MDIRE_mml(
#   X = sim_m1$X,
#   Y = sim_m1$Y,
#   parTab = sim_m1$parTab,
#   main = sim_m1$main,
#   pca_type = "cov",
#   use_residual = TRUE,
#   var_threshold = 0.90,
#   calcCor = TRUE
# )
#
#
# # ------------------------------------------------------------------------
# # 4. Draw DIRE plausible values
# # ------------------------------------------------------------------------
#
# dire_PVs <- MDIRE_drawPVs(
#   object = dire_fit,
#   npv = 10L,
#   group_vars = sim_m1$main
# )
#
#
# # ------------------------------------------------------------------------
# # 5. Sort both PV datasets into the original sim_m1 subject order
# # ------------------------------------------------------------------------
#
# N <- nrow(sim_m1$X)
# D <- ncol(sim_m1$theta)
# subject_order <- as.character(seq_len(N))
#
# mfarlr_dat <- mfarlr_PVs$datPVs
# mfarlr_order <- match(
#   subject_order,
#   as.character(mfarlr_dat$id)
# )
#
# if (anyNA(mfarlr_order)) {
#   stop("Some subjects are missing from the MFARLR PV results.")
# }
#
# mfarlr_dat <- mfarlr_dat[mfarlr_order, , drop = FALSE]
# rownames(mfarlr_dat) <- NULL
#
# dire_dat <- dire_PVs$data
# dire_order <- match(
#   subject_order,
#   as.character(dire_dat$id)
# )
#
# if (anyNA(dire_order)) {
#   stop("Some subjects are missing from the DIRE PV results.")
# }
#
# dire_dat <- dire_dat[dire_order, , drop = FALSE]
# rownames(dire_dat) <- NULL
#
# stopifnot(
#   identical(as.character(mfarlr_dat$id), subject_order),
#   identical(as.character(dire_dat$id), subject_order)
# )
#
#
# # ------------------------------------------------------------------------
# # 6. Extract PV1 from every dimension
# # ------------------------------------------------------------------------
#
# mfarlr_first_PV <- sapply(
#   seq_len(D),
#   function(dimension) {
#     mfarlr_dat[[paste0("_mfarlr_D", dimension, "_PV1")]]
#   }
# )
#
# dire_first_PV <- sapply(
#   seq_len(D),
#   function(dimension) {
#     subtest <- dire_fit$subtest_names[dimension]
#     dire_dat[[paste0(subtest, "_dire1")]]
#   }
# )
#
# dimension_names <- paste0("Dimension", seq_len(D))
# colnames(mfarlr_first_PV) <- dimension_names
# colnames(dire_first_PV) <- dimension_names
#
#
# # ------------------------------------------------------------------------
# # 7. Plot MFARLR and DIRE together
# # ------------------------------------------------------------------------
#
# plot_MFARLR_MDIRE <- function(
#     theta,
#     mfarlr_first_PV,
#     dire_first_PV,
#     file = NULL
# ) {
#   theta <- as.matrix(theta)
#   mfarlr_first_PV <- as.matrix(mfarlr_first_PV)
#   dire_first_PV <- as.matrix(dire_first_PV)
#
#   if (!identical(dim(theta), dim(mfarlr_first_PV)) ||
#       !identical(dim(theta), dim(dire_first_PV))) {
#     stop("theta and both PV matrices must have identical dimensions.")
#   }
#
#   D <- ncol(theta)
#   total_panels <- D + 1L
#   panel_columns <- ceiling(sqrt(total_panels))
#   panel_rows <- ceiling(total_panels / panel_columns)
#
#   save_plot <- !is.null(file)
#   if (save_plot) {
#     grDevices::png(
#       filename = file,
#       width = 3000,
#       height = 2100,
#       res = 300
#     )
#   }
#
#   old_par <- graphics::par(no.readonly = TRUE)
#   on.exit({
#     graphics::par(old_par)
#     if (save_plot && grDevices::dev.cur() > 1L) {
#       grDevices::dev.off()
#     }
#   }, add = TRUE)
#
#   graphics::par(
#     mfrow = c(panel_rows, panel_columns),
#     mar = c(4.5, 4.5, 3, 1),
#     oma = c(0, 0, 5, 0),
#     cex.axis = 1.05,
#     cex.lab = 1.1
#   )
#
#   for (dimension in seq_len(D)) {
#     true_theta <- theta[, dimension]
#     mfarlr_pv <- mfarlr_first_PV[, dimension]
#     dire_pv <- dire_first_PV[, dimension]
#
#     plot_range <- range(
#       true_theta,
#       mfarlr_pv,
#       dire_pv,
#       finite = TRUE
#     )
#
#     graphics::plot(
#       NA_real_,
#       NA_real_,
#       xlim = plot_range,
#       ylim = plot_range,
#       xlab = paste0("True theta, dimension ", dimension),
#       ylab = "First plausible value",
#       main = paste0("Dimension ", dimension)
#     )
#
#     graphics::grid(
#       col = "gray90",
#       lty = 1
#     )
#
#     graphics::points(
#       true_theta,
#       dire_pv,
#       pch = 16,
#       cex = 0.45,
#       col = grDevices::adjustcolor("steelblue", alpha.f = 0.30)
#     )
#
#     graphics::points(
#       true_theta,
#       mfarlr_pv,
#       pch = 16,
#       cex = 0.45,
#       col = grDevices::adjustcolor("darkorange", alpha.f = 0.30)
#     )
#
#     graphics::abline(
#       stats::lm(dire_pv ~ true_theta),
#       col = "steelblue",
#       lwd = 3
#     )
#
#     graphics::abline(
#       stats::lm(mfarlr_pv ~ true_theta),
#       col = "darkorange",
#       lwd = 3
#     )
#
#     graphics::abline(
#       0,
#       1,
#       col = "gray35",
#       lwd = 2,
#       lty = 2
#     )
#   }
#
#   # Use one panel for a shared legend.
#   graphics::par(mar = c(0, 0, 0, 0))
#   graphics::plot.new()
#   graphics::legend(
#     "center",
#     title = "Method",
#     legend = c(
#       "DIRE regression",
#       "MFARLR regression",
#       "Identity"
#     ),
#     col = c("steelblue", "darkorange", "gray35"),
#     lty = c(1, 1, 2),
#     lwd = c(3, 3, 2),
#     bty = "n",
#     cex = 1.05
#   )
#
#   # Fill any unused layout panels after the shared legend.
#   unused_panels <- panel_rows * panel_columns - total_panels
#   if (unused_panels > 0L) {
#     for (panel in seq_len(unused_panels)) {
#       graphics::plot.new()
#     }
#   }
#
#   graphics::mtext(
#     "Comparison of MFARLR and DIRE Plausible Values",
#     side = 3,
#     outer = TRUE,
#     line = 2,
#     cex = 1.5,
#     font = 2
#   )
#
#   invisible(file)
# }
#
#
# # Display the comparison plot.
# plot_MFARLR_MDIRE(
#   theta = sim_m1$theta,
#   mfarlr_first_PV = mfarlr_first_PV,
#   dire_first_PV = dire_first_PV
# )
#
#
# # To save a larger PNG, use:
# plot_file <- file.path(getwd(), "MFARLR_DIRE_PV_comparison.png")
# plot_MFARLR_MDIRE(
#   theta = sim_m1$theta,
#   mfarlr_first_PV = mfarlr_first_PV,
#   dire_first_PV = dire_first_PV,
#   file = plot_file
# )
