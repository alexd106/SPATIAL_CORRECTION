# ============================================================================
# Generate figures and summary data for Report Appendix
# Tests joint mgcv model under two simple gradient patterns:
#   Type 2 — row gradient (half-sine, low at row 1, high at last row)
#   Type 3 — column gradient (half-sine, low at col 1, high at last col)
#
# Each pattern uses 4 benches with intensities 15, 25, 40, 60.
# Compared against: uncorrected means and per-bench SpATS BLUPs.
#
# Produces:
#   docs/figures/sa_true_spatial.png
#   docs/figures/sa_mgcv_spatial.png
#   docs/figures/sa_scatter.png
#   output/appendix/summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("scripts/simulate_spatial_data.R")
source("scripts/fit_spatial_models.R")

fig_dir <- "docs/figures"
out_dir <- "output/appendix"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

bench_intensities <- c(15, 25, 40, 60)
n_bench           <- length(bench_intensities)
patterns          <- list(
  list(type = 2, label = "Row gradient"),
  list(type = 3, label = "Column gradient"),
  list(type = 7, label = "Localised hotspot")
)

# ============================================================================
# Helper: run one pattern, return metrics + spatial estimates
# ============================================================================

run_pattern <- function(spatial_type, label) {
  cat(sprintf("\n=== Pattern: %s (type %d) ===\n", label, spatial_type))

  sim <- simulate_field_trial(
    n_bench                     = n_bench,
    n_rows                      = 10,
    n_cols                      = 30,
    mean_pheno                  = 55,
    sd_geno                     = 12,
    spatial_type                = spatial_type,
    spatial_intensity           = 30,
    spatial_type_per_bench      = rep(spatial_type, n_bench),
    spatial_intensity_per_bench = bench_intensities,
    spatial_scale               = 1,
    sd_error                    = 4,
    pheno_min                   = 0,
    pheno_max                   = 100,
    pheno_name                  = "BNI",
    save_csv                    = FALSE,
    seed                        = 42
  )

  dfr          <- sim$data
  true_effects <- sim$true_effects
  pheno        <- sim$pheno_name
  benches      <- sort(unique(dfr$Bench))
  true_col     <- paste0("True_", pheno)

  # --- Approach 1: Uncorrected means ---
  unc       <- aggregate(dfr[[pheno]], list(Genotype = dfr$Genotype), mean)
  colnames(unc)[2] <- "Est"
  unc_m     <- merge(unc, true_effects, by = "Genotype")
  r_unc     <- cor(unc_m$Est, unc_m[[true_col]])
  rmse_unc  <- sqrt(mean((unc_m$Est - unc_m[[true_col]])^2))
  cat(sprintf("  Uncorrected:        r = %.4f, RMSE = %.2f\n", r_unc, rmse_unc))

  # --- Approach 2: Per-bench SpATS BLUPs averaged ---
  cat("  Running per-bench SpATS BLUPs...\n")
  blup_list <- list()
  for (bn in benches) {
    r <- tryCatch(
      run_SpATS_bench(dfr[dfr$Bench == bn, ], pheno = pheno,
                      gt_col = "Genotype", row_col = "Row",
                      col_col = "Col", output_type = "BLUPs"),
      error = function(e) list(blups = NULL)
    )
    if (!is.null(r$blups)) blup_list[[bn]] <- r$blups
  }
  if (length(blup_list) > 0) {
    avg_b   <- aggregate(BLUP ~ Genotype, do.call(rbind, blup_list), mean)
    blup_m  <- merge(avg_b, true_effects, by = "Genotype")
    r_blup  <- cor(blup_m$BLUP, blup_m[[true_col]])
    rmse_b  <- sqrt(mean((blup_m$BLUP - blup_m[[true_col]])^2))
    cat(sprintf("  SpATS BLUPs (avg):  r = %.4f, RMSE = %.2f\n", r_blup, rmse_b))
  } else {
    r_blup <- NA; rmse_b <- NA
    blup_m <- NULL
  }

  # --- Approach 3: Joint mgcv BLUEs ---
  cat("  Running joint mgcv...\n")
  mgcv_res <- tryCatch(
    run_mgcv_joint(dfr, pheno = pheno, gt_col = "Genotype",
                   row_col = "Row", col_col = "Col", bench_col = "Bench"),
    error = function(e) { message("  mgcv error: ", e$message); list(blues = NULL) }
  )
  if (!is.null(mgcv_res$blues)) {
    mgcv_m  <- merge(mgcv_res$blues, true_effects, by = "Genotype")
    r_mgcv  <- cor(mgcv_m$BLUE, mgcv_m[[true_col]])
    rmse_mg <- sqrt(mean((mgcv_m$BLUE - mgcv_m[[true_col]])^2))
    cat(sprintf("  Joint mgcv BLUEs:   r = %.4f, RMSE = %.2f\n", r_mgcv, rmse_mg))
    dfr$mgcv_spatial <- mgcv_res$spatial
  } else {
    r_mgcv <- NA; rmse_mg <- NA; mgcv_m <- NULL
  }

  list(
    dfr          = dfr,
    benches      = benches,
    true_effects = true_effects,
    true_col     = true_col,
    mgcv_m       = mgcv_m,
    blup_m       = blup_m,
    unc_m        = unc_m,
    metrics      = data.frame(
      Pattern  = label,
      Approach = c("Uncorrected", "SpATS BLUPs (avg)", "Joint mgcv BLUEs"),
      r        = round(c(r_unc, r_blup, r_mgcv), 4),
      RMSE     = round(c(rmse_unc, rmse_b, rmse_mg), 2)
    )
  )
}

# ============================================================================
# Run both patterns
# ============================================================================

results <- lapply(patterns, function(p) run_pattern(p$type, p$label))

# ============================================================================
# Figure A1: True spatial patterns (both pattern types, shared scale per type)
# ============================================================================

cat("\n=== Generating Figure A1: true spatial patterns ===\n")

true_panels <- lapply(seq_along(patterns), function(pi) {
  res     <- results[[pi]]
  dfr     <- res$dfr
  benches <- res$benches
  lim     <- range(dfr$True_Spatial)

  lapply(seq_along(benches), function(i) {
    plot_heatmap(dfr[dfr$Bench == benches[i], ], "Row", "Col", "True_Spatial",
                 title  = paste0(benches[i], " (intensity = ", bench_intensities[i], ")"),
                 limits = lim)
  })
})

p_true <- (
  wrap_plots(true_panels[[1]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Row gradient",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) / (
  wrap_plots(true_panels[[2]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Column gradient",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) / (
  wrap_plots(true_panels[[3]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Localised hotspot",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) +
  plot_annotation(
    title    = "True Spatial Patterns",
    subtitle = "Intensities 15, 25, 40, 60 (left to right). Shared colour scale within each pattern type.",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

ggsave(file.path(fig_dir, "sa_true_spatial.png"), p_true,
       width = 16, height = 11, dpi = 150)

# ============================================================================
# Figure A2: mgcv estimated spatial trends (both pattern types)
# ============================================================================

cat("=== Generating Figure A2: mgcv spatial trends ===\n")

mgcv_panels <- lapply(seq_along(patterns), function(pi) {
  res     <- results[[pi]]
  dfr     <- res$dfr
  benches <- res$benches
  lim     <- range(dfr$mgcv_spatial, na.rm = TRUE)

  lapply(seq_along(benches), function(i) {
    plot_heatmap(dfr[dfr$Bench == benches[i], ], "Row", "Col", "mgcv_spatial",
                 title  = paste0(benches[i], " (intensity = ", bench_intensities[i], ")"),
                 limits = lim)
  })
})

p_mgcv <- (
  wrap_plots(mgcv_panels[[1]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Row gradient",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) / (
  wrap_plots(mgcv_panels[[2]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Column gradient",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) / (
  wrap_plots(mgcv_panels[[3]], ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Localised hotspot",
                    theme = theme(plot.title = element_text(size = 11, face = "bold")))
) +
  plot_annotation(
    title    = "Joint mgcv: Estimated Spatial Trends",
    subtitle = "Shared colour scale within each pattern type.",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

ggsave(file.path(fig_dir, "sa_mgcv_spatial.png"), p_mgcv,
       width = 16, height = 11, dpi = 150)

# ============================================================================
# Figure A3: Scatter plots — all three approaches for both pattern types
# ============================================================================

cat("=== Generating Figure A3: scatter plots ===\n")

make_scatter <- function(x, y, xlabel, ylabel, title) {
  r    <- round(cor(x, y, use = "complete.obs"), 3)
  rmse <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
  lims <- range(c(x, y), na.rm = TRUE)
  df   <- data.frame(x = x, y = y)
  ggplot(df, aes(x, y)) +
    geom_point(alpha = 0.3, size = 0.8, colour = "steelblue") +
    geom_abline(colour = "red", linetype = "dashed") +
    coord_equal(xlim = lims, ylim = lims) +
    annotate("text", x = lims[1], y = lims[2],
             label = paste0("r = ", r, "\nRMSE = ", rmse),
             hjust = 0, vjust = 1, size = 3, colour = "grey30") +
    labs(title = title, x = xlabel, y = ylabel) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"))
}

scatter_panels <- lapply(seq_along(patterns), function(pi) {
  res  <- results[[pi]]
  lab  <- patterns[[pi]]$label
  tc   <- res$true_col

  p1 <- make_scatter(res$unc_m[[tc]], res$unc_m$Est,
                     "True", "Uncorrected mean", paste0(lab, "\nUncorrected"))
  p2 <- if (!is.null(res$blup_m))
    make_scatter(res$blup_m[[tc]], res$blup_m$BLUP,
                 "True", "SpATS BLUP", paste0(lab, "\nSpATS BLUPs (avg)"))
  else ggplot() + theme_void() + labs(title = "SpATS BLUPs: unavailable")
  p3 <- if (!is.null(res$mgcv_m))
    make_scatter(res$mgcv_m[[tc]], res$mgcv_m$BLUE,
                 "True", "mgcv BLUE", paste0(lab, "\nJoint mgcv BLUEs"))
  else ggplot() + theme_void() + labs(title = "mgcv BLUEs: unavailable")

  list(p1, p2, p3)
})

p_scatter <- wrap_plots(
  c(scatter_panels[[1]], scatter_panels[[2]], scatter_panels[[3]]),
  ncol = 3
) +
  plot_annotation(
    title    = "Estimated vs True Genotype Values",
    subtitle = "Red dashed line = 1:1. Rows: row gradient, column gradient, localised hotspot.",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 9, colour = "grey40"))
  )

ggsave(file.path(fig_dir, "sa_scatter.png"), p_scatter,
       width = 12, height = 11, dpi = 150)

# ============================================================================
# Summary table
# ============================================================================

cat("=== Writing summary CSV ===\n")

summary_df <- do.call(rbind, lapply(results, function(r) r$metrics))
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(out_dir, "appendix_summary.csv"), row.names = FALSE)
cat("\n=== Summary ===\n")
print(summary_df)
cat("\n=== Appendix figures and data generated ===\n")
cat("Figures in:", fig_dir, "\nData in:", out_dir, "\n")
