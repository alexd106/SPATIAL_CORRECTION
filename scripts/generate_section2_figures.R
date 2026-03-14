# ============================================================================
# Generate figures and summary data for Report Section 2
# (Joint multi-bench models for unreplicated designs)
#
# Uses spatial type 6 (edge effects, centred) with 4 benches spanning
# very weak to very strong spatial intensity.
#
# Shared colour scales used across bench panels in spatial heatmaps.
#
# Produces:
#   docs/figures/s2_true_spatial.png
#   docs/figures/s2_mgcv_joint_spatial.png
#   docs/figures/s2_sommer_joint_spatial.png
#   docs/figures/s2_scatter_estimates_vs_true.png
#   docs/figures/s2_blup_shrinkage.png
#   output/section2/comparison_summary.csv
#   output/section2/residual_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("scripts/simulate_spatial_data.R")
source("scripts/fit_spatial_models.R")

fig_dir <- "docs/figures"
out_dir <- "output/section2"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ============================================================================
# 1. Simulation parameters — 4 benches, very low → very high intensity
# ============================================================================

# After centring, type-6 range ≈ ±0.38 × intensity:
#   intensity  3 → ±1.1  (very low, ~9% of genotype SD 12)
#   intensity 15 → ±5.7  (intermediate-low)
#   intensity 35 → ±13   (intermediate-high)
#   intensity 60 → ±23   (very high, ~2× genotype SD)
bench_intensities <- c(15, 25, 40, 60)
n_bench <- length(bench_intensities)

cat("=== Simulating", n_bench, "-bench data with edge-effect spatial pattern ===\n")

sim <- simulate_field_trial(
  n_bench           = n_bench,
  n_rows            = 10,
  n_cols            = 30,
  mean_pheno        = 55,
  sd_geno           = 12,
  spatial_type      = 6,
  spatial_intensity = 30,
  spatial_type_per_bench      = rep(6, n_bench),
  spatial_intensity_per_bench = bench_intensities,
  spatial_scale     = 1,
  sd_error          = 4,
  pheno_min         = 0,
  pheno_max         = 100,
  pheno_name        = "BNI",
  save_csv          = FALSE,
  seed              = 42
)

dfr          <- sim$data
true_effects <- sim$true_effects
pheno        <- sim$pheno_name
benches      <- sort(unique(dfr$Bench))

cat("Rows:", nrow(dfr), "  Genotypes:", length(unique(dfr$Genotype)),
    "  Benches:", paste(benches, collapse = ", "), "\n")

for (i in seq_along(benches)) {
  rng <- range(dfr$True_Spatial[dfr$Bench == benches[i]])
  cat(sprintf("  %s (intensity %d): True_Spatial [%.1f, %.1f]\n",
              benches[i], bench_intensities[i], rng[1], rng[2]))
}

# Shared spatial scale across all benches (for comparison figures)
shared_spatial_limits <- range(dfr$True_Spatial, na.rm = TRUE)
cat(sprintf("\nShared spatial scale: [%.1f, %.1f]\n",
            shared_spatial_limits[1], shared_spatial_limits[2]))

# ============================================================================
# 2. Figure: True spatial patterns (shared colour scale)
# ============================================================================

cat("\n=== Generating true spatial pattern figure ===\n")

true_plots <- lapply(seq_along(benches), function(i) {
  bn <- benches[i]
  bd <- dfr[dfr$Bench == bn, ]
  plot_heatmap(bd, "Row", "Col", "True_Spatial",
               title  = paste0(bn, "  (intensity = ", bench_intensities[i], ")"),
               limits = shared_spatial_limits)
})
p_true <- wrap_plots(true_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "True Spatial Patterns by Bench (Type 6: Edge Effects)",
    subtitle = paste0("Shared colour scale across benches. Intensity increases from ",
                      bench_intensities[1], " (very low) to ",
                      bench_intensities[n_bench], " (very high)."),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )
ggsave(file.path(fig_dir, "s2_true_spatial.png"), p_true,
       width = 14, height = 8, dpi = 150)

# ============================================================================
# 3. Approach 1: Uncorrected means
# ============================================================================

cat("\n=== Approach 1: Uncorrected means ===\n")

uncorrected <- aggregate(dfr[[pheno]], by = list(Genotype = dfr$Genotype), FUN = mean)
colnames(uncorrected) <- c("Genotype", "Mean")

merged_unc <- merge(uncorrected, true_effects, by = "Genotype")
r_unc    <- cor(merged_unc$Mean,  merged_unc$True_BNI)
rmse_unc <- sqrt(mean((merged_unc$Mean - merged_unc$True_BNI)^2))
cat(sprintf("  r = %.4f, RMSE = %.2f\n", r_unc, rmse_unc))

# ============================================================================
# 4. Approach 2: Per-bench SpATS BLUPs averaged
# ============================================================================

cat("\n=== Approach 2: Per-bench SpATS BLUPs averaged ===\n")

blup_list <- list()
for (bn in benches) {
  cat(sprintf("  Processing %s...\n", bn))
  bench_data <- dfr[dfr$Bench == bn, ]
  r <- tryCatch(
    run_SpATS_bench(bench_data, pheno = pheno, gt_col = "Genotype",
                    row_col = "Row", col_col = "Col", output_type = "BLUPs"),
    error = function(e) { message("  SpATS error: ", e$message); list(blups = NULL) }
  )
  if (!is.null(r$blups)) {
    r$blups$Bench <- bn
    blup_list[[bn]] <- r$blups
  }
}

if (length(blup_list) > 0) {
  all_blups    <- do.call(rbind, blup_list)
  avg_blups    <- aggregate(BLUP ~ Genotype, data = all_blups, FUN = mean)
  merged_blup  <- merge(avg_blups, true_effects, by = "Genotype")
  r_blup    <- cor(merged_blup$BLUP, merged_blup$True_BNI)
  rmse_blup <- sqrt(mean((merged_blup$BLUP - merged_blup$True_BNI)^2))
  cat(sprintf("  r = %.4f, RMSE = %.2f\n", r_blup, rmse_blup))
} else {
  r_blup <- NA;  rmse_blup <- NA
}

# ============================================================================
# 5. BLUP shrinkage figure
# ============================================================================

cat("\n=== Generating BLUP shrinkage figure ===\n")

if (length(blup_list) > 0) {
  all_blups_merged <- merge(all_blups, true_effects, by = "Genotype")
  all_blups_merged$bench_label <- factor(
    all_blups_merged$Bench,
    levels = benches,
    labels = paste0(benches, "\n(intensity = ", bench_intensities, ")")
  )
  true_sd <- round(sd(true_effects$True_BNI), 1)
  blup_sd <- round(sd(all_blups$BLUP), 1)

  p_shrink <- ggplot(all_blups_merged, aes(x = True_BNI, y = BLUP)) +
    geom_point(alpha = 0.35, size = 1, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    facet_wrap(~ bench_label, ncol = 2) +
    labs(
      title    = "Per-Bench SpATS BLUPs vs True Genotype Values",
      subtitle = paste0("BLUPs are compressed toward the mean regardless of spatial severity ",
                        "(true SD = ", true_sd, ", BLUP SD = ", blup_sd, ")"),
      x = "True Genotype Value", y = "SpATS BLUP"
    ) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 9))
  ggsave(file.path(fig_dir, "s2_blup_shrinkage.png"), p_shrink,
         width = 10, height = 8, dpi = 150)
}

# ============================================================================
# 6. Approach 3: Joint mgcv BLUEs
# ============================================================================

cat("\n=== Approach 3: Joint mgcv BLUEs ===\n")

mgcv_res <- tryCatch(
  run_mgcv_joint(dfr, pheno = pheno, gt_col = "Genotype",
                 row_col = "Row", col_col = "Col", bench_col = "Bench"),
  error = function(e) { message("  Joint mgcv error: ", e$message); list(blues = NULL) }
)

if (!is.null(mgcv_res$blues)) {
  merged_mgcv  <- merge(mgcv_res$blues, true_effects, by = "Genotype")
  r_mgcv    <- cor(merged_mgcv$BLUE, merged_mgcv$True_BNI)
  rmse_mgcv <- sqrt(mean((merged_mgcv$BLUE - merged_mgcv$True_BNI)^2))
  cat(sprintf("  r = %.4f, RMSE = %.2f\n", r_mgcv, rmse_mgcv))
} else { r_mgcv <- NA;  rmse_mgcv <- NA }

# ============================================================================
# 7. Approach 4: Joint sommer BLUEs
# ============================================================================

cat("\n=== Approach 4: Joint sommer BLUEs ===\n")

sommer_res <- tryCatch(
  run_sommer_joint(dfr, pheno = pheno, gt_col = "Genotype",
                   row_col = "Row", col_col = "Col", bench_col = "Bench"),
  error = function(e) { message("  Joint sommer error: ", e$message); list(blues = NULL) }
)

if (!is.null(sommer_res$blues)) {
  merged_som  <- merge(sommer_res$blues, true_effects, by = "Genotype")
  r_som    <- cor(merged_som$BLUE, merged_som$True_BNI)
  rmse_som <- sqrt(mean((merged_som$BLUE - merged_som$True_BNI)^2))
  cat(sprintf("  r = %.4f, RMSE = %.2f\n", r_som, rmse_som))
} else { r_som <- NA;  rmse_som <- NA }

# ============================================================================
# 8. Figure: mgcv joint spatial trends (shared scale)
# ============================================================================

cat("\n=== Generating joint spatial trend figures ===\n")

# Shared scale: combine estimated spatial from both models
all_spatial <- c(mgcv_res$spatial, sommer_res$spatial)
shared_est_limits <- range(all_spatial, na.rm = TRUE)

if (!is.null(mgcv_res$spatial)) {
  dfr$mgcv_spatial <- mgcv_res$spatial
  mgcv_plots <- lapply(seq_along(benches), function(i) {
    bn <- benches[i]
    bd <- dfr[dfr$Bench == bn, ]
    plot_heatmap(bd, "Row", "Col", "mgcv_spatial",
                 title  = paste0(bn, "  (intensity = ", bench_intensities[i], ")"),
                 limits = shared_est_limits)
  })
  p_mgcv <- wrap_plots(mgcv_plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "mgcv Joint Model: Estimated Spatial Trends",
      subtitle = "Shared colour scale across benches",
      theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                       plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )
  ggsave(file.path(fig_dir, "s2_mgcv_joint_spatial.png"), p_mgcv,
         width = 14, height = 8, dpi = 150)
}

# ============================================================================
# 9. Figure: sommer joint spatial trends (shared scale)
# ============================================================================

if (!is.null(sommer_res$spatial)) {
  dfr$sommer_spatial <- sommer_res$spatial
  sommer_plots <- lapply(seq_along(benches), function(i) {
    bn <- benches[i]
    bd <- dfr[dfr$Bench == bn, ]
    plot_heatmap(bd, "Row", "Col", "sommer_spatial",
                 title  = paste0(bn, "  (intensity = ", bench_intensities[i], ")"),
                 limits = shared_est_limits)
  })
  p_som <- wrap_plots(sommer_plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "sommer Joint Model: Estimated Spatial Trends",
      subtitle = "Shared colour scale across benches",
      theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                       plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )
  ggsave(file.path(fig_dir, "s2_sommer_joint_spatial.png"), p_som,
         width = 14, height = 8, dpi = 150)
}

# ============================================================================
# 10. Figure: Scatter plots — estimates vs true values
# ============================================================================

cat("\n=== Generating scatter comparison figure ===\n")

scatter_list <- list()
merged_unc$Method <- "Uncorrected Mean"
scatter_list[["unc"]] <- data.frame(
  Genotype = merged_unc$Genotype, Estimate = merged_unc$Mean,
  True = merged_unc$True_BNI, Method = "Uncorrected Mean", stringsAsFactors = FALSE)

if (!is.na(r_blup))
  scatter_list[["spats"]] <- data.frame(
    Genotype = merged_blup$Genotype, Estimate = merged_blup$BLUP,
    True = merged_blup$True_BNI, Method = "Per-Bench SpATS BLUPs (avg)", stringsAsFactors = FALSE)

if (!is.na(r_mgcv))
  scatter_list[["mgcv"]] <- data.frame(
    Genotype = merged_mgcv$Genotype, Estimate = merged_mgcv$BLUE,
    True = merged_mgcv$True_BNI, Method = "Joint mgcv BLUEs", stringsAsFactors = FALSE)

if (!is.na(r_som))
  scatter_list[["sommer"]] <- data.frame(
    Genotype = merged_som$Genotype, Estimate = merged_som$BLUE,
    True = merged_som$True_BNI, Method = "Joint sommer BLUEs", stringsAsFactors = FALSE)

all_scatter <- do.call(rbind, scatter_list)
rownames(all_scatter) <- NULL
all_scatter$Method <- factor(all_scatter$Method,
  levels = c("Uncorrected Mean", "Per-Bench SpATS BLUPs (avg)",
             "Joint mgcv BLUEs", "Joint sommer BLUEs"))

scatter_p <- ggplot(all_scatter, aes(x = True, y = Estimate)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  facet_wrap(~ Method, scales = "free_y") +
  labs(title = "Estimated vs True Genotype Values",
       x = "True Genotype Value", y = "Estimate") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "s2_scatter_estimates_vs_true.png"), scatter_p,
       width = 10, height = 8, dpi = 150)

# ============================================================================
# 11. Per-bench residual SD and summary
# ============================================================================

cat("\n=== Per-bench residual SD ===\n")

resid_summary <- data.frame(Bench = character(), Intensity = integer(),
                             Method = character(), Residual_SD = numeric(),
                             stringsAsFactors = FALSE)

for (i in seq_along(benches)) {
  bn  <- benches[i]
  idx <- dfr$Bench == bn
  if (!is.null(mgcv_res$residuals))
    resid_summary <- rbind(resid_summary, data.frame(
      Bench = bn, Intensity = bench_intensities[i],
      Method = "mgcv_joint", Residual_SD = sd(mgcv_res$residuals[idx], na.rm = TRUE)))
  if (!is.null(sommer_res$residuals))
    resid_summary <- rbind(resid_summary, data.frame(
      Bench = bn, Intensity = bench_intensities[i],
      Method = "sommer_joint", Residual_SD = sd(sommer_res$residuals[idx], na.rm = TRUE)))
}

if (nrow(resid_summary) > 0) {
  print(resid_summary, row.names = FALSE)
  write.csv(resid_summary, file.path(out_dir, "residual_summary.csv"), row.names = FALSE)
}

cat("\n=== Summary ===\n")
summary_df <- data.frame(
  Approach = c("Uncorrected mean", "Per-bench SpATS BLUPs (avg)",
               "Joint mgcv BLUEs", "Joint sommer BLUEs"),
  r_true = c(r_unc, r_blup, r_mgcv, r_som),
  RMSE   = c(rmse_unc, rmse_blup, rmse_mgcv, rmse_som),
  stringsAsFactors = FALSE
)
print(summary_df, row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "comparison_summary.csv"), row.names = FALSE)

cat("\n=== Section 2 figures and data generated ===\n")
cat("Figures in:", fig_dir, "\n")
cat("Data in:",    out_dir, "\n")
