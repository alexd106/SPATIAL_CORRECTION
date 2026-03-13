# ============================================================================
# Generate figures and summary data for Report Section 1
# (Per-bench spatial correction with replicated data — wheatdata)
#
# Produces:
#   docs/figures/s1_wheatdata_spats_diag.png
#   docs/figures/s1_wheatdata_method_comparison.png
#   docs/figures/s1_mgcv_psre_diag.png
#   output/section1/BLUEs.csv
#   output/section1/mgcv_variants_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("scripts/SpatialCorrectionFlexible.R")

fig_dir <- "docs/figures"
out_dir <- "output/section1"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ============================================================================
# 1. Load wheatdata
# ============================================================================

dfr <- read_input("data/wheatdata.rda", "wheatdata")
pheno   <- "yield"
gt_col  <- "geno"
row_col <- "row"
col_col <- "col"

cat("wheatdata loaded:", nrow(dfr), "rows,", length(unique(dfr[[gt_col]])), "genotypes\n")

# Prepare factor columns
dfr$row_f <- as.factor(dfr[[row_col]])
dfr$col_f <- as.factor(dfr[[col_col]])

# ============================================================================
# 2. Run all three base methods (SpATS, mgcv, sommer)
# ============================================================================

cat("\n=== Running SpATS ===\n")
res_spats <- run_SpATS_bench(dfr, pheno = pheno, gt_col = gt_col,
                              row_col = row_col, col_col = col_col,
                              output_type = "BLUEs")

cat("\n=== Running mgcv (default) ===\n")
res_mgcv <- run_mgcv_bench(dfr, pheno = pheno, gt_col = gt_col,
                            row_col = row_col, col_col = col_col,
                            output_type = "BLUEs", smoother_type = "te_default")

cat("\n=== Running sommer ===\n")
res_sommer <- run_sommer_bench(dfr, pheno = pheno, gt_col = gt_col,
                                row_col = row_col, col_col = col_col,
                                output_type = "BLUEs")

# ============================================================================
# 3. SpATS diagnostic figure
# ============================================================================

cat("\n=== Generating SpATS diagnostic ===\n")
p_spats_diag <- plot_diagnostics(dfr, res_spats, pheno, "SpATS",
                                  row_col, col_col)
ggsave(file.path(fig_dir, "s1_wheatdata_spats_diag.png"),
       p_spats_diag, width = 14, height = 9, dpi = 150)

# ============================================================================
# 4. Three-method comparison figure
# ============================================================================

cat("\n=== Generating method comparison ===\n")
results_list <- list(SpATS = res_spats, mgcv = res_mgcv, sommer = res_sommer)
p_compare <- plot_comparison(dfr, results_list, pheno, row_col, col_col)
ggsave(file.path(fig_dir, "s1_wheatdata_method_comparison.png"),
       p_compare, width = 14, height = 10, dpi = 150)

# ============================================================================
# 5. Run all mgcv variants + SpATS + sommer (7 methods)
# ============================================================================

cat("\n=== Running all mgcv variants ===\n")
mgcv_variants <- c("te_default", "te_ps", "te_ps_re", "te_highk", "te_ad")
variant_names <- c("mgcv", "mgcv_ps", "mgcv_ps_re", "mgcv_highk", "mgcv_ad")

all_results <- list(SpATS = res_spats)

for (i in seq_along(mgcv_variants)) {
  cat(sprintf("  Running %s...\n", variant_names[i]))
  all_results[[variant_names[i]]] <- tryCatch(
    run_mgcv_bench(dfr, pheno = pheno, gt_col = gt_col,
                   row_col = row_col, col_col = col_col,
                   output_type = "BLUEs", smoother_type = mgcv_variants[i]),
    error = function(e) { message("  FAILED: ", e$message); NULL }
  )
}

all_results[["sommer"]] <- res_sommer

# ============================================================================
# 6. mgcv_ps_re diagnostic figure
# ============================================================================

cat("\n=== Generating mgcv_ps_re diagnostic ===\n")
p_psre_diag <- plot_diagnostics(dfr, all_results[["mgcv_ps_re"]], pheno,
                                 "mgcv_ps_re", row_col, col_col)
ggsave(file.path(fig_dir, "s1_mgcv_psre_diag.png"),
       p_psre_diag, width = 14, height = 9, dpi = 150)

# ============================================================================
# 7. Summary metrics table for all 7 variants
# ============================================================================

cat("\n=== Computing summary metrics ===\n")
spats_blues <- res_spats$blues
spats_vals  <- spats_blues$BLUE
names(spats_vals) <- spats_blues$Genotype

summary_rows <- list()
for (mn in names(all_results)) {
  r <- all_results[[mn]]
  if (is.null(r)) next

  # Residual SD
  rsd <- sd(r$residuals, na.rm = TRUE)

  # Correlation with SpATS BLUEs
  if (mn == "SpATS") {
    r_spats <- NA
  } else {
    blues_mn <- r$blues
    merged <- merge(spats_blues[, c("Genotype", "BLUE")],
                    blues_mn[, c("Genotype", "BLUE")],
                    by = "Genotype", suffixes = c("_SpATS", paste0("_", mn)))
    r_spats <- cor(merged$BLUE_SpATS, merged[[ncol(merged)]], use = "complete.obs")
  }

  summary_rows <- c(summary_rows, list(data.frame(
    Method = mn, Residual_SD = round(rsd, 2),
    r_with_SpATS = round(r_spats, 4),
    stringsAsFactors = FALSE
  )))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "mgcv_variants_summary.csv"),
          row.names = FALSE)

cat("\nSummary:\n")
print(summary_df, row.names = FALSE)

# ============================================================================
# 8. Save BLUEs from all 3 base methods
# ============================================================================

blues_out <- res_spats$blues[, c("Genotype", "BLUE")]
names(blues_out)[2] <- "BLUE_SpATS"
blues_out$BLUE_mgcv   <- res_mgcv$blues$BLUE[match(blues_out$Genotype, res_mgcv$blues$Genotype)]
blues_out$BLUE_sommer <- res_sommer$blues$BLUE[match(blues_out$Genotype, res_sommer$blues$Genotype)]
write.csv(blues_out, file.path(out_dir, "BLUEs.csv"), row.names = FALSE)

cat("\n=== Section 1 figures and data generated ===\n")
cat("Figures in:", fig_dir, "\n")
cat("Data in:", out_dir, "\n")
