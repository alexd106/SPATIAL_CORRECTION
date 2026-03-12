### R script: Spatial correction using mgcv::gamm() with tensor product P-splines
#
# Literature basis:
#   Rodriguez-Alvarez et al. (2018, TAG) — SpATS / SAP model
#   Piepho et al. (2015) — P-spline spatial models for field trials
#
# Two modelling modes (set separate_smoothers in CONFIG):
#
#   separate_smoothers = TRUE  (default, recommended)
#     One model per bench x trait:
#       pheno ~ te(row, col, bs=c("ps","ps"), k=k_dims)   [fixed spatial smooth]
#       + random: Genotype, rowf, colf                     [via gamm / nlme]
#     Each bench gets its own spatial surface — handles heterogeneous patterns.
#
#   separate_smoothers = FALSE
#     One pooled model per trait across all benches:
#       pheno ~ Bench + te(row, col, bs=c("ps","ps"), k=k_dims)
#       + random: Genotype, bench_rowf, bench_colf
#     Assumes the same spatial gradient shape on every bench; bench differences
#     absorbed by a fixed effect. Produces one cross-bench BLUP directly.
#
# Key CONFIG options:
#   bench_col          — column with bench identifiers; NULL = single field
#   bench_filter       — character vector to process a bench subset; NULL = all
#   separate_smoothers — TRUE / FALSE (see above)
#
# Outputs:
#   GAMM_plots.pdf            — 6-panel diagnostic + smooth heatmap per trait x bench
#   GAMM_fitted_values.csv    — plot-level: observed/fitted/residual/corrected/spatial trend
#   GAMM_BLUPs.csv            — genotype BLUPs per bench + cross-bench mean (separate)
#                               or single pooled BLUPs (pooled)
#   Outliers_GAMM_report.txt  — outlier log

library(mgcv)
library(nlme)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)

# ==============================================================================
# Helper functions
# ==============================================================================

#' Auto-compute k_dims for te() based on bench grid size
#' @param mydfr  data frame for a single bench
#' @param rn     name of row column
#' @param cn     name of column column
#' @return integer vector length 2: c(k_row, k_col)
auto_k_dims <- function(mydfr, rn, cn) {
  n_r <- length(unique(mydfr[[rn]]))
  n_c <- length(unique(mydfr[[cn]]))
  c(max(3L, n_r %/% 2L), max(3L, n_c %/% 2L))
}


#' Fit GAMM spatial model for one bench x one trait
#' @param mydfr   data frame (single bench)
#' @param pheno   response column name
#' @param gt      genotype column name
#' @param rn      row column name
#' @param cn      column column name
#' @param k_dims  integer vector c(k_row, k_col) for te() knots
#' @return gamm() model object (list with $gam and $lme)
fit_gamm_spatial <- function(mydfr, pheno, gt, rn, cn, k_dims) {
  mydfr$rowf      <- as.factor(mydfr[[rn]])
  mydfr$colf      <- as.factor(mydfr[[cn]])
  mydfr[[gt]]     <- as.factor(mydfr[[gt]])

  formula_obj <- as.formula(
    paste0(pheno, " ~ te(", rn, ", ", cn, ", bs = c('ps','ps'), k = c(",
           k_dims[1], ",", k_dims[2], "))")
  )

  model <- gamm(
    formula_obj,
    random = setNames(list(~1, ~1, ~1), c(gt, "rowf", "colf")),
    data   = mydfr,
    method = "REML"
  )

  # Attach prepared data to model for downstream use
  attr(model, "model_data") <- mydfr
  model
}


#' Detect outliers from GAMM conditional residuals
#' @param model   gamm() object
#' @param mydfr   data frame used to fit model
#' @param rn      row column name
#' @param cn      column column name
#' @param method  "IQR" or "SD"
#' @param k       multiplier for threshold
#' @return logical vector (TRUE = outlier), same length as nrow(mydfr)
detect_outliers_residuals <- function(model, mydfr, rn, cn,
                                       method = "IQR", k = 3) {
  # Handle both gamm() (list with $lme) and gam() objects
  resids <- if (!is.null(model$lme)) residuals(model$lme, level = 1) else residuals(model)

  if (method == "IQR") {
    Q1  <- quantile(resids, 0.25, na.rm = TRUE)
    Q3  <- quantile(resids, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower   <- Q1 - k * IQR_val
    upper   <- Q3 + k * IQR_val
  } else if (method == "SD") {
    mn    <- mean(resids, na.rm = TRUE)
    sd_r  <- sd(resids, na.rm = TRUE)
    lower <- mn - k * sd_r
    upper <- mn + k * sd_r
  } else {
    stop("method must be 'IQR' or 'SD'")
  }

  is_outlier <- !is.na(resids) & (resids < lower | resids > upper)

  n_out <- sum(is_outlier)
  cat("  Outliers detected:", n_out, "\n")
  if (n_out > 0) {
    cat("  Outlier residual values:", round(resids[is_outlier], 3), "\n")
    cat("  Outlier row positions:", mydfr[[rn]][is_outlier], "\n")
    cat("  Outlier col positions:", mydfr[[cn]][is_outlier], "\n")
  } else {
    cat("  No outliers found.\n")
  }

  is_outlier
}


#' Extract genotype BLUPs from a fitted gamm() object
#' @param model  gamm() object
#' @param mydfr  data frame used to fit (needs gt column as factor)
#' @param gt     genotype column name
#' @param n_rep  number of replicates per genotype per bench (default 1)
#' @return data.frame with columns: Genotype, predicted, deBLUP, StdErr, reliability
extract_blups_gamm <- function(model, mydfr, gt, n_rep = 1) {
  gt_ranef   <- ranef(model$lme)[[gt]]
  # gamm() prefixes fixed-effect names with "X"; intercept is the first element
  fe         <- fixef(model$lme)
  grand_mean <- fe[grep("Intercept", names(fe))[1]]

  ranef_dev  <- gt_ranef[, "(Intercept)"]   # deviation from grand mean
  adj_means  <- grand_mean + ranef_dev

  vc         <- VarCorr(model$lme)
  gt_hdr     <- which(rownames(vc) == paste0(gt, " ="))
  sigma2_g   <- as.numeric(vc[gt_hdr + 1L, "Variance"])
  sigma2_e   <- model$lme$sigma^2

  reliability <- sigma2_g / (sigma2_g + sigma2_e / n_rep)
  blup_se     <- sqrt(sigma2_g * (1 - reliability))

  # De-regressed BLUP: removes shrinkage so the value is approximately unbiased.
  # BLUP deviation = ranef_dev ≈ reliability * true_effect
  # => deBLUP = grand_mean + ranef_dev / reliability
  # Equivalent to BLUE from a fixed-genotype model in the balanced case.
  deblup <- as.numeric(grand_mean + ranef_dev / reliability)

  # ranef rownames can be "outer/level" — strip the "outer/" prefix if present
  geno_names <- sub("^[^/]+/", "", rownames(gt_ranef))

  data.frame(
    Genotype    = geno_names,
    predicted   = as.numeric(adj_means),
    deBLUP      = deblup,
    StdErr      = blup_se,
    reliability = reliability,
    stringsAsFactors = FALSE
  )
}


#' Compute empirical variogram of residuals
#' @param resids    numeric vector of residuals
#' @param rows      numeric vector of row positions
#' @param cols      numeric vector of column positions
#' @param n_bins    number of distance bins
#' @param max_pairs maximum number of random pairs to evaluate
#' @return data.frame with columns: dist_bin_mid, semivariance, n_pairs
compute_variogram <- function(resids, rows, cols, n_bins = 15, max_pairs = 5000) {
  n   <- length(resids)
  idx <- which(!is.na(resids))
  n_v <- length(idx)

  # All pairs
  pairs <- combn(n_v, 2)  # 2 x n_pairs matrix

  if (ncol(pairs) > max_pairs) {
    sel   <- sample(ncol(pairs), max_pairs)
    pairs <- pairs[, sel]
  }

  i <- pairs[1, ]
  j <- pairs[2, ]

  dists <- sqrt((rows[idx[i]] - rows[idx[j]])^2 + (cols[idx[i]] - cols[idx[j]])^2)
  gamma <- 0.5 * (resids[idx[i]] - resids[idx[j]])^2

  bins     <- cut(dists, breaks = n_bins)
  bin_mid  <- tapply(dists, bins, mean)
  semi_var <- tapply(gamma, bins, mean)
  n_pairs  <- tapply(gamma, bins, length)

  keep <- !is.na(bin_mid) & !is.na(semi_var)
  data.frame(
    dist_bin_mid = as.numeric(bin_mid[keep]),
    semivariance = as.numeric(semi_var[keep]),
    n_pairs      = as.integer(n_pairs[keep])
  )
}


#' Build 6-panel diagnostic plot page for one trait x one bench
#' @param model        gamm() object
#' @param mydfr        data frame used for fitting
#' @param blups        data frame from extract_blups_gamm()
#' @param pheno        response column name
#' @param bench_label  label for plot titles
#' @param rn           row column name
#' @param cn           column column name
#' @param gt           genotype column name
#' @return patchwork ggplot object (6 panels)
plot_spatial_diagnostics <- function(model, mydfr, blups, pheno, bench_label,
                                     rn, cn, gt) {
  # ---- Residuals and fitted values ----
  cond_resids <- residuals(model$lme, level = 1)
  fitted_vals <- fitted(model$lme)
  spatial_raw <- predict(model$gam)

  mydfr$.resid   <- cond_resids
  mydfr$.fitted  <- fitted_vals
  mydfr$.spatial <- spatial_raw

  row_var <- mydfr[[rn]]
  col_var <- mydfr[[cn]]

  base_theme <- theme_bw() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7))

  # Panel 1: Observed heatmap
  p1 <- ggplot(mydfr, aes(x = .data[[cn]], y = .data[[rn]], fill = .data[[pheno]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(name = pheno, na.value = "grey70") +
    scale_y_reverse() +
    labs(title = paste0("Observed: ", pheno, " [", bench_label, "]"),
         x = cn, y = rn) +
    base_theme

  # Panel 2: Fitted heatmap
  p2 <- ggplot(mydfr, aes(x = .data[[cn]], y = .data[[rn]], fill = .fitted)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(name = "Fitted", option = "plasma") +
    scale_y_reverse() +
    labs(title = paste0("Fitted: ", pheno, " [", bench_label, "]"),
         x = cn, y = rn) +
    base_theme

  # Panel 3: Residuals heatmap
  resid_lim <- max(abs(mydfr$.resid), na.rm = TRUE)
  p3 <- ggplot(mydfr, aes(x = .data[[cn]], y = .data[[rn]], fill = .resid)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colours = c("#2166ac", "white", "#d6604d"),
      limits  = c(-resid_lim, resid_lim),
      name    = "Resid"
    ) +
    scale_y_reverse() +
    labs(title = paste0("Residuals: ", pheno, " [", bench_label, "]"),
         x = cn, y = rn) +
    base_theme

  # Panel 4: Spatial surface (interpolated grid via predict.gam on fine grid)
  rn_seq  <- seq(min(row_var), max(row_var), length.out = 50)
  cn_seq  <- seq(min(col_var), max(col_var), length.out = 50)
  pred_grid <- expand.grid(setNames(list(rn_seq, cn_seq), c(rn, cn)))
  pred_grid$.spatial_pred <- predict(model$gam, newdata = pred_grid)

  p4 <- ggplot(pred_grid, aes(x = .data[[cn]], y = .data[[rn]],
                               fill = .spatial_pred)) +
    geom_tile() +
    geom_contour(aes(z = .spatial_pred), colour = "white", alpha = 0.4,
                 linewidth = 0.3) +
    scale_fill_viridis_c(name = "Spatial\ntrend", option = "magma") +
    scale_y_reverse() +
    labs(title = paste0("Spatial surface [", bench_label, "]"),
         x = cn, y = rn) +
    base_theme

  # Panel 5: BLUP ranked dotplot
  blups_ord <- blups[order(blups$predicted), ]
  blups_ord$rank <- seq_len(nrow(blups_ord))
  p5 <- ggplot(blups_ord, aes(x = rank, y = predicted)) +
    geom_point(size = 1.2, colour = "#2c7bb6") +
    geom_errorbar(aes(ymin = predicted - StdErr, ymax = predicted + StdErr),
                  width = 0, colour = "#2c7bb6", alpha = 0.6) +
    labs(title = paste0("BLUPs ranked [", bench_label, "]"),
         x = "Genotype rank", y = paste0("Adjusted mean: ", pheno)) +
    base_theme

  # Panel 6: Empirical variogram
  vgm <- compute_variogram(cond_resids, row_var, col_var)
  p6 <- ggplot(vgm, aes(x = dist_bin_mid, y = semivariance)) +
    geom_point(aes(size = n_pairs), colour = "#d7301f") +
    geom_line(colour = "#d7301f", linewidth = 0.5) +
    scale_size_continuous(name = "n pairs", range = c(1, 4)) +
    labs(title = paste0("Variogram of residuals [", bench_label, "]"),
         x = "Distance", y = "Semivariance") +
    base_theme

  # Compose with patchwork
  (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_annotation(
      title   = paste0(pheno, " — ", bench_label),
      theme   = theme(plot.title = element_text(size = 11, face = "bold"))
    )
}


#' Plot the fitted spatial smooth as a discrete heatmap on the actual plot grid
#'
#' The surface is centred at zero so that positive values (red) indicate plots
#' elevated by the spatial environment and negative values (blue) indicate
#' plots suppressed. This is the correction applied: plots with red tiles had
#' their phenotype adjusted downward; blue tiles upward.
#'
#' @param model        gamm() object (per-bench or pooled)
#' @param mydfr        data frame used to fit the model
#' @param rn           row column name
#' @param cn           column column name
#' @param bench_label  label shown in the plot title
#' @param pheno_label  fill legend label (default "Spatial\ntrend")
#' @return ggplot object — diverging blue–white–red scale, centred at 0
plot_smooth_heatmap <- function(model, mydfr, rn, cn,
                                 bench_label  = "",
                                 pheno_label  = "Spatial\ntrend") {
  # Handle both gamm() (list with $gam) and plain gam() objects
  gam_obj         <- if (!is.null(model$gam)) model$gam else model
  raw_pred        <- predict(gam_obj)
  spatial_centred <- raw_pred - mean(raw_pred, na.rm = TRUE)
  mydfr$.smooth   <- as.numeric(spatial_centred)

  lim <- max(abs(mydfr$.smooth), na.rm = TRUE)
  if (lim == 0) lim <- 1

  title_str <- if (nchar(bench_label) > 0)
    paste0("Fitted spatial smooth [", bench_label, "]")
  else
    "Fitted spatial smooth"

  ggplot(mydfr, aes(x = .data[[cn]], y = .data[[rn]], fill = .smooth)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = c("#2166ac", "white", "#d6604d"),
      limits  = c(-lim, lim),
      name    = pheno_label
    ) +
    scale_y_reverse(breaks = sort(unique(mydfr[[rn]]))) +
    scale_x_continuous(breaks = sort(unique(mydfr[[cn]]))) +
    labs(
      title    = title_str,
      subtitle = "Blue = suppressed by environment  |  Red = elevated by environment",
      x = cn, y = rn
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      axis.title    = element_text(size = 8),
      axis.text     = element_text(size = 7),
      legend.title  = element_text(size = 8),
      legend.text   = element_text(size = 7),
      panel.grid    = element_blank()
    )
}


# ==============================================================================
# Main pipeline: one trait x one bench (separate smoother)
# ==============================================================================

#' Run full spatial GAMM pipeline for one trait x one bench
#' @param mydfr         data frame (single bench, all plots)
#' @param pheno         response column name
#' @param gt            genotype column name
#' @param rn            row column name
#' @param cn            column column name
#' @param bench_label   bench identifier string
#' @param pmain         plot title prefix
#' @param outlier_method "IQR" or "SD"
#' @param outlier_k     multiplier for outlier threshold
#' @return list(dfrF = plot-level data frame, dfrP = BLUP data frame)
run_spatial_gamm <- function(mydfr, pheno, gt, rn, cn, bench_label, pmain,
                              outlier_method = "IQR", outlier_k = 3) {

  cat("\n--- Bench:", bench_label, "| Trait:", pheno, "---\n")

  # Step 1: Auto k_dims
  k_dims <- auto_k_dims(mydfr, rn, cn)
  cat("  k_dims:", k_dims[1], "x", k_dims[2], "\n")

  # Step 2: Initial fit
  cat("  Fitting initial GAMM...\n")
  model_init <- tryCatch(
    fit_gamm_spatial(mydfr, pheno, gt, rn, cn, k_dims),
    error = function(e) {
      cat("  ERROR in initial GAMM fit:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(model_init)) return(NULL)

  # Step 3: Detect outliers on conditional residuals
  cat("  Detecting outliers (method =", outlier_method, ", k =", outlier_k, ")...\n")
  is_out <- detect_outliers_residuals(model_init, mydfr, rn, cn,
                                       method = outlier_method,
                                       k      = outlier_k)

  # Log to report (already captured by sink in main block)
  cat("[REPORT] Bench:", bench_label, "| Trait:", pheno, "\n")
  cat("[REPORT] Method:", outlier_method, "| k:", outlier_k,
      "| Outliers:", sum(is_out), "\n")
  if (any(is_out)) {
    cat("[REPORT] Outlier row positions  :", mydfr[[rn]][is_out], "\n")
    cat("[REPORT] Outlier col positions  :", mydfr[[cn]][is_out], "\n")
    cat("[REPORT] Outlier observed values:",
        round(mydfr[[pheno]][is_out], 3), "\n")
  }
  cat("[REPORT] ---\n")

  # Step 4: Replace outliers with NA and re-fit
  mydfr_clean <- mydfr
  mydfr_clean[[pheno]][is_out] <- NA

  cat("  Fitting final GAMM on cleaned data...\n")
  model_final <- tryCatch(
    fit_gamm_spatial(mydfr_clean, pheno, gt, rn, cn, k_dims),
    error = function(e) {
      cat("  ERROR in final GAMM fit:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(model_final)) return(NULL)

  # Step 5: Extract BLUPs
  n_rep  <- 1   # unreplicated design
  blup_df <- tryCatch(
    extract_blups_gamm(model_final, mydfr_clean, gt, n_rep),
    error = function(e) {
      cat("  ERROR extracting BLUPs:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(blup_df)) return(NULL)

  # Step 6: Build per-plot output frame
  # Spatial components
  spatial_raw   <- predict(model_final$gam)
  spatial_trend <- spatial_raw - mean(spatial_raw, na.rm = TRUE)  # centred

  row_re_tbl <- ranef(model_final$lme)$rowf
  col_re_tbl <- ranef(model_final$lme)$colf

  # Match row/col random effects back to each plot
  mydfr_clean$rowf <- as.factor(mydfr_clean[[rn]])
  mydfr_clean$colf <- as.factor(mydfr_clean[[cn]])

  row_re <- row_re_tbl[match(as.character(mydfr_clean$rowf),
                              rownames(row_re_tbl)), 1]
  col_re <- col_re_tbl[match(as.character(mydfr_clean$colf),
                              rownames(col_re_tbl)), 1]
  row_re[is.na(row_re)] <- 0
  col_re[is.na(col_re)] <- 0

  corrected <- mydfr_clean[[pheno]] - spatial_trend - row_re - col_re
  fitted_v  <- fitted(model_final$lme)

  dfrF <- data.frame(
    Genotype = mydfr_clean[[gt]],
    Bench    = bench_label,
    row      = mydfr_clean[[rn]],
    col      = mydfr_clean[[cn]],
    stringsAsFactors = FALSE
  )
  dfrF[[paste0(pheno, "_Observed")]]     <- mydfr_clean[[pheno]]
  dfrF[[paste0(pheno, "_Fitted")]]       <- fitted_v
  dfrF[[paste0(pheno, "_Residual")]]     <- residuals(model_final$lme, level = 1)
  dfrF[[paste0(pheno, "_Corrected")]]    <- corrected
  dfrF[[paste0(pheno, "_SpatialTrend")]] <- spatial_trend

  # Mark original outlier positions as NA in Observed column
  dfrF[[paste0(pheno, "_Observed")]][is_out] <- NA

  # Step 7: BLUP data frame with bench label
  # BLUE = spatially corrected observation (genotype treated as fixed, no shrinkage)
  # For the unreplicated per-bench design (1 rep), BLUE = corrected phenotype value.
  blue_by_geno <- tapply(corrected, as.character(mydfr_clean[[gt]]), mean, na.rm = TRUE)

  dfrP <- data.frame(
    Genotype = blup_df$Genotype,
    Bench    = bench_label,
    stringsAsFactors = FALSE
  )
  dfrP[[pheno]]                         <- blup_df$predicted
  dfrP[[paste0(pheno, "_BLUE")]]        <- as.numeric(blue_by_geno[blup_df$Genotype])
  dfrP[[paste0(pheno, "_deBLUP")]]      <- blup_df$deBLUP
  dfrP[[paste0(pheno, "_StdErr")]]      <- blup_df$StdErr
  dfrP[[paste0(pheno, "_Reliability")]] <- blup_df$reliability

  # Step 8: Diagnostic plots (drawn to current PDF device)
  cat("  Generating diagnostic plots...\n")
  pg <- tryCatch(
    plot_spatial_diagnostics(model_final, mydfr_clean, blup_df, pheno,
                              bench_label, rn, cn, gt),
    error = function(e) {
      cat("  WARNING: Plot generation failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(pg)) print(pg)

  # Step 9: Spatial smooth heatmap (discrete, on actual plot grid)
  sh <- tryCatch(
    plot_smooth_heatmap(model_final, mydfr_clean, rn, cn, bench_label),
    error = function(e) {
      cat("  WARNING: Smooth heatmap failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(sh)) print(sh)

  list(dfrF = dfrF, dfrP = dfrP)
}


# ==============================================================================
# Pooled pipeline: one trait, shared smooth across all benches
# ==============================================================================

#' Fit GAMM with a single shared spatial smooth across multiple benches
#'
#' Bench is included as a fixed effect to absorb between-bench mean differences.
#' Row and column random effects are bench-specific to avoid confounding plots
#' that share the same grid coordinates across benches.
#'
#' Use this when you expect the spatial gradient to have the same shape and
#' intensity on every bench (e.g. identical equipment, same orientation).
#' When bench-specific patterns are likely, use fit_gamm_spatial() per bench.
#'
#' @param dfr       full data frame (multiple benches)
#' @param pheno     response column name
#' @param gt        genotype column name
#' @param rn        row column name
#' @param cn        column column name
#' @param bench_col bench column name
#' @param k_dims    integer vector c(k_row, k_col)
#' @return gamm() object
fit_gamm_pooled <- function(dfr, pheno, gt, rn, cn, bench_col, k_dims) {
  # Use gam() with bs="re" smooths — more numerically stable than gamm() when
  # pooling multiple benches with identical row/col coordinate grids.
  dfr$rowf         <- as.factor(dfr[[rn]])
  dfr$colf         <- as.factor(dfr[[cn]])
  dfr[[gt]]        <- as.factor(dfr[[gt]])
  dfr[[bench_col]] <- as.factor(dfr[[bench_col]])

  formula_obj <- as.formula(
    paste0(pheno, " ~ ", bench_col,
           " + te(", rn, ", ", cn, ", bs = c('ps','ps'), k = c(",
           k_dims[1], ",", k_dims[2], "))",
           " + s(", gt, ", bs = 're')",
           " + s(rowf, bs = 're') + s(colf, bs = 're')")
  )

  model <- gam(formula_obj, data = dfr, method = "REML")
  attr(model, "model_data") <- dfr
  model
}


#' Extract genotype BLUPs from a pooled GAMM (shared smooth, multiple benches)
#'
#' The grand mean is computed as the intercept plus the average of all bench
#' fixed effects so that BLUPs are not tied to any single reference bench.
#'
#' @param model     gamm() object from fit_gamm_pooled()
#' @param dfr       data frame used for fitting (needs gt and bench_col columns)
#' @param gt        genotype column name
#' @param bench_col bench column name
#' @return data.frame: Genotype, predicted, StdErr, reliability
extract_blups_gamm_pooled <- function(model, dfr, gt, bench_col) {
  # model is a gam() object with s(gt, bs='re') smooth
  terms_pred <- predict(model, type = "terms")
  gt_col     <- paste0("s(", gt, ")")

  geno_re    <- terms_pred[, gt_col]

  # Grand mean: fixed intercept + average bench fixed effect
  fe               <- coef(model)
  bench_idx        <- grep(paste0("^", bench_col), names(fe))
  mean_bench       <- mean(c(0, fe[bench_idx]))
  grand_mean       <- fe["(Intercept)"] + mean_bench

  # Per-genotype BLUP = grand mean + genotype random effect
  geno_factor <- as.character(dfr[[gt]])
  geno_ranef  <- tapply(geno_re, geno_factor, function(x) x[1])

  # Approximate SE and reliability from smooth variance components
  sigma2_g    <- var(geno_ranef, na.rm = TRUE)
  sigma2_e    <- var(residuals(model), na.rm = TRUE)
  n_bench     <- length(unique(dfr[[bench_col]]))
  reliability <- sigma2_g / (sigma2_g + sigma2_e / n_bench)
  blup_se     <- sqrt(sigma2_g * (1 - reliability))

  deblup <- as.numeric(grand_mean + geno_ranef / reliability)

  data.frame(
    Genotype    = names(geno_ranef),
    predicted   = as.numeric(grand_mean + geno_ranef),
    deBLUP      = deblup,
    StdErr      = blup_se,
    reliability = reliability,
    stringsAsFactors = FALSE
  )
}


#' Run full spatial GAMM pipeline with a SHARED smooth across all benches
#'
#' Fits one model pooling all benches. The te(row, col) surface is assumed
#' identical in shape across benches; bench differences are absorbed by a
#' fixed bench effect. Produces one smooth heatmap, one BLUP dotplot, and
#' one variogram across all pooled residuals.
#'
#' @param dfr            data frame (all benches combined)
#' @param pheno          response column name
#' @param gt             genotype column name
#' @param rn             row column name
#' @param cn             column column name
#' @param bench_col      bench column name
#' @param outlier_method "IQR" or "SD"
#' @param outlier_k      multiplier for outlier threshold
#' @return list(dfrF = plot-level data frame, dfrP = BLUP data frame with Bench="Pooled")
run_spatial_gamm_pooled <- function(dfr, pheno, gt, rn, cn, bench_col,
                                     outlier_method = "IQR", outlier_k = 3) {
  benches <- unique(dfr[[bench_col]])
  cat("\n--- Pooled model | Trait:", pheno,
      "| Benches:", paste(benches, collapse = ", "), "---\n")

  # k_dims from the first bench (assumes all benches share the same grid size)
  bench1 <- dfr[dfr[[bench_col]] == benches[1], ]
  k_dims <- auto_k_dims(bench1, rn, cn)
  cat("  k_dims:", k_dims[1], "x", k_dims[2], "\n")

  # Step 1: Initial pooled fit
  cat("  Fitting initial pooled GAMM...\n")
  model_init <- tryCatch(
    fit_gamm_pooled(dfr, pheno, gt, rn, cn, bench_col, k_dims),
    error = function(e) {
      cat("  ERROR in initial pooled GAMM fit:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(model_init)) return(NULL)

  # Step 2: Detect outliers across all benches jointly
  cat("  Detecting outliers (method =", outlier_method, ", k =", outlier_k, ")...\n")
  is_out <- detect_outliers_residuals(model_init, dfr, rn, cn,
                                       method = outlier_method, k = outlier_k)
  cat("[REPORT] Pooled | Trait:", pheno, "\n")
  cat("[REPORT] Method:", outlier_method, "| k:", outlier_k,
      "| Outliers:", sum(is_out), "\n")
  if (any(is_out)) {
    cat("[REPORT] Outlier bench          :", as.character(dfr[[bench_col]][is_out]), "\n")
    cat("[REPORT] Outlier row positions  :", dfr[[rn]][is_out], "\n")
    cat("[REPORT] Outlier col positions  :", dfr[[cn]][is_out], "\n")
    cat("[REPORT] Outlier observed values:", round(dfr[[pheno]][is_out], 3), "\n")
  }
  cat("[REPORT] ---\n")

  # Step 3: Replace outliers and re-fit
  dfr_clean          <- dfr
  dfr_clean[[pheno]][is_out] <- NA

  cat("  Fitting final pooled GAMM on cleaned data...\n")
  model_final <- tryCatch(
    fit_gamm_pooled(dfr_clean, pheno, gt, rn, cn, bench_col, k_dims),
    error = function(e) {
      cat("  ERROR in final pooled GAMM fit:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(model_final)) return(NULL)

  # Step 4: Extract BLUPs
  blup_df <- tryCatch(
    extract_blups_gamm_pooled(model_final, dfr_clean, gt, bench_col),
    error = function(e) {
      cat("  ERROR extracting pooled BLUPs:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(blup_df)) return(NULL)

  # Step 5: Per-plot corrections using the shared spatial surface (gam model)
  terms_mat     <- predict(model_final, type = "terms")
  te_col        <- grep("^te\\(", colnames(terms_mat), value = TRUE)[1]
  row_re_col    <- "s(rowf)"
  col_re_col    <- "s(colf)"

  spatial_raw   <- terms_mat[, te_col]
  spatial_trend <- spatial_raw - mean(spatial_raw, na.rm = TRUE)
  row_re        <- if (row_re_col %in% colnames(terms_mat)) terms_mat[, row_re_col] else 0
  col_re        <- if (col_re_col %in% colnames(terms_mat)) terms_mat[, col_re_col] else 0

  corrected <- dfr_clean[[pheno]] - spatial_trend - row_re - col_re
  fitted_v  <- fitted(model_final)

  dfrF <- data.frame(
    Genotype = dfr_clean[[gt]],
    Bench    = dfr_clean[[bench_col]],
    row      = dfr_clean[[rn]],
    col      = dfr_clean[[cn]],
    stringsAsFactors = FALSE
  )
  dfrF[[paste0(pheno, "_Observed")]]     <- dfr_clean[[pheno]]
  dfrF[[paste0(pheno, "_Fitted")]]       <- fitted_v
  dfrF[[paste0(pheno, "_Residual")]]     <- residuals(model_final)
  dfrF[[paste0(pheno, "_Corrected")]]    <- corrected
  dfrF[[paste0(pheno, "_SpatialTrend")]] <- spatial_trend
  dfrF[[paste0(pheno, "_Observed")]][is_out] <- NA

  # Step 6: BLUPs — one row per genotype, labelled "Pooled"
  blue_by_geno_p <- tapply(corrected, as.character(dfr_clean[[gt]]), mean, na.rm = TRUE)

  dfrP <- data.frame(
    Genotype = blup_df$Genotype,
    Bench    = "Pooled",
    stringsAsFactors = FALSE
  )
  dfrP[[pheno]]                         <- blup_df$predicted
  dfrP[[paste0(pheno, "_BLUE")]]        <- as.numeric(blue_by_geno_p[blup_df$Genotype])
  dfrP[[paste0(pheno, "_deBLUP")]]      <- blup_df$deBLUP
  dfrP[[paste0(pheno, "_StdErr")]]      <- blup_df$StdErr
  dfrP[[paste0(pheno, "_Reliability")]] <- blup_df$reliability

  # Step 7: Diagnostics
  cat("  Generating pooled diagnostic plots...\n")

  # 7a: Shared smooth heatmap (displayed on first-bench grid as representative)
  sh <- tryCatch(
    plot_smooth_heatmap(model_final, bench1, rn, cn,
                        bench_label = "shared across all benches"),
    error = function(e) {
      cat("  WARNING: Smooth heatmap failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(sh)) print(sh)

  # 7b: BLUP ranked dotplot
  blups_ord      <- blup_df[order(blup_df$predicted), ]
  blups_ord$rank <- seq_len(nrow(blups_ord))
  p_blup <- ggplot(blups_ord, aes(x = rank, y = predicted)) +
    geom_point(size = 1.2, colour = "#2c7bb6") +
    geom_errorbar(aes(ymin = predicted - StdErr, ymax = predicted + StdErr),
                  width = 0, colour = "#2c7bb6", alpha = 0.6) +
    labs(title = paste0("Pooled BLUPs [", pheno, "]"),
         x = "Genotype rank",
         y = paste0("Adjusted mean: ", pheno)) +
    theme_bw()
  print(p_blup)

  # 7c: Global empirical variogram of pooled residuals
  pool_resids <- residuals(model_final)
  vgm <- tryCatch(
    compute_variogram(pool_resids, dfr_clean[[rn]], dfr_clean[[cn]]),
    error = function(e) NULL
  )
  if (!is.null(vgm)) {
    p_vgm <- ggplot(vgm, aes(x = dist_bin_mid, y = semivariance)) +
      geom_point(aes(size = n_pairs), colour = "#d7301f") +
      geom_line(colour = "#d7301f", linewidth = 0.5) +
      scale_size_continuous(name = "n pairs", range = c(1, 4)) +
      labs(title = paste0("Variogram of pooled residuals [", pheno, "]"),
           x = "Distance", y = "Semivariance") +
      theme_bw()
    print(p_vgm)
  }

  list(dfrF = dfrF, dfrP = dfrP)
}


# ==============================================================================
# CONFIG — user edits this section
# ==============================================================================

fn             <- "Phenotype_data.csv"
fnpdf          <- "GAMM_plots.pdf"
report_file    <- "Outliers_GAMM_report.txt"
fnFitted       <- "GAMM_fitted_values.csv"
fnPredicted    <- "GAMM_BLUPs.csv"

gt_col         <- "Genotype"
rn_col         <- "Row"
cn_col         <- "Column"
bench_col      <- "Bench"       # set to NULL for single-field (no bench column)
trFc           <- 5             # column index of first trait
study          <- "MyStudy"
year           <- "2025"

outlier_method <- "IQR"         # "IQR" or "SD"
outlier_k      <- 3             # multiplier

bench_filter       <- NULL      # NULL = all benches; character vector to process a subset, e.g. c("Bench1","Bench2")
separate_smoothers <- TRUE      # TRUE = separate te() smooth per bench (default)
                                # FALSE = one shared te() smooth across all benches (assumes same spatial pattern everywhere)


# ==============================================================================
# MAIN EXECUTION
# Run this block directly with Rscript; skipped when sourced or knitted.
# ==============================================================================

if (!isTRUE(getOption("knitr.in.progress")) && sys.nframe() == 0) {

# ---- Read data ---------------------------------------------------------------
dfr <- read.csv(fn)
cat("Data loaded:", nrow(dfr), "rows,", ncol(dfr), "columns\n")
cat("Columns:", paste(colnames(dfr), collapse = ", "), "\n")

# Add unique plot identifier (genotype + position)
dfr$genotype_pos <- paste(dfr[[gt_col]], dfr[[rn_col]], dfr[[cn_col]], sep = "_")

# Identify traits (columns trFc to second-to-last; last is genotype_pos)
trLc   <- ncol(dfr) - 1
traits <- colnames(dfr)[trFc:trLc]
cat("Traits:", paste(traits, collapse = ", "), "\n")

# Identify benches
if (!is.null(bench_col) && bench_col %in% colnames(dfr)) {
  benches <- unique(dfr[[bench_col]])
} else {
  benches    <- "AllData"
  bench_col  <- NULL
}

# Apply bench filter
if (!is.null(bench_filter) && !is.null(bench_col)) {
  benches <- intersect(benches, bench_filter)
  if (length(benches) == 0) stop("bench_filter matches no benches found in the data.")
}
cat("Benches:", paste(benches, collapse = ", "), "\n")
cat("Mode:   ", if (is.null(bench_col) || !separate_smoothers) "shared smooth (pooled)"
               else "separate smooth per bench", "\n")

# ---- Open output devices ------------------------------------------------------
pdf(fnpdf, width = 14, height = 8)
report_con <- file(report_file, open = "w")
sink(report_con, type = "output")
sink(report_con, type = "message")

cat("=== GAMM Spatial Correction Report ===\n")
cat("Study:", study, "| Year:", year, "\n")
cat("Mode:", if (!is.null(bench_col) && separate_smoothers)
              "separate spatial smooth per bench"
            else
              "shared spatial smooth across all benches", "\n")
cat("Outlier method:", outlier_method, "| k:", outlier_k, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- Accumulate helper -------------------------------------------------------
accumulate_fitted <- function(combined, res_dfrF, tr) {
  if (is.null(combined)) return(res_dfrF)
  obs_col <- paste0(tr, "_Observed");  fit_col <- paste0(tr, "_Fitted")
  res_col <- paste0(tr, "_Residual"); cor_col <- paste0(tr, "_Corrected")
  spt_col <- paste0(tr, "_SpatialTrend")
  keep <- intersect(c("Genotype", "Bench", "row", "col",
                       obs_col, fit_col, res_col, cor_col, spt_col),
                    colnames(res_dfrF))
  merge(combined, res_dfrF[, keep, drop = FALSE],
        by = c("Genotype", "Bench", "row", "col"), all = TRUE)
}

# ---- Main loop ---------------------------------------------------------------
combined_fitted_data <- NULL
blup_list            <- list()

if (!is.null(bench_col) && separate_smoothers) {
  # ---- Mode A: separate spatial smooth per bench (default) ------------------
  for (bench in benches) {
    bench_dfr <- dfr[dfr[[bench_col]] == bench, ]

    for (tr in traits) {
      keep_cols <- unique(c(tr, gt_col, rn_col, cn_col, bench_col))
      mydfr     <- bench_dfr[, keep_cols, drop = FALSE]
      mydfr[[tr]] <- as.numeric(mydfr[[tr]])
      pmain <- paste(study, year, sub("_.*", "", tr), bench)

      res <- tryCatch(
        run_spatial_gamm(
          mydfr = mydfr, pheno = tr, gt = gt_col, rn = rn_col, cn = cn_col,
          bench_label = bench, pmain = pmain,
          outlier_method = outlier_method, outlier_k = outlier_k
        ),
        error = function(e) {
          cat("ERROR [bench", bench, "| trait", tr, "]:", conditionMessage(e), "\n")
          NULL
        }
      )
      if (is.null(res)) next
      combined_fitted_data <- accumulate_fitted(combined_fitted_data, res$dfrF, tr)
      blup_list <- c(blup_list, list(res$dfrP))
    }
  }

} else {
  # ---- Mode B: shared smooth across all benches (or single field) -----------
  for (tr in traits) {
    keep_cols <- unique(c(tr, gt_col, rn_col, cn_col,
                          if (!is.null(bench_col)) bench_col))
    mydfr_all <- dfr[if (!is.null(bench_col)) dfr[[bench_col]] %in% benches
                     else rep(TRUE, nrow(dfr)), keep_cols, drop = FALSE]
    mydfr_all[[tr]] <- as.numeric(mydfr_all[[tr]])

    if (is.null(bench_col)) {
      # Single field — no bench structure
      res <- tryCatch(
        run_spatial_gamm(
          mydfr = mydfr_all, pheno = tr, gt = gt_col, rn = rn_col, cn = cn_col,
          bench_label = "AllData", pmain = paste(study, year, tr),
          outlier_method = outlier_method, outlier_k = outlier_k
        ),
        error = function(e) {
          cat("ERROR [trait", tr, "]:", conditionMessage(e), "\n")
          NULL
        }
      )
    } else {
      # Pooled multi-bench model
      res <- tryCatch(
        run_spatial_gamm_pooled(
          dfr = mydfr_all, pheno = tr, gt = gt_col, rn = rn_col, cn = cn_col,
          bench_col = bench_col,
          outlier_method = outlier_method, outlier_k = outlier_k
        ),
        error = function(e) {
          cat("ERROR [pooled | trait", tr, "]:", conditionMessage(e), "\n")
          NULL
        }
      )
    }
    if (is.null(res)) next
    combined_fitted_data <- accumulate_fitted(combined_fitted_data, res$dfrF, tr)
    blup_list <- c(blup_list, list(res$dfrP))
  }
}

# ---- Close output devices ----------------------------------------------------
sink(type = "message")
sink(type = "output")
close(report_con)
dev.off()

cat("PDF written to:", fnpdf, "\n")
cat("Report written to:", report_file, "\n")

# ==============================================================================
# POST-PROCESSING
# ==============================================================================

# ---- Build merged fitted data ------------------------------------------------
combined_fitted_data$genotype_pos <- paste(
  combined_fitted_data$Genotype,
  combined_fitted_data$row,
  combined_fitted_data$col,
  sep = "_"
)

# Merge back any extra original columns (e.g. Replicate) by genotype_pos
extra_orig_cols <- setdiff(
  colnames(dfr),
  c(traits, gt_col, rn_col, cn_col,
    if (!is.null(bench_col)) bench_col,
    "genotype_pos")
)
if (length(extra_orig_cols) > 0) {
  subdfr <- dfr[, c("genotype_pos", extra_orig_cols), drop = FALSE]
  merged_fitted_data <- merge(subdfr, combined_fitted_data,
                               by = "genotype_pos", all = TRUE)
} else {
  merged_fitted_data <- combined_fitted_data
}

# Remove helper column and move Genotype first
merged_fitted_data <- merged_fitted_data[,
  !names(merged_fitted_data) %in% "genotype_pos"]
merged_fitted_data <- merged_fitted_data[,
  c("Genotype", setdiff(names(merged_fitted_data), "Genotype"))]

write.csv(merged_fitted_data, fnFitted, row.names = FALSE, na = "")
cat("Fitted values written to:", fnFitted, "\n")

# ---- Build per-bench BLUP table from list ------------------------------------
# Merge all per-(bench x trait) dfrP data frames:
# Within each bench, merge traits by Genotype + Bench; then rbind benches.
if (length(blup_list) > 0) {
  # Split by bench
  bench_labels <- sapply(blup_list, function(x) x$Bench[1])
  blup_by_bench <- split(blup_list, bench_labels)

  newdfr <- do.call(rbind, lapply(blup_by_bench, function(bench_blups) {
    Reduce(
      function(a, b) merge(a, b, by = c("Genotype", "Bench"), all = TRUE),
      bench_blups
    )
  }))
  rownames(newdfr) <- NULL
} else {
  newdfr <- NULL
}

# ---- Build per-bench BLUP table and cross-bench summary ---------------------
if (!is.null(newdfr)) {
  blup_traits <- traits

  # Round BLUP and StdErr columns to 3 significant figures
  round_cols <- c(blup_traits,
                  paste0(blup_traits, "_StdErr"),
                  paste0(blup_traits, "_Reliability"))
  round_cols <- intersect(round_cols, colnames(newdfr))
  newdfr[, round_cols] <- lapply(newdfr[, round_cols, drop = FALSE],
                                  function(x) signif(as.numeric(x), 3))

  is_pooled <- !is.null(bench_col) && !separate_smoothers

  if (is_pooled) {
    # Pooled model: BLUPs are already cross-bench (Bench = "Pooled")
    # No per-bench breakdown to average; just output as-is
    finaldfrP <- newdfr
  } else {
    # Separate smoothers: add cross-bench mean rows
    numeric_cols <- intersect(round_cols, colnames(newdfr))
    mean_rows <- newdfr %>%
      group_by(Genotype) %>%
      summarise(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(Bench = "Mean")
    finaldfrP <- bind_rows(newdfr, mean_rows)
  }

  # Ensure Genotype and Bench are first two columns
  other_cols <- setdiff(colnames(finaldfrP), c("Genotype", "Bench"))
  finaldfrP  <- finaldfrP[, c("Genotype", "Bench", other_cols)]

  write.csv(finaldfrP, fnPredicted, row.names = FALSE, na = "")
  cat("BLUPs written to:", fnPredicted, "\n")

  # Summary
  cat("\n=== Summary ===\n")
  cat("Genotypes:", length(unique(finaldfrP$Genotype)), "\n")
  cat("Benches:", paste(unique(finaldfrP$Bench), collapse = ", "), "\n")
  cat("Traits:", paste(blup_traits, collapse = ", "), "\n")

  print(head(finaldfrP))
}

cat("\nDone.\n")

# ---- Render all Rmd files in docs/ -------------------------------------------
docs_dir <- normalizePath(file.path(dirname(fn), "..", "docs"), mustWork = FALSE)
if (!dir.exists(docs_dir)) docs_dir <- normalizePath("docs", mustWork = FALSE)

if (requireNamespace("rmarkdown", quietly = TRUE) && dir.exists(docs_dir)) {
  rmd_files <- list.files(docs_dir, pattern = "[.]Rmd$", full.names = TRUE)
  if (length(rmd_files) == 0) {
    cat("(no Rmd files found in", docs_dir, "— skipping render)\n")
  } else {
    for (rmd in rmd_files) {
      rmd_abs <- normalizePath(rmd, mustWork = FALSE)
      cat("\nRendering:", rmd_abs, "\n")
      tryCatch(
        rmarkdown::render(
          input         = rmd_abs,
          output_dir    = dirname(rmd_abs),
          knit_root_dir = dirname(rmd_abs),
          quiet         = TRUE
        ),
        error = function(e) cat("  WARNING: render failed:", conditionMessage(e), "\n")
      )
      html_out <- sub("[.]Rmd$", ".html", rmd_abs)
      cat("  -> ", html_out, "\n")
    }
  }
} else {
  cat("(rmarkdown not available or docs/ not found — skipping render)\n")
}

} # end if (run directly)
