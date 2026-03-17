# ============================================================================
# spatial_correct_gam.R
# mgcv-based spatial correction for field trial data
#
# Best-performing variant: te_ps_re
#   - 2D P-spline tensor product spatial surface  [te(..., bs='ps')]
#   - Row and column random effects               [s(row_f, bs='re') etc.]
#   - BLUEs as primary output (genotype fixed)
#
# Single-bench:  per-bench gam() with te(row, col, by=NULL)
# Multi-bench:   joint gam() with te(row, col, by=bench_f) + bench fixed effect
#
# Usage:
#   source("scripts/spatial_correct_gam.R")
#   out <- run_spatial_gam(
#     data_file  = "data/mydata.csv",
#     pheno_cols = c("yield"),
#     geno_col   = "Genotype",
#     row_col    = "Row",
#     col_col    = "Col",
#     bench_col  = "Bench",         # NULL for single-bench
#     output_dir = "output/gam"
#   )
#
# Outputs written to output_dir:
#   BLUEs.csv                        - one row per genotype, one col per trait
#   BLUPs.csv                        - same for BLUPs (if requested)
#   spatial_trends.csv               - per-plot spatial effect
#   model_summary.csv                - per-trait/bench fit statistics
#   diagnostics_{pheno}_{bench}.png  - 6-panel diagnostic plots
#   blues_distribution.png           - BLUE density per trait
#   spatial_surfaces.png             - spatial trend panels
#
# Returns (invisibly): list(blues, blups, spatial_trends, model_summary)
# ============================================================================


# ============================================================================
# Section 1: Packages
# ============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(ggplot2)
  library(patchwork)
})


# ============================================================================
# Section 2: Helper functions
# ============================================================================

#' Choose number of B-spline basis functions adaptively
#' Rule: roughly half the unique positions, bounded to [5, 20].
#' @param n_unique  Number of unique positions along one axis
#' @return integer
adaptive_nseg <- function(n_unique) {
  as.integer(min(max(5L, floor(n_unique / 2L)), 20L))
}


#' Read CSV or R binary (.rda/.RData) file
#' @param fn          File path
#' @param rda_object  Object name when fn is .rda; NULL = first data frame
#' @return data.frame
read_input <- function(fn, rda_object = NULL) {
  if (!file.exists(fn)) stop("File not found: ", fn)
  ext <- tolower(tools::file_ext(fn))
  if (ext == "csv") {
    return(read.csv(fn, stringsAsFactors = FALSE))
  }
  if (ext %in% c("rda", "rdata")) {
    e <- new.env(parent = emptyenv())
    load(fn, envir = e)
    if (!is.null(rda_object)) {
      if (!exists(rda_object, envir = e))
        stop("Object '", rda_object, "' not found in ", fn)
      return(as.data.frame(get(rda_object, envir = e)))
    }
    objs <- ls(e)
    dfs  <- Filter(function(nm) is.data.frame(get(nm, envir = e)), objs)
    if (length(dfs) == 0) stop("No data frames found in ", fn)
    if (length(dfs) > 1)
      message("Multiple data frames in ", fn, "; using '", dfs[1],
              "'. Set rda_object= to choose.")
    return(as.data.frame(get(dfs[1], envir = e)))
  }
  stop("Unsupported file type '", ext, "'. Use .csv, .rda, or .RData.")
}


#' Compute empirical semivariogram
#' @param data      data.frame with row, col, and value columns
#' @param row_col   Row coordinate column name
#' @param col_col   Column coordinate column name
#' @param value_col Value column name
#' @param n_bins    Number of distance bins
#' @return data.frame(dist, sv) or NULL when < 10 non-NA observations
empirical_semivariogram <- function(data, row_col, col_col, value_col,
                                    n_bins = 15) {
  ok  <- !is.na(data[[value_col]])
  d   <- data[ok, ]
  if (nrow(d) < 10) return(NULL)

  rows <- d[[row_col]]
  cols <- d[[col_col]]
  vals <- d[[value_col]]

  row_diff <- outer(rows, rows, "-")
  col_diff <- outer(cols, cols, "-")
  dists    <- sqrt(row_diff^2 + col_diff^2)
  sv_mat   <- (outer(vals, vals, "-"))^2 / 2

  up     <- upper.tri(dists)
  d_vec  <- dists[up]
  sv_vec <- sv_mat[up]

  max_d  <- quantile(d_vec, 0.6)
  keep   <- d_vec <= max_d & d_vec > 0
  d_vec  <- d_vec[keep]
  sv_vec <- sv_vec[keep]
  if (length(d_vec) == 0) return(NULL)

  breaks <- seq(0, max(d_vec), length.out = n_bins + 1)
  bin    <- cut(d_vec, breaks = breaks, include.lowest = TRUE)
  agg    <- aggregate(cbind(sv = sv_vec, dist = d_vec) ~ bin, FUN = mean)
  agg[order(agg$dist), c("dist", "sv")]
}


#' Plot a field heatmap
#' @param data      data.frame with coordinates and fill column
#' @param row_col   Row coordinate column name
#' @param col_col   Column coordinate column name
#' @param value_col Fill column name
#' @param title     Plot title
#' @param limits    c(min, max) for fill scale; NULL = auto
#' @return ggplot object
plot_heatmap <- function(data, row_col, col_col, value_col,
                          title = "", limits = NULL) {
  ok_rows <- !is.na(data[[row_col]]) & !is.na(data[[col_col]])
  d       <- data[ok_rows, ]
  if (is.null(limits)) limits <- range(d[[value_col]], na.rm = TRUE)

  ggplot(d, aes(x = .data[[col_col]], y = .data[[row_col]],
                fill = .data[[value_col]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(
      name     = value_col,
      limits   = limits,
      option   = "viridis",
      na.value = "grey70"
    ) +
    scale_y_reverse(breaks = sort(unique(d[[row_col]]))) +
    scale_x_continuous(breaks = sort(unique(d[[col_col]]))) +
    labs(title = title, x = col_col, y = row_col) +
    theme_bw(base_size = 8) +
    theme(
      plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid   = element_blank(),
      axis.text    = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 6)
    )
}


#' Plot empirical semivariogram
#' @param vgram  data.frame(dist, sv) from empirical_semivariogram()
#' @param title  Plot title
#' @return ggplot object
plot_variogram <- function(vgram, title = "") {
  if (is.null(vgram) || nrow(vgram) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
                 size = 3, colour = "grey50") +
        theme_void() +
        labs(title = title) +
        theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
  }
  ggplot(vgram, aes(x = dist, y = sv)) +
    geom_point(colour = "#2c7bb6", size = 1.5) +
    geom_line(colour = "#2c7bb6", linewidth = 0.5) +
    labs(title = title, x = "Distance (plot units)", y = "Semivariance") +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6)
    )
}


#' Plot row and column random effect magnitudes
#'
#' Returns a two-panel (row | col) strip plot showing the estimated random
#' effect for each row and column level. Useful for diagnosing systematic
#' edge or lane effects that the smooth surface cannot capture.
#'
#' @param bench_data  data.frame for the bench
#' @param result      Return value from fit_mgcv_bench() or fit_mgcv_joint()
#' @param row_col     Row coordinate column name
#' @param col_col     Column coordinate column name
#' @param bench_label Character label for title
#' @return patchwork ggplot or NULL if no RE data available
plot_row_col_re <- function(bench_data, result, row_col, col_col,
                            bench_label = "") {
  if (is.null(result$row_re) && is.null(result$col_re)) return(NULL)

  plots <- list()

  if (!is.null(result$row_re)) {
    rd <- data.frame(level = bench_data[[row_col]], effect = result$row_re)
    rd <- aggregate(effect ~ level, data = rd, FUN = mean)
    plots$row <- ggplot(rd, aes(x = factor(level), y = effect)) +
      geom_col(fill = "#2c7bb6", alpha = 0.7, width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(title = "Row random effects", x = row_col, y = "Effect") +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 7),
            axis.text  = element_text(size = 6))
  }

  if (!is.null(result$col_re)) {
    cd <- data.frame(level = bench_data[[col_col]], effect = result$col_re)
    cd <- aggregate(effect ~ level, data = cd, FUN = mean)
    plots$col <- ggplot(cd, aes(x = factor(level), y = effect)) +
      geom_col(fill = "#d7191c", alpha = 0.7, width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(title = "Column random effects", x = col_col, y = "Effect") +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 7),
            axis.text  = element_text(size = 6))
  }

  if (length(plots) == 0) return(NULL)

  page_title <- paste0("Row / column random effects",
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  composite <- wrap_plots(plots, ncol = length(plots)) +
    plot_annotation(
      title = page_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
  composite
}


#' Compose 6-panel diagnostic plot for one bench x trait
#'
#' Panels: observed | spatial trend | residuals (top row)
#'         variogram raw | variogram residuals | QQ residuals (bottom row)
#'
#' @param bench_data  data.frame for the bench (includes observed pheno)
#' @param result      Return value from fit_mgcv_bench() or fit_mgcv_joint()
#' @param pheno       Phenotype column name in bench_data
#' @param row_col     Row coordinate column name
#' @param col_col     Column coordinate column name
#' @param bench_label Character label for title
#' @return patchwork ggplot
make_diagnostic_plots <- function(bench_data, result, pheno,
                                   row_col, col_col, bench_label = "") {
  pd <- bench_data[, c(row_col, col_col), drop = FALSE]
  pd$obs     <- bench_data[[pheno]]
  pd$spatial <- if (!is.null(result$spatial_smooth)) result$spatial_smooth
                else if (!is.null(result$spatial))   result$spatial
                else NA_real_
  pd$resid   <- if (!is.null(result$residuals)) result$residuals else NA_real_

  safe_range <- function(...) {
    r <- range(..., na.rm = TRUE)
    if (any(!is.finite(r))) c(0, 1) else r
  }

  p1 <- plot_heatmap(pd, row_col, col_col, "obs",
                     title = "Observed", limits = safe_range(pd$obs))
  p2 <- plot_heatmap(pd, row_col, col_col, "spatial",
                     title = "Spatial trend", limits = safe_range(pd$spatial))
  p3 <- plot_heatmap(pd, row_col, col_col, "resid",
                     title = "Residuals", limits = safe_range(pd$resid))

  vg_obs  <- empirical_semivariogram(pd, row_col, col_col, "obs")
  vg_res  <- empirical_semivariogram(pd, row_col, col_col, "resid")
  p4 <- plot_variogram(vg_obs, title = "Variogram: observed")
  p5 <- plot_variogram(vg_res, title = "Variogram: residuals")

  # QQ plot of corrected residuals
  qq_df <- data.frame(resid = pd$resid[!is.na(pd$resid)])
  p6 <- if (nrow(qq_df) > 2) {
    ggplot(qq_df, aes(sample = resid)) +
      stat_qq(size = 1, colour = "#2c7bb6") +
      stat_qq_line(colour = "red", linewidth = 0.5) +
      labs(title = "QQ: corrected residuals", x = "Theoretical", y = "Sample") +
      theme_bw(base_size = 8) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text  = element_text(size = 6)
      )
  } else {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
               size = 3, colour = "grey50") +
      theme_void() +
      labs(title = "QQ: corrected residuals") +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  }

  page_title <- paste0(pheno,
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_annotation(
      title = page_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
}


# ============================================================================
# Section 3: Core fitting functions
# ============================================================================

#' mgcv te_ps_re spatial correction for a single bench
#'
#' Formula (BLUEs):
#'   pheno ~ 0 + Genotype + te(row, col, bs='ps') + s(row_f, bs='re') + s(col_f, bs='re')
#'
#' @param bench_data    data.frame for one bench
#' @param pheno         Phenotype column name
#' @param geno_col      Genotype column name
#' @param row_col       Row coordinate column name
#' @param col_col       Column coordinate column name
#' @param k_row         Basis dimension for rows (NULL = auto)
#' @param k_col         Basis dimension for columns (NULL = auto)
#' @param estimate_type "BLUEs", "BLUPs", or "both"
#' @return list(blues, blups, fitted, residuals, spatial_smooth, spatial_total,
#'              row_re, col_re, converged, edf_spatial)
fit_mgcv_bench <- function(bench_data, pheno, geno_col, row_col, col_col,
                            k_row = NULL, k_col = NULL,
                            estimate_type = "BLUEs") {
  res <- list(blues = NULL, blups = NULL, fitted = NULL,
              residuals = NULL, spatial_smooth = NULL, spatial_total = NULL,
              row_re = NULL, col_re = NULL,
              converged = NA, edf_spatial = NA)

  bench_data[[geno_col]] <- as.factor(bench_data[[geno_col]])
  bench_data$row_f       <- as.factor(bench_data[[row_col]])
  bench_data$col_f       <- as.factor(bench_data[[col_col]])
  geno_levels            <- levels(bench_data[[geno_col]])

  n_row <- length(unique(bench_data[[row_col]]))
  n_col <- length(unique(bench_data[[col_col]]))
  kr    <- if (!is.null(k_row)) k_row else adaptive_nseg(n_row)
  kc    <- if (!is.null(k_col)) k_col else adaptive_nseg(n_col)

  te_term  <- paste0("te(", row_col, ", ", col_col,
                     ", bs=c('ps','ps'), k=c(", kr, ",", kc, "))")
  re_extra <- " + s(row_f, bs='re') + s(col_f, bs='re')"

  # ---- BLUEs ---------------------------------------------------------------
  if (estimate_type %in% c("BLUEs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ 0 + ", geno_col, " + ", te_term, re_extra))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { warning("mgcv BLUEs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      res$converged <- m$converged
      if (!isTRUE(m$converged))
        warning("GAM did not converge for BLUEs (", pheno, ")")

      coef_names <- paste0(geno_col, geno_levels)
      blues_vals <- coef(m)[coef_names]
      vcov_diag  <- diag(vcov(m))
      blues_se   <- sqrt(vcov_diag[coef_names])

      res$blues <- data.frame(
        Genotype = geno_levels,
        BLUE     = as.numeric(blues_vals),
        BLUE_SE  = as.numeric(blues_se),
        stringsAsFactors = FALSE
      )
      names(res$blues)[1] <- geno_col

      terms_pred     <- predict(m, type = "terms")
      te_col         <- grep("^te\\(", colnames(terms_pred), value = TRUE)
      row_re_col     <- grep("^s\\(row_f\\)", colnames(terms_pred), value = TRUE)
      col_re_col     <- grep("^s\\(col_f\\)", colnames(terms_pred), value = TRUE)

      smooth_val     <- as.numeric(terms_pred[, te_col[1]])
      row_re_val     <- if (length(row_re_col)) as.numeric(terms_pred[, row_re_col[1]]) else rep(0, nrow(bench_data))
      col_re_val     <- if (length(col_re_col)) as.numeric(terms_pred[, col_re_col[1]]) else rep(0, nrow(bench_data))

      res$spatial_smooth <- smooth_val
      res$spatial_total  <- smooth_val + row_re_val + col_re_val
      res$row_re         <- row_re_val
      res$col_re         <- col_re_val
      res$fitted         <- as.numeric(fitted(m))
      res$residuals      <- as.numeric(residuals(m))
      res$edf_spatial    <- sum(m$edf[grep("^te\\(", names(m$edf))])
    }
  }

  # ---- BLUPs ---------------------------------------------------------------
  if (estimate_type %in% c("BLUPs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ s(", geno_col, ", bs='re') + ",
                            te_term, re_extra))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { warning("mgcv BLUPs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      if (is.na(res$converged)) res$converged <- m$converged
      if (!isTRUE(m$converged))
        warning("GAM did not converge for BLUPs (", pheno, ")")

      re_pattern <- paste0("s\\(", geno_col, "\\)")
      re_idx     <- grep(re_pattern, names(coef(m)))
      blup_vals  <- coef(m)[re_idx]
      intercept  <- coef(m)["(Intercept)"]
      blup_total <- as.numeric(intercept) + as.numeric(blup_vals)

      # SE via prediction at mean spatial location for each genotype
      mean_row <- mean(bench_data[[row_col]], na.rm = TRUE)
      mean_col <- mean(bench_data[[col_col]], na.rm = TRUE)
      newdat   <- setNames(
        data.frame(geno_levels, mean_row, mean_col),
        c(geno_col, row_col, col_col)
      )
      newdat[[geno_col]] <- factor(newdat[[geno_col]], levels = geno_levels)
      newdat$row_f <- factor(round(mean_row), levels = levels(bench_data$row_f))
      newdat$col_f <- factor(round(mean_col), levels = levels(bench_data$col_f))
      tp <- tryCatch(
        predict(m, newdata = newdat, type = "terms", se.fit = TRUE),
        error = function(e) NULL
      )
      re_col_name <- if (!is.null(tp))
        grep(re_pattern, colnames(tp$fit), value = TRUE)[1] else NA
      blup_se <- if (!is.null(tp) && !is.na(re_col_name))
        as.numeric(tp$se.fit[, re_col_name]) else rep(NA_real_, length(geno_levels))

      res$blups <- data.frame(
        Genotype = geno_levels,
        BLUP     = blup_total,
        BLUP_SE  = blup_se,
        stringsAsFactors = FALSE
      )
      names(res$blups)[1] <- geno_col

      if (is.null(res$fitted)) {
        terms_pred    <- predict(m, type = "terms")
        te_col        <- grep("^te\\(", colnames(terms_pred), value = TRUE)
        res$spatial   <- as.numeric(terms_pred[, te_col[1]])
        res$fitted    <- as.numeric(fitted(m))
        res$residuals <- as.numeric(residuals(m))
        if (is.na(res$edf_spatial))
          res$edf_spatial <- sum(m$edf[grep("^te\\(", names(m$edf))])
      }
    }
  }

  res
}


#' Joint multi-bench mgcv::gam spatial correction
#'
#' Fits a single GAM across all benches:
#'   pheno ~ 0 + Genotype + bench_f +
#'           te(row, col, bs='ps', by=bench_f) +
#'           s(row_f, bench_f, bs='re') + s(col_f, bench_f, bs='re')
#'
#' @param data          Full data.frame (all benches)
#' @param pheno         Phenotype column name
#' @param geno_col      Genotype column name
#' @param row_col       Row coordinate column name
#' @param col_col       Column coordinate column name
#' @param bench_col     Bench column name
#' @param k_row         Basis dimension for rows (NULL = auto)
#' @param k_col         Basis dimension for columns (NULL = auto)
#' @param estimate_type Currently only "BLUEs" (multi-bench BLUPs not implemented)
#' @return list(blues, fitted, residuals, spatial, converged, edf_spatial)
fit_mgcv_joint <- function(data, pheno, geno_col, row_col, col_col, bench_col,
                            k_row = NULL, k_col = NULL,
                            estimate_type = "BLUEs") {
  res <- list(blues = NULL, fitted = NULL, residuals = NULL,
              spatial_smooth = NULL, spatial_total = NULL,
              row_re = NULL, col_re = NULL,
              converged = NA, edf_spatial = NA)

  data[[geno_col]] <- as.factor(data[[geno_col]])
  data$bench_f     <- as.factor(data[[bench_col]])
  data$row_f       <- as.factor(data[[row_col]])
  data$col_f       <- as.factor(data[[col_col]])
  geno_levels      <- levels(data[[geno_col]])
  bench_levels     <- levels(data$bench_f)

  # Adaptive k: use minimum dimensions across benches
  nr_vec <- sapply(bench_levels, function(b)
    length(unique(data[[row_col]][data$bench_f == b])))
  nc_vec <- sapply(bench_levels, function(b)
    length(unique(data[[col_col]][data$bench_f == b])))
  kr <- if (!is.null(k_row)) k_row else adaptive_nseg(min(nr_vec))
  kc <- if (!is.null(k_col)) k_col else adaptive_nseg(min(nc_vec))

  fm <- as.formula(paste0(
    pheno, " ~ 0 + ", geno_col, " + bench_f + ",
    "te(", row_col, ", ", col_col,
    ", bs=c('ps','ps'), k=c(", kr, ",", kc, "), by=bench_f) + ",
    "s(row_f, bench_f, bs='re') + s(col_f, bench_f, bs='re')"
  ))

  m <- tryCatch(
    gam(fm, data = data, method = "REML"),
    error = function(e) { warning("mgcv joint error: ", e$message); NULL }
  )
  if (is.null(m)) return(res)

  res$converged <- m$converged
  if (!isTRUE(m$converged)) warning("Joint GAM did not converge (", pheno, ")")

  # BLUEs
  coef_names  <- paste0(geno_col, geno_levels)
  blues_vals  <- coef(m)[coef_names]
  vcov_diag   <- diag(vcov(m))
  blues_se    <- sqrt(vcov_diag[coef_names])

  res$blues <- data.frame(
    Genotype = geno_levels,
    BLUE     = as.numeric(blues_vals),
    BLUE_SE  = as.numeric(blues_se),
    stringsAsFactors = FALSE
  )
  names(res$blues)[1] <- geno_col

  res$fitted    <- as.numeric(fitted(m))
  res$residuals <- as.numeric(residuals(m))

  terms_pred      <- predict(m, type = "terms")
  te_cols         <- grep("^te\\(", colnames(terms_pred), value = TRUE)
  row_re_cols     <- grep("^s\\(row_f", colnames(terms_pred), value = TRUE)
  col_re_cols     <- grep("^s\\(col_f", colnames(terms_pred), value = TRUE)

  smooth_val <- if (length(te_cols) > 0)
    as.numeric(rowSums(terms_pred[, te_cols, drop = FALSE]))
  else
    rep(0, nrow(data))
  row_re_val <- if (length(row_re_cols) > 0)
    as.numeric(rowSums(terms_pred[, row_re_cols, drop = FALSE]))
  else
    rep(0, nrow(data))
  col_re_val <- if (length(col_re_cols) > 0)
    as.numeric(rowSums(terms_pred[, col_re_cols, drop = FALSE]))
  else
    rep(0, nrow(data))

  res$spatial_smooth <- smooth_val
  res$spatial_total  <- smooth_val + row_re_val + col_re_val
  res$row_re         <- row_re_val
  res$col_re         <- col_re_val
  res$edf_spatial    <- sum(m$edf[grep("^te\\(", names(m$edf))])

  res
}


# ============================================================================
# Section 4: Main orchestrator function
# ============================================================================

#' Run mgcv-based spatial correction end to end
#'
#' Reads data, validates columns, fits the spatial model (single- or
#' multi-bench), writes CSV outputs and diagnostic plots to output_dir,
#' and returns the key results invisibly.
#'
#' @param data_file     Path to input file (.csv, .rda, or .RData)
#' @param pheno_cols    Character vector of phenotype column names
#' @param output_dir    Directory for all outputs (created if absent)
#' @param geno_col      Genotype identifier column name
#' @param row_col       Row coordinate column name (numeric)
#' @param col_col       Column coordinate column name (numeric)
#' @param bench_col     Bench/trial column name; NULL triggers single-bench model
#' @param estimate_type One of "BLUEs" (default), "BLUPs", or "both"
#' @param k_row         Basis dimension for row axis; NULL = auto
#' @param k_col         Basis dimension for column axis; NULL = auto
#' @param rda_object    Object name to extract from .rda file; NULL = first data frame
#'
#' @return Invisibly: list with elements blues, blups, spatial_trends, model_summary
run_spatial_gam <- function(
  data_file,
  pheno_cols,
  output_dir    = "output/gam",
  geno_col      = "geno",
  row_col       = "row",
  col_col       = "col",
  bench_col     = NULL,
  estimate_type = "BLUEs",
  k_row         = NULL,
  k_col         = NULL,
  rda_object    = NULL
) {

  # --------------------------------------------------------------------------
  # Data loading & validation
  # --------------------------------------------------------------------------
  cat("-- Loading data ----------------------------------------------------------\n")
  data <- read_input(data_file, rda_object = rda_object)
  cat("  File:", data_file, "\n")
  cat("  Dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

  required_cols <- c(geno_col, row_col, col_col, pheno_cols)
  if (!is.null(bench_col)) required_cols <- c(required_cols, bench_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("Required column(s) not found in data: ",
         paste(missing_cols, collapse = ", "), "\n",
         "  Available columns: ", paste(names(data), collapse = ", "))

  for (coord in c(row_col, col_col)) {
    if (!is.numeric(data[[coord]]))
      stop("Column '", coord, "' must be numeric (found: ",
           class(data[[coord]]), ")")
    if (anyNA(data[[coord]]))
      stop("Column '", coord, "' contains NA values; remove or impute before running.")
  }

  for (ph in pheno_cols) {
    if (!is.numeric(data[[ph]]))
      stop("Phenotype column '", ph, "' must be numeric.")
  }

  n_unique_row <- length(unique(data[[row_col]]))
  n_unique_col <- length(unique(data[[col_col]]))
  if (n_unique_row < 5)
    warning("Only ", n_unique_row, " unique row values -- spatial surface may not be identifiable.")
  if (n_unique_col < 5)
    warning("Only ", n_unique_col, " unique column values -- spatial surface may not be identifiable.")

  is_multibench <- !is.null(bench_col)
  if (is_multibench) {
    bench_ids <- sort(unique(data[[bench_col]]))
    n_bench   <- length(bench_ids)
    cat("  Design: multi-bench (", n_bench, "benches )\n")
    cat("  Benches:", paste(bench_ids, collapse = ", "), "\n")
  } else {
    n_bench   <- 1L
    bench_ids <- "single"
    cat("  Design: single bench\n")
  }

  n_geno <- length(unique(data[[geno_col]]))
  cat("  Genotypes:", n_geno, "\n")
  cat("  Grid:", n_unique_row, "rows x", n_unique_col, "cols\n")
  cat("  Traits:", paste(pheno_cols, collapse = ", "), "\n")
  cat("  Estimate type:", estimate_type, "\n")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  Created output directory:", output_dir, "\n")
  }

  # --------------------------------------------------------------------------
  # Model fitting
  # --------------------------------------------------------------------------
  cat("\n-- Fitting models --------------------------------------------------------\n")

  all_blues_list <- list()
  all_blups_list <- list()
  spatial_rows   <- list()
  summary_rows   <- list()

  for (pheno in pheno_cols) {
    cat("  Trait:", pheno, "\n")

    if (!is_multibench) {
      # Single-bench path
      result <- fit_mgcv_bench(
        bench_data    = data,
        pheno         = pheno,
        geno_col      = geno_col,
        row_col       = row_col,
        col_col       = col_col,
        k_row         = k_row,
        k_col         = k_col,
        estimate_type = estimate_type
      )

      if (!is.null(result$blues)) {
        df <- result$blues
        names(df)[names(df) == "BLUE"]    <- pheno
        names(df)[names(df) == "BLUE_SE"] <- paste0(pheno, "_SE")
        all_blues_list[[pheno]] <- df
      }
      if (!is.null(result$blups)) {
        df <- result$blups
        names(df)[names(df) == "BLUP"]    <- pheno
        names(df)[names(df) == "BLUP_SE"] <- paste0(pheno, "_SE")
        all_blups_list[[pheno]] <- df
      }

      if (!is.null(result$spatial_smooth)) {
        spatial_rows[[length(spatial_rows) + 1]] <- data.frame(
          data[, c(geno_col, row_col, col_col), drop = FALSE],
          bench          = "single",
          pheno          = pheno,
          spatial_smooth = result$spatial_smooth,
          spatial_total  = result$spatial_total,
          stringsAsFactors = FALSE
        )
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        pheno       = pheno,
        bench       = "single",
        resid_sd    = if (!is.null(result$residuals)) sd(result$residuals, na.rm = TRUE) else NA,
        edf_spatial = result$edf_spatial,
        converged   = result$converged,
        stringsAsFactors = FALSE
      )

      cat("    Saving diagnostics...\n")
      diag_plot <- tryCatch(
        make_diagnostic_plots(data, result, pheno, row_col, col_col, ""),
        error = function(e) { warning("Diagnostic plot failed: ", e$message); NULL }
      )
      if (!is.null(diag_plot)) {
        out_png <- file.path(output_dir,
                             paste0("diagnostics_", pheno, "_single.png"))
        ggsave(out_png, diag_plot, width = 12, height = 8, dpi = 150)
        cat("    Saved:", out_png, "\n")
      }

      re_plot <- tryCatch(
        plot_row_col_re(data, result, row_col, col_col, ""),
        error = function(e) { warning("RE plot failed: ", e$message); NULL }
      )
      if (!is.null(re_plot)) {
        out_png <- file.path(output_dir,
                             paste0("row_col_re_", pheno, "_single.png"))
        ggsave(out_png, re_plot, width = 8, height = 4, dpi = 150)
        cat("    Saved:", out_png, "\n")
      }

    } else {
      # Multi-bench path
      result <- fit_mgcv_joint(
        data          = data,
        pheno         = pheno,
        geno_col      = geno_col,
        row_col       = row_col,
        col_col       = col_col,
        bench_col     = bench_col,
        k_row         = k_row,
        k_col         = k_col,
        estimate_type = estimate_type
      )

      if (!is.null(result$blues)) {
        df <- result$blues
        names(df)[names(df) == "BLUE"]    <- pheno
        names(df)[names(df) == "BLUE_SE"] <- paste0(pheno, "_SE")
        all_blues_list[[pheno]] <- df
      }

      if (!is.null(result$spatial_smooth)) {
        spatial_rows[[length(spatial_rows) + 1]] <- data.frame(
          data[, c(geno_col, row_col, col_col, bench_col), drop = FALSE],
          pheno          = pheno,
          spatial_smooth = result$spatial_smooth,
          spatial_total  = result$spatial_total,
          stringsAsFactors = FALSE
        )
      }

      resid_sd_val <- if (!is.null(result$residuals))
        sd(result$residuals, na.rm = TRUE) else NA

      for (b in bench_ids) {
        b_mask  <- data[[bench_col]] == b
        b_resid <- if (!is.null(result$residuals)) result$residuals[b_mask] else NULL
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          pheno       = pheno,
          bench       = as.character(b),
          resid_sd    = if (!is.null(b_resid)) sd(b_resid, na.rm = TRUE) else NA,
          edf_spatial = result$edf_spatial,
          converged   = result$converged,
          stringsAsFactors = FALSE
        )
      }

      for (b in bench_ids) {
        b_mask  <- data[[bench_col]] == b
        b_data  <- data[b_mask, , drop = FALSE]
        b_result <- list(
          spatial_smooth = if (!is.null(result$spatial_smooth)) result$spatial_smooth[b_mask] else NULL,
          residuals      = if (!is.null(result$residuals))      result$residuals[b_mask]      else NULL,
          fitted         = if (!is.null(result$fitted))         result$fitted[b_mask]          else NULL,
          row_re         = if (!is.null(result$row_re))         result$row_re[b_mask]          else NULL,
          col_re         = if (!is.null(result$col_re))         result$col_re[b_mask]          else NULL
        )
        diag_plot <- tryCatch(
          make_diagnostic_plots(b_data, b_result, pheno, row_col, col_col,
                                paste0("Bench: ", b)),
          error = function(e) { warning("Diagnostic plot failed (bench ", b, "): ",
                                        e$message); NULL }
        )
        if (!is.null(diag_plot)) {
          safe_bench <- gsub("[^A-Za-z0-9_.-]", "_", as.character(b))
          out_png <- file.path(output_dir,
                               paste0("diagnostics_", pheno, "_", safe_bench, ".png"))
          ggsave(out_png, diag_plot, width = 12, height = 8, dpi = 150)
          cat("    Saved:", out_png, "\n")
        }

        re_plot <- tryCatch(
          plot_row_col_re(b_data, b_result, row_col, col_col,
                          paste0("Bench: ", b)),
          error = function(e) { warning("RE plot failed (bench ", b, "): ",
                                        e$message); NULL }
        )
        if (!is.null(re_plot)) {
          out_png <- file.path(output_dir,
                               paste0("row_col_re_", pheno, "_", safe_bench, ".png"))
          ggsave(out_png, re_plot, width = 8, height = 4, dpi = 150)
          cat("    Saved:", out_png, "\n")
        }
      }
    }
  }

  # --------------------------------------------------------------------------
  # CSV outputs
  # --------------------------------------------------------------------------
  cat("\n-- Writing CSV outputs ---------------------------------------------------\n")

  blues_df <- NULL
  blups_df <- NULL
  sp_df    <- NULL
  sum_df   <- NULL

  if (length(all_blues_list) > 0) {
    blues_df <- Reduce(function(a, b) merge(a, b, by = geno_col, all = TRUE),
                       all_blues_list)
    out_path <- file.path(output_dir, "BLUEs.csv")
    write.csv(blues_df, out_path, row.names = FALSE)
    cat("  BLUEs:", out_path, "(", nrow(blues_df), "genotypes )\n")
  }

  if (length(all_blups_list) > 0) {
    blups_df <- Reduce(function(a, b) merge(a, b, by = geno_col, all = TRUE),
                       all_blups_list)
    out_path <- file.path(output_dir, "BLUPs.csv")
    write.csv(blups_df, out_path, row.names = FALSE)
    cat("  BLUPs:", out_path, "(", nrow(blups_df), "genotypes )\n")
  }

  if (length(spatial_rows) > 0) {
    sp_df    <- do.call(rbind, spatial_rows)
    out_path <- file.path(output_dir, "spatial_trends.csv")
    write.csv(sp_df, out_path, row.names = FALSE)
    cat("  Spatial trends:", out_path, "(", nrow(sp_df), "rows )\n")
  }

  if (length(summary_rows) > 0) {
    sum_df   <- do.call(rbind, summary_rows)
    out_path <- file.path(output_dir, "model_summary.csv")
    write.csv(sum_df, out_path, row.names = FALSE)
    cat("  Model summary:", out_path, "\n")
    print(sum_df)
  }

  # --------------------------------------------------------------------------
  # Summary plots
  # --------------------------------------------------------------------------
  cat("\n-- Saving summary plots --------------------------------------------------\n")

  if (length(all_blues_list) > 0) {
    blues_long <- do.call(rbind, lapply(names(all_blues_list), function(ph) {
      df <- all_blues_list[[ph]]
      data.frame(trait = ph, value = df[[ph]], stringsAsFactors = FALSE)
    }))
    blues_long <- blues_long[!is.na(blues_long$value), ]

    p_dens <- ggplot(blues_long, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30, fill = "#2c7bb6", alpha = 0.6) +
      geom_density(colour = "red", linewidth = 0.7) +
      facet_wrap(~ trait, scales = "free") +
      labs(title = "BLUE distribution per trait",
           x = "BLUE value", y = "Density") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    out_png <- file.path(output_dir, "blues_distribution.png")
    ggsave(out_png, p_dens, width = max(6, 4 * length(all_blues_list)), height = 5,
           dpi = 150)
    cat("  BLUE distribution:", out_png, "\n")
  }

  if (!is.null(sp_df) && !is.null(sum_df)) {
    bench_col_name <- if (is_multibench) bench_col else "bench"

    sp_plots <- lapply(seq_len(nrow(sum_df)), function(i) {
      ph <- sum_df$pheno[i]
      b  <- sum_df$bench[i]
      if (is_multibench) {
        sub <- sp_df[sp_df$pheno == ph & sp_df[[bench_col]] == b, ]
      } else {
        sub <- sp_df[sp_df$pheno == ph, ]
      }
      if (nrow(sub) == 0) return(NULL)
      sub_plot <- sub[, c(row_col, col_col, "spatial_smooth"), drop = FALSE]
      plot_heatmap(sub_plot, row_col, col_col, "spatial_smooth",
                   title = paste0(ph, " | ", b))
    })
    sp_plots <- Filter(Negate(is.null), sp_plots)

    if (length(sp_plots) > 0) {
      n_col_grid <- min(length(sp_plots), 3L)
      composite  <- wrap_plots(sp_plots, ncol = n_col_grid) +
        plot_annotation(
          title = "Spatial smooth surfaces (te only)",
          theme = theme(plot.title = element_text(size = 13, face = "bold"))
        )
      out_png <- file.path(output_dir, "spatial_surfaces.png")
      n_rows_grid <- ceiling(length(sp_plots) / n_col_grid)
      ggsave(out_png, composite,
             width  = 5 * n_col_grid,
             height = 4 * n_rows_grid,
             dpi = 150)
      cat("  Spatial surfaces:", out_png, "\n")
    }
  }

  cat("\n-- Done ------------------------------------------------------------------\n")
  cat("  Results written to:", output_dir, "\n")

  invisible(list(
    blues          = blues_df,
    blups          = blups_df,
    spatial_trends = sp_df,
    model_summary  = sum_df
  ))
}
