# ============================================================================
# fit_spatial_models.R
# Multi-method spatial correction framework for field trial data
#
# Methods:
#   SpATS  — 2D P-spline via SAP(); adaptive nseg; BLUEs or BLUPs
#   mgcv   — gam() with te(row, col); automatic smoothness via REML;
#             BLUEs (0 + Genotype fixed) or BLUPs (s(Genotype, bs="re"))
#             Smoother variants: mgcv (default cr), mgcv_ps (P-spline),
#             mgcv_ps_re (P-spline + row/col RE), mgcv_highk, mgcv_ad (adaptive)
#   sommer — mmes() with spl2Dc(); spatial as random effect;
#             BLUEs (Genotype fixed) or BLUPs (Genotype random)
#
# Outputs per run:
#   spatial_diagnostics.pdf   — 6-panel diagnostic page per trait × method;
#                               comparison page when method="all"
#   BLUEs.csv / BLUPs.csv     — one row per genotype (× bench); all original
#                               columns preserved; SE columns appended
#   fitted_values.csv         — one row per plot; all original columns + fitted,
#                               residual, and spatial-trend columns
#   outlier_report.txt        — per-trait outlier count and values
#
# Usage:
#   source("scripts/fit_spatial_models.R")
#   run_spatial_correction(
#     fn         = "data/wheatdata.rda",
#     rda_object = "wheatdata",
#     trait_cols = "yield",
#     gt_col     = "geno",
#     row_col    = "row",
#     col_col    = "col",
#     method     = "all",
#     output_type = "BLUEs",
#     output_dir = "output"
#   )
# ============================================================================

suppressPackageStartupMessages({
  library(SpATS)
  library(mgcv)
  library(sommer)
  library(ggplot2)
  library(patchwork)
})

# Load shared utility and mgcv fitting functions from the canonical
# standalone script.
source("scripts/spatial_correct_gam.R")


# ============================================================================
# Section 1: Utility functions
# ============================================================================

#' Replace values outside 1.5 × IQR fences with NA
#'
#' @param x  Numeric vector
#' @return   Named list: $x (cleaned vector), $n_outliers, $outlier_values
replace_outliers_with_na <- function(x) {
  Q1  <- quantile(x, 0.25, na.rm = TRUE)
  Q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- Q3 - Q1
  lo  <- Q1 - 1.5 * iqr
  hi  <- Q3 + 1.5 * iqr
  is_out <- !is.na(x) & (x < lo | x > hi)
  list(
    x             = replace(x, is_out, NA),
    n_outliers    = sum(is_out),
    outlier_values = x[is_out]
  )
}



# ============================================================================
# Section 2: Per-bench correction functions
# ============================================================================
# Each function returns a named list:
#   $blues      data.frame(Genotype, BLUE, BLUE_SE)     or NULL
#   $blups      data.frame(Genotype, BLUP, BLUP_SE)     or NULL
#   $fitted     numeric vector (length = nrow(bench_data))
#   $residuals  numeric vector
#   $spatial    numeric vector (spatial trend at each plot)
# ============================================================================


#' SpATS spatial correction for one bench
#'
#' Uses SAP() smoother with adaptive nseg. Row/col random effects included.
#' BLUEs: genotype.as.random = FALSE; BLUPs: genotype.as.random = TRUE.
#'
#' @param bench_data  data.frame for a single bench
#' @param pheno       Name of the phenotype column
#' @param gt_col      Name of the genotype column
#' @param row_col     Name of the row coordinate column
#' @param col_col     Name of the column coordinate column
#' @param output_type "BLUEs", "BLUPs", or "both"
#' @return standardised result list
run_SpATS_bench <- function(bench_data, pheno, gt_col, row_col, col_col,
                             output_type = "BLUEs") {
  res <- list(blues = NULL, blups = NULL, fitted = NULL,
              residuals = NULL, spatial = NULL)

  # SpATS needs factor row/col variables named rnf / cnf in the data
  bench_data[[gt_col]] <- as.factor(bench_data[[gt_col]])
  bench_data$rnf       <- as.factor(bench_data[[row_col]])
  bench_data$cnf       <- as.factor(bench_data[[col_col]])

  n_row    <- length(unique(bench_data[[row_col]]))
  n_col    <- length(unique(bench_data[[col_col]]))
  nseg_row <- adaptive_nseg(n_row)
  nseg_col <- adaptive_nseg(n_col)

  spatial_fmla <- as.formula(
    paste0("~ SAP(", col_col, ", ", row_col,
           ", nseg = c(", nseg_col, ", ", nseg_row, "))")
  )

  # Helper: run SpATS and extract genotype predictions
  extract_preds <- function(model, type) {
    preds <- predict(model, which = gt_col)
    col_pred <- grep("predicted", names(preds), value = TRUE)
    col_se   <- grep("^se$|errors|^se_", names(preds), value = TRUE)
    if (length(col_pred) == 0) col_pred <- names(preds)[2]
    if (length(col_se)   == 0) col_se   <- names(preds)[3]
    data.frame(
      Genotype = as.character(preds[[gt_col]]),
      Value    = preds[[col_pred]],
      SE       = preds[[col_se]],
      stringsAsFactors = FALSE
    )
  }

  primary_model <- NULL

  if (output_type %in% c("BLUEs", "both")) {
    m <- tryCatch(
      SpATS(response = pheno, spatial = spatial_fmla, genotype = gt_col,
            random = ~ rnf + cnf, data = bench_data,
            control = list(tolerance = 1e-04),
            genotype.as.random = FALSE),
      error = function(e) { message("SpATS BLUEs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      primary_model <- m
      p <- extract_preds(m, "BLUE")
      res$blues <- data.frame(Genotype = p$Genotype, BLUE = p$Value,
                               BLUE_SE = p$SE, stringsAsFactors = FALSE)
    }
  }

  if (output_type %in% c("BLUPs", "both")) {
    m <- tryCatch(
      SpATS(response = pheno, spatial = spatial_fmla, genotype = gt_col,
            random = ~ rnf + cnf, data = bench_data,
            control = list(tolerance = 1e-04),
            genotype.as.random = TRUE),
      error = function(e) { message("SpATS BLUPs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      if (is.null(primary_model)) primary_model <- m
      p <- extract_preds(m, "BLUP")
      res$blups <- data.frame(Genotype = p$Genotype, BLUP = p$Value,
                               BLUP_SE = p$SE, stringsAsFactors = FALSE)
    }
  }

  if (!is.null(primary_model)) {
    res$fitted    <- as.numeric(primary_model$fitted)
    res$residuals <- as.numeric(primary_model$residuals)

    # Spatial trend = fitted minus per-plot genotype contribution
    geno_preds <- predict(primary_model, which = gt_col)
    col_pred   <- grep("predicted", names(geno_preds), value = TRUE)
    if (length(col_pred) == 0) col_pred <- names(geno_preds)[2]
    geno_effect <- geno_preds[[col_pred]][
      match(as.character(bench_data[[gt_col]]),
            as.character(geno_preds[[gt_col]]))
    ]
    res$spatial <- res$fitted - geno_effect

    # SpATS's PSANOVA surface is not constrained to be zero-mean at the combined
    # level, so mean(spatial) can be substantially non-zero. This causes
    # predict(model, which = gt_col) to return BLUEs/BLUPs that are
    # ~mean(spatial) below the phenotype scale. Centre the spatial surface and
    # add its mean back to BLUEs/BLUPs so all methods are on the same absolute
    # scale.
    sp_mean <- mean(res$spatial, na.rm = TRUE)
    res$spatial <- res$spatial - sp_mean
    if (!is.null(res$blues)) res$blues$BLUE <- res$blues$BLUE + sp_mean
    if (!is.null(res$blups)) res$blups$BLUP <- res$blups$BLUP + sp_mean
  }

  res
}


#' mgcv::gam spatial correction for one bench
#'
#' Uses te(row, col) tensor-product smooth (anisotropic, smoothness via REML).
#' BLUEs: 0 + Genotype as parametric fixed; BLUPs: s(Genotype, bs="re").
#'
#' The best-performing variant (te_ps_re) delegates to fit_mgcv_bench() from
#' spatial_correct_gam.R. Other smoother variants are handled inline.
#'
#' @inheritParams run_SpATS_bench
#' @param smoother_type One of "te_default", "te_ps", "te_ps_re", "te_highk",
#'   "te_ad". Controls the spatial smoother basis and structure.
#' @return standardised result list
run_mgcv_bench <- function(bench_data, pheno, gt_col, row_col, col_col,
                            output_type = "BLUEs",
                            smoother_type = "te_default") {

  # te_ps_re: delegate to the canonical implementation in spatial_correct_gam.R
  if (smoother_type == "te_ps_re") {
    r <- fit_mgcv_bench(bench_data, pheno, gt_col, row_col, col_col,
                        estimate_type = output_type)
    # Normalize first column name to "Genotype" (internal convention)
    if (!is.null(r$blues)) names(r$blues)[1] <- "Genotype"
    if (!is.null(r$blups)) names(r$blups)[1] <- "Genotype"
    return(r)
  }

  res <- list(blues = NULL, blups = NULL, fitted = NULL,
              residuals = NULL, spatial = NULL)

  bench_data[[gt_col]] <- as.factor(bench_data[[gt_col]])
  geno_levels          <- levels(bench_data[[gt_col]])
  n_geno               <- length(geno_levels)

  n_row    <- length(unique(bench_data[[row_col]]))
  n_col    <- length(unique(bench_data[[col_col]]))
  nr       <- adaptive_nseg(n_row)
  nc       <- adaptive_nseg(n_col)

  te_term <- switch(smoother_type,
    te_default = paste0("te(", row_col, ", ", col_col, ")"),
    te_ps      = paste0("te(", row_col, ", ", col_col,
                        ", bs=c('ps','ps'), k=c(", nr, ",", nc, "))"),
    te_highk   = paste0("te(", row_col, ", ", col_col,
                        ", k=c(", nr, ",", nc, "))"),
    te_ad      = paste0("t2(", row_col, ", ", col_col,
                        ", bs=c('tp','tp'), k=c(", nr, ",", nc, "), full=TRUE)"),
    stop("Unknown smoother_type: ", smoother_type)
  )

  # ---- BLUEs ----------------------------------------------------------
  if (output_type %in% c("BLUEs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ 0 + ", gt_col, " + ", te_term))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { message("mgcv BLUEs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      coef_names <- paste0(gt_col, geno_levels)
      blues_vals <- coef(m)[coef_names]
      vcov_diag  <- diag(vcov(m))
      blues_se   <- sqrt(vcov_diag[coef_names])

      res$blues <- data.frame(
        Genotype = geno_levels,
        BLUE     = as.numeric(blues_vals),
        BLUE_SE  = as.numeric(blues_se),
        stringsAsFactors = FALSE
      )

      if (is.null(res$fitted)) {
        res$fitted    <- as.numeric(fitted(m))
        res$residuals <- as.numeric(residuals(m))
        terms_pred  <- predict(m, type = "terms")
        te_col      <- grep("^t[e2]\\(", colnames(terms_pred), value = TRUE)
        res$spatial <- as.numeric(terms_pred[, te_col[1]])
      }
    }
  }

  # ---- BLUPs ----------------------------------------------------------
  if (output_type %in% c("BLUPs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ s(", gt_col, ", bs='re') + ", te_term))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { message("mgcv BLUPs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      re_pattern <- paste0("s\\(", gt_col, "\\)")
      re_idx     <- grep(re_pattern, names(coef(m)))
      blup_vals  <- coef(m)[re_idx]

      mean_row <- mean(bench_data[[row_col]], na.rm = TRUE)
      mean_col <- mean(bench_data[[col_col]], na.rm = TRUE)
      newdat   <- setNames(
        data.frame(geno_levels, mean_row, mean_col),
        c(gt_col, row_col, col_col)
      )
      newdat[[gt_col]] <- factor(newdat[[gt_col]], levels = geno_levels)
      tp <- tryCatch(
        predict(m, newdata = newdat, type = "terms", se.fit = TRUE),
        error = function(e) NULL
      )
      re_col_name <- if (!is.null(tp))
        grep(re_pattern, colnames(tp$fit), value = TRUE)[1] else NA

      blup_se <- if (!is.null(tp) && !is.na(re_col_name))
        as.numeric(tp$se.fit[, re_col_name]) else rep(NA_real_, n_geno)

      intercept  <- coef(m)["(Intercept)"]
      blup_total <- as.numeric(intercept) + as.numeric(blup_vals)

      res$blups <- data.frame(
        Genotype = geno_levels,
        BLUP     = blup_total,
        BLUP_SE  = blup_se,
        stringsAsFactors = FALSE
      )

      if (is.null(res$fitted)) {
        res$fitted    <- as.numeric(fitted(m))
        res$residuals <- as.numeric(residuals(m))
        terms_pred    <- predict(m, type = "terms")
        te_col        <- grep("^te\\(", colnames(terms_pred), value = TRUE)
        res$spatial   <- as.numeric(terms_pred[, te_col[1]])
      }
    }
  }

  res
}


#' sommer::mmes spatial correction for one bench
#'
#' Uses spl2Dc(col, row) 2-D P-spline as a random effect (PSANOVA type).
#' BLUEs: Genotype fixed; BLUPs: Genotype random.
#' Note: spl2Dc() requires sommer::mmes() (not mmer()).
#'
#' @inheritParams run_SpATS_bench
#' @return standardised result list
run_sommer_bench <- function(bench_data, pheno, gt_col, row_col, col_col,
                              output_type = "BLUEs") {
  res <- list(blues = NULL, blups = NULL, fitted = NULL,
              residuals = NULL, spatial = NULL)

  bench_data[[gt_col]] <- as.character(bench_data[[gt_col]])
  geno_sorted          <- sort(unique(bench_data[[gt_col]]))
  n_geno               <- length(geno_sorted)

  n_row    <- length(unique(bench_data[[row_col]]))
  n_col    <- length(unique(bench_data[[col_col]]))
  nseg_row <- adaptive_nseg(n_row)
  nseg_col <- adaptive_nseg(n_col)

  # Row/col factor variables for sommer random terms
  bench_data$rnf <- as.factor(bench_data[[row_col]])
  bench_data$cnf <- as.factor(bench_data[[col_col]])

  spl_term  <- paste0("spl2Dc(", col_col, ", ", row_col,
                      ", nsegments = c(", nseg_col, ", ", nseg_row, "))")
  # rand_base starts WITHOUT "~" so it can be combined with optional gt_col term
  rand_base <- "rnf + cnf + "

  # ---- BLUEs ----------------------------------------------------------
  if (output_type %in% c("BLUEs", "both")) {
    fm_fixed  <- as.formula(paste0(pheno, " ~ ", gt_col))
    fm_random <- as.formula(paste0("~ ", rand_base, spl_term))
    m <- tryCatch(
      suppressMessages(
        mmes(fixed = fm_fixed, random = fm_random,
             rcov = ~ units, data = bench_data,
             verbose = FALSE, dateWarning = FALSE)
      ),
      error = function(e) { message("sommer BLUEs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      # mmes uses treatment contrasts: b[1] = intercept (= BLUE for reference
      # genotype), b[2:n] = contrasts. Recover absolute BLUEs via X * b.
      uniq_df <- data.frame(
        setNames(list(factor(geno_sorted, levels = geno_sorted)), gt_col)
      )
      X_uniq    <- model.matrix(as.formula(paste0("~ ", gt_col)), data = uniq_df)
      blue_vals <- as.numeric(X_uniq %*% m$b)

      n_b <- length(m$b)
      blue_se <- if (!is.null(m$Ci) && nrow(m$Ci) >= n_b) {
        fixed_cov <- m$Ci[seq_len(n_b), seq_len(n_b), drop = FALSE]
        sqrt(diag(X_uniq %*% fixed_cov %*% t(X_uniq)))
      } else {
        rep(NA_real_, n_geno)
      }

      res$blues <- data.frame(
        Genotype = geno_sorted,
        BLUE     = blue_vals,
        BLUE_SE  = blue_se,
        stringsAsFactors = FALSE
      )

      if (is.null(res$fitted)) {
        # fitted(m) = Xb + Zu (full predictions); m$fitted = Xb only
        res$fitted    <- as.numeric(fitted(m))
        res$residuals <- bench_data[[pheno]] - res$fitted
        # Spatial trend = full fitted minus per-plot genotype (fixed) contribution
        geno_match   <- match(bench_data[[gt_col]], geno_sorted)
        geno_contrib <- blue_vals[geno_match]
        res$spatial  <- res$fitted - geno_contrib
      }
    }
  }

  # ---- BLUPs ----------------------------------------------------------
  if (output_type %in% c("BLUPs", "both")) {
    fm_fixed  <- as.formula(paste0(pheno, " ~ 1"))
    fm_random <- as.formula(paste0("~ ", gt_col, " + ", rand_base, spl_term))
    m <- tryCatch(
      suppressMessages(
        mmes(fixed = fm_fixed, random = fm_random,
             rcov = ~ units, data = bench_data,
             verbose = FALSE, dateWarning = FALSE)
      ),
      error = function(e) { message("sommer BLUPs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      intercept <- as.numeric(m$b[1])
      # Identify genotype random effect component
      geno_key <- grep(paste0("ism.*", gt_col, "|", gt_col, ".*ism"),
                       names(m$uList), value = TRUE)[1]
      if (is.na(geno_key)) geno_key <- names(m$uList)[1]

      blup_mat   <- m$uList[[geno_key]]
      blup_names <- rownames(blup_mat)
      blup_vals  <- as.numeric(blup_mat[, 1]) + intercept
      names(blup_vals) <- blup_names

      # PEV-based SE
      pev <- m$uPevList[[geno_key]]
      if (!is.null(pev)) {
        blup_se_raw <- sqrt(diag(pev))
        # pev may be scalar (balanced design) or full matrix
        if (length(blup_se_raw) == 1)
          blup_se <- rep(blup_se_raw, length(blup_vals))
        else
          blup_se <- blup_se_raw[blup_names]
      } else {
        blup_se <- rep(NA_real_, length(blup_vals))
      }

      res$blups <- data.frame(
        Genotype = blup_names,
        BLUP     = blup_vals,
        BLUP_SE  = blup_se,
        stringsAsFactors = FALSE
      )

      if (is.null(res$fitted)) {
        # fitted(m) = Xb + Zu (full predictions); m$fitted = Xb only
        res$fitted    <- as.numeric(fitted(m))
        res$residuals <- bench_data[[pheno]] - res$fitted
        # Spatial trend = full fitted minus per-plot genotype (random) contribution
        geno_match   <- match(bench_data[[gt_col]], blup_names)
        geno_contrib <- blup_vals[geno_match]
        res$spatial  <- res$fitted - geno_contrib
      }
    }
  }

  res
}


# ============================================================================
# Section 2b: Joint multi-bench correction functions
# ============================================================================
# Each function takes the full multi-bench dataset and returns:
#   $blues      data.frame(Genotype, BLUE, BLUE_SE)  — one row per genotype
#   $fitted     numeric vector (length = nrow(data))
#   $residuals  numeric vector
#   $spatial    numeric vector (spatial trend at each plot)
# ============================================================================


#' Joint multi-bench mgcv::gam spatial correction
#'
#' Delegates to fit_mgcv_joint() from spatial_correct_gam.R (the canonical
#' implementation). Returns list(blues, fitted, residuals, spatial) with the
#' genotype column renamed to "Genotype" for internal consistency.
#'
#' @param data       data.frame with all benches
#' @param pheno      Name of the phenotype column
#' @param gt_col     Name of the genotype column
#' @param row_col    Name of the row coordinate column
#' @param col_col    Name of the column coordinate column
#' @param bench_col  Name of the bench column
#' @return list(blues, fitted, residuals, spatial)
run_mgcv_joint <- function(data, pheno, gt_col, row_col, col_col, bench_col) {
  r <- fit_mgcv_joint(data, pheno, gt_col, row_col, col_col, bench_col)
  if (!is.null(r$blues)) names(r$blues)[1] <- "Genotype"
  r
}


#' Joint multi-bench sommer::mmes spatial correction
#'
#' Fits a single mixed model across all benches with:
#'   - Genotype + Bench as fixed effects
#'   - Effectively per-bench spl2Dc() via column-offset trick (each bench's
#'     columns are shifted to a non-overlapping range so the single spl2Dc
#'     surface is effectively independent per bench)
#'   - Row and column random effects nested within bench
#'
#' BLUEs recovered via X_uniq %*% m$b, averaging over bench levels.
#'
#' @inheritParams run_mgcv_joint
#' @return list(blues, fitted, residuals, spatial)
run_sommer_joint <- function(data, pheno, gt_col, row_col, col_col, bench_col) {
  res <- list(blues = NULL, fitted = NULL, residuals = NULL, spatial = NULL)

  data[[gt_col]]  <- as.factor(data[[gt_col]])
  data$Bench_f    <- as.factor(data[[bench_col]])
  data$rnf        <- as.factor(data[[row_col]])
  data$cnf        <- as.factor(data[[col_col]])
  geno_sorted     <- sort(levels(data[[gt_col]]))
  bench_levels    <- levels(data$Bench_f)

  # Adaptive nseg per bench
  nr_vec <- sapply(bench_levels, function(b) {
    length(unique(data[[row_col]][data$Bench_f == b]))
  })
  nc_vec <- sapply(bench_levels, function(b) {
    length(unique(data[[col_col]][data$Bench_f == b]))
  })
  nseg_row <- adaptive_nseg(min(nr_vec))
  nseg_col_per <- adaptive_nseg(min(nc_vec))

  # Column-offset trick: shift each bench's column coordinates to non-overlapping

  # ranges separated by a large gap, so a single spl2Dc effectively fits
  # independent surfaces per bench
  col_max <- max(data[[col_col]])
  gap     <- col_max * 10   # 10× gap prevents cross-bench smoothing
  data$col_offset <- data[[col_col]]
  for (i in seq_along(bench_levels)) {
    mask <- data$Bench_f == bench_levels[i]
    data$col_offset[mask] <- data[[col_col]][mask] + (i - 1) * gap
  }

  # Total nseg: per-bench segments × number of benches
  nseg_col_total <- nseg_col_per * length(bench_levels)

  # Fixed formula
  fm_fixed <- as.formula(paste0(pheno, " ~ ", gt_col, " + Bench_f"))

  # Random formula with offset columns and bench-nested row/col RE
  fm_random <- as.formula(paste0(
    "~ spl2Dc(col_offset, ", row_col,
    ", nsegments = c(", nseg_col_total, ", ", nseg_row, "))",
    " + vsm(dsm(Bench_f), ism(rnf))",
    " + vsm(dsm(Bench_f), ism(cnf))"
  ))

  m <- tryCatch(
    suppressMessages(
      mmes(fixed = fm_fixed, random = fm_random,
           rcov = ~ units, data = data,
           verbose = FALSE, dateWarning = FALSE)
    ),
    error = function(e) {
      message("sommer joint error: ", e$message)
      NULL
    }
  )

  if (is.null(m)) return(res)

  # Extract BLUEs: recover absolute values via X_uniq %*% m$b
  # Average over bench levels: predict at each bench, then average per genotype
  n_bench <- length(bench_levels)
  n_geno  <- length(geno_sorted)

  # Build prediction design matrix: one row per geno, averaged over benches
  pred_list <- lapply(bench_levels, function(b) {
    data.frame(
      setNames(list(factor(geno_sorted, levels = levels(data[[gt_col]])),
                    factor(b, levels = bench_levels)),
               c(gt_col, "Bench_f"))
    )
  })
  pred_df <- do.call(rbind, pred_list)
  X_pred  <- model.matrix(as.formula(paste0("~ ", gt_col, " + Bench_f")), data = pred_df)

  # Predicted fixed values for all geno × bench combinations
  pred_vals <- as.numeric(X_pred %*% m$b)
  # Reshape to n_geno × n_bench and average across benches
  pred_mat  <- matrix(pred_vals, nrow = n_geno, ncol = n_bench)
  blue_vals <- rowMeans(pred_mat)

  # SE via variance of averaged predictions
  n_b <- length(m$b)
  if (!is.null(m$Ci) && nrow(m$Ci) >= n_b) {
    fixed_cov <- m$Ci[seq_len(n_b), seq_len(n_b), drop = FALSE]
    # Averaging matrix: for each genotype, average the X rows across benches
    X_avg <- matrix(0, nrow = n_geno, ncol = ncol(X_pred))
    for (i in seq_len(n_bench)) {
      X_avg <- X_avg + X_pred[((i - 1) * n_geno + 1):(i * n_geno), , drop = FALSE]
    }
    X_avg <- X_avg / n_bench
    blue_se <- sqrt(diag(X_avg %*% fixed_cov %*% t(X_avg)))
  } else {
    blue_se <- rep(NA_real_, n_geno)
  }

  res$blues <- data.frame(
    Genotype = geno_sorted,
    BLUE     = blue_vals,
    BLUE_SE  = blue_se,
    stringsAsFactors = FALSE
  )

  # Fitted values and residuals
  res$fitted    <- as.numeric(fitted(m))
  res$residuals <- data[[pheno]] - res$fitted

  # Spatial trend = fitted minus per-plot genotype + bench fixed contribution
  X_data     <- model.matrix(as.formula(paste0("~ ", gt_col, " + Bench_f")), data = data)
  fixed_pred <- as.numeric(X_data %*% m$b)
  res$spatial <- res$fitted - fixed_pred

  res
}


# ============================================================================
# Section 3: Visualisation functions
# ============================================================================
# plot_heatmap(), plot_variogram(), and empirical_semivariogram() are defined
# in spatial_correct_gam.R (sourced above).


#' Compose 6-panel diagnostic page for one bench × trait × method
#'
#' Panels: raw heatmap | fitted heatmap | spatial-trend heatmap |
#'         residuals heatmap | raw variogram | residual variogram
#'
#' @param bench_data  data.frame for the bench including the observed pheno
#' @param result      Return value from run_*_bench()
#' @param pheno       Name of the phenotype column in bench_data
#' @param method      Character label for the page title
#' @param row_col     Name of row coordinate column
#' @param col_col     Name of column coordinate column
#' @param bench_label Character label for bench; used in title
#' @return patchwork ggplot
plot_diagnostics <- function(bench_data, result, pheno, method,
                              row_col, col_col, bench_label = "") {
  obs <- bench_data[[pheno]]

  # Build a working data.frame with coordinates and derived quantities
  pd <- bench_data[, c(row_col, col_col), drop = FALSE]
  pd$obs     <- obs
  pd$fitted  <- if (!is.null(result$fitted))  result$fitted  else NA_real_
  pd$spatial <- if (!is.null(result$spatial)) result$spatial else NA_real_
  pd$resid   <- if (!is.null(result$residuals)) result$residuals else obs - pd$fitted

  # Shared colour limits for obs and fitted (comparable scale)
  # Guard against all-NA ranges (model failed to produce values)
  safe_range <- function(...) {
    r <- range(..., na.rm = TRUE)
    if (any(!is.finite(r))) c(0, 1) else r
  }
  lim_obs <- safe_range(pd$obs, pd$fitted)
  lim_sp  <- safe_range(pd$spatial)
  lim_res <- safe_range(pd$resid)

  p1 <- plot_heatmap(pd, row_col, col_col, "obs",
                     title = "Observed", limits = lim_obs)
  p2 <- plot_heatmap(pd, row_col, col_col, "fitted",
                     title = "Fitted", limits = lim_obs)
  p3 <- plot_heatmap(pd, row_col, col_col, "spatial",
                     title = "Spatial trend", limits = lim_sp)
  p4 <- plot_heatmap(pd, row_col, col_col, "resid",
                     title = "Residuals", limits = lim_res)

  vg_obs  <- empirical_semivariogram(pd, row_col, col_col, "obs")
  vg_res  <- empirical_semivariogram(pd, row_col, col_col, "resid")
  p5 <- plot_variogram(vg_obs,  title = "Variogram: observed")
  p6 <- plot_variogram(vg_res,  title = "Variogram: residuals")

  page_title <- paste0(pheno, "  |  Method: ", method,
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_annotation(
      title = page_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
}


#' Comparison page for method = "all": spatial trends and BLUE scatter plots
#'
#' @param bench_data    data.frame for the bench
#' @param results_list  Named list of results; names = method labels
#' @param pheno         Phenotype name
#' @param row_col       Row coordinate column
#' @param col_col       Column coordinate column
#' @param bench_label   Character label for the bench
#' @return patchwork ggplot
plot_comparison <- function(bench_data, results_list, pheno, row_col, col_col,
                             bench_label = "") {
  methods <- names(results_list)
  n_meth  <- length(methods)

  # Spatial-trend heatmaps (top row)
  sp_plots <- lapply(methods, function(mn) {
    r   <- results_list[[mn]]
    pd  <- bench_data[, c(row_col, col_col), drop = FALSE]
    pd$spatial <- if (!is.null(r$spatial)) r$spatial else NA_real_
    lim <- range(pd$spatial, na.rm = TRUE)
    plot_heatmap(pd, row_col, col_col, "spatial",
                 title = paste0("Spatial trend: ", mn), limits = lim)
  })

  # Pairwise BLUE scatter plots (lower section)
  # Collect all BLUEs into one data.frame
  blue_list <- lapply(methods, function(mn) {
    r <- results_list[[mn]]
    if (!is.null(r$blues))
      data.frame(Genotype = r$blues$Genotype, v = r$blues$BLUE, method = mn,
                 stringsAsFactors = FALSE)
    else if (!is.null(r$blups))
      data.frame(Genotype = r$blups$Genotype, v = r$blups$BLUP, method = mn,
                 stringsAsFactors = FALSE)
    else NULL
  })
  blue_list <- Filter(Negate(is.null), blue_list)

  scatter_plots <- list()
  if (length(blue_list) >= 2) {
    # All unique method pairs
    pairs <- combn(seq_along(blue_list), 2, simplify = FALSE)
    scatter_plots <- lapply(pairs, function(pr) {
      da <- blue_list[[pr[1]]]
      db <- blue_list[[pr[2]]]
      merged <- merge(da, db, by = "Genotype", suffixes = c("_x", "_y"))
      r_val  <- if (nrow(merged) > 2)
        round(cor(merged$v_x, merged$v_y, use = "complete.obs"), 3) else NA
      ggplot(merged, aes(x = v_x, y = v_y)) +
        geom_point(size = 0.8, alpha = 0.6, colour = "#2c7bb6") +
        geom_abline(slope = 1, intercept = 0, colour = "red",
                    linetype = "dashed", linewidth = 0.4) +
        labs(
          title = paste0(da$method[1], " vs ", db$method[1]),
          subtitle = paste0("r = ", r_val),
          x = paste0("BLUE (", da$method[1], ")"),
          y = paste0("BLUE (", db$method[1], ")")
        ) +
        theme_bw(base_size = 8) +
        theme(
          plot.title    = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 7, hjust = 0.5),
          axis.title    = element_text(size = 7),
          axis.text     = element_text(size = 6)
        )
    })
  }

  page_title <- paste0(pheno, "  |  Method comparison",
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  top_row  <- wrap_plots(sp_plots, nrow = 1)
  if (length(scatter_plots) > 0) {
    bot_row <- wrap_plots(scatter_plots, nrow = 1)
    combined <- top_row / bot_row
  } else {
    combined <- top_row
  }

  combined + plot_annotation(
    title = page_title,
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )
}


# ============================================================================
# Section 4: Main orchestration function
# ============================================================================

#' Run spatial correction across all traits (and benches)
#'
#' @param fn          Path to the data file (.csv, .rda, .RData)
#' @param trait_cols  Character vector of phenotype column names
#' @param bench_col   Name of bench/field column; NULL = whole dataset is one bench
#' @param gt_col      Name of genotype column (default "Genotype")
#' @param row_col     Name of row coordinate column (default "row")
#' @param col_col     Name of column coordinate column (default "col")
#' @param method      "SpATS", "mgcv", "mgcv_ps", "mgcv_ps_re", "mgcv_highk",
#'   "mgcv_ad", "sommer", "all", or "all_mgcv"
#' @param output_type "BLUEs", "BLUPs", or "both"
#' @param outlier_iqr Logical; replace IQR outliers with NA before modelling
#' @param output_dir  Directory for output files (created if needed)
#' @param rda_object  Object name when fn is an .rda file (NULL = auto-detect)
#' @return Invisibly, a named list of all per-bench per-trait results
run_spatial_correction <- function(fn,
                                    trait_cols,
                                    bench_col   = NULL,
                                    gt_col      = "Genotype",
                                    row_col     = "row",
                                    col_col     = "col",
                                    method      = "SpATS",
                                    output_type = "BLUEs",
                                    outlier_iqr = TRUE,
                                    output_dir  = "output",
                                    rda_object  = NULL) {

  # ------------------------------------------------------------------
  # 0. Setup
  # ------------------------------------------------------------------
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  methods_to_run <- if (method == "all") {
    c("SpATS", "mgcv", "sommer")
  } else if (method == "all_mgcv") {
    c("SpATS", "mgcv", "mgcv_ps", "mgcv_ps_re", "mgcv_highk", "mgcv_ad", "sommer")
  } else {
    method
  }

  # ------------------------------------------------------------------
  # 1. Load data
  # ------------------------------------------------------------------
  dfr <- read_input(fn, rda_object)

  # Validate required columns
  required <- c(gt_col, row_col, col_col, trait_cols)
  if (!is.null(bench_col)) required <- c(required, bench_col)
  missing_cols <- setdiff(required, names(dfr))
  if (length(missing_cols) > 0)
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))

  # Convert coordinates to numeric
  dfr[[row_col]] <- as.numeric(dfr[[row_col]])
  dfr[[col_col]] <- as.numeric(dfr[[col_col]])

  # Identify benches
  if (is.null(bench_col)) {
    dfr$.bench_internal <- "AllData"
    bench_col_use       <- ".bench_internal"
  } else {
    bench_col_use <- bench_col
  }
  benches <- sort(unique(dfr[[bench_col_use]]))

  # ------------------------------------------------------------------
  # 2. Outlier detection and reporting
  # ------------------------------------------------------------------
  report_path <- file.path(output_dir, "outlier_report.txt")
  rcon        <- file(report_path, "w")
  writeLines(paste0("Outlier report  |  ", Sys.time()), rcon)
  writeLines(paste0("File: ", fn), rcon)
  writeLines("", rcon)

  for (tr in trait_cols) {
    if (outlier_iqr) {
      cat(sprintf("Checking outliers: %s\n", tr))
      out <- replace_outliers_with_na(dfr[[tr]])
      dfr[[tr]] <- out$x
      writeLines(paste0("Trait: ", tr), rcon)
      writeLines(paste0("  Outliers replaced: ", out$n_outliers), rcon)
      if (out$n_outliers > 0)
        writeLines(paste0("  Values: ",
                           paste(round(out$outlier_values, 3), collapse = ", ")),
                   rcon)
      writeLines("", rcon)
    }
  }
  close(rcon)

  # ------------------------------------------------------------------
  # 3. Prepare output containers
  # ------------------------------------------------------------------
  # fitted_values.csv: start from full original data frame
  fitted_out <- dfr

  # BLUEs / BLUPs: accumulated across benches and traits
  blues_out <- NULL
  blups_out <- NULL

  # Nested results list: results[[bench]][[trait]][[method]]
  all_results <- list()

  # ------------------------------------------------------------------
  # 4. Open PDF
  # ------------------------------------------------------------------
  pdf_path <- file.path(output_dir, "spatial_diagnostics.pdf")
  grDevices::pdf(pdf_path, width = 14, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)

  # ------------------------------------------------------------------
  # 5. Main loop: bench × trait × method
  # ------------------------------------------------------------------
  for (bn in benches) {
    bench_data_full <- dfr[dfr[[bench_col_use]] == bn, , drop = FALSE]
    bench_label     <- if (bench_col_use == ".bench_internal") "" else
                         paste0("Bench: ", bn)
    cat(sprintf("\n=== %s ===\n", if (bench_label == "") "Processing data" else bench_label))

    all_results[[as.character(bn)]] <- list()

    for (tr in trait_cols) {
      cat(sprintf("  Trait: %s\n", tr))
      all_results[[as.character(bn)]][[tr]] <- list()

      # Working data for this bench: only rows with non-NA pheno
      bench_data <- bench_data_full
      bench_data[[tr]] <- as.numeric(bench_data[[tr]])

      # Track which rows have valid phenotype values so model outputs
      # (which exclude NA rows) can be written back to the correct positions
      complete_mask <- !is.na(bench_data[[tr]])
      bench_data_complete <- bench_data[complete_mask, , drop = FALSE]

      results_this <- list()

      # ---- Run each method ------------------------------------------
      for (mn in methods_to_run) {
        cat(sprintf("    Method: %s\n", mn))

        # Map method name to function and optional smoother_type
        mgcv_smoother_map <- c(
          mgcv       = "te_default",
          mgcv_ps    = "te_ps",
          mgcv_ps_re = "te_ps_re",
          mgcv_highk = "te_highk",
          mgcv_ad    = "te_ad"
        )

        if (mn %in% names(mgcv_smoother_map)) {
          run_fn       <- run_mgcv_bench
          smoother_arg <- mgcv_smoother_map[[mn]]
        } else {
          run_fn       <- switch(mn,
            SpATS  = run_SpATS_bench,
            sommer = run_sommer_bench
          )
          smoother_arg <- NULL
        }

        extra_args <- if (!is.null(smoother_arg))
          list(smoother_type = smoother_arg) else list()

        r <- tryCatch(
          do.call(run_fn, c(list(bench_data_complete, pheno = tr, gt_col = gt_col,
                                 row_col = row_col, col_col = col_col,
                                 output_type = output_type), extra_args)),
          error = function(e) {
            message("  ERROR in ", mn, " for ", tr, ": ", e$message)
            list(blues = NULL, blups = NULL, fitted = NULL,
                 residuals = NULL, spatial = NULL)
          }
        )

        results_this[[mn]] <- r
        all_results[[as.character(bn)]][[tr]][[mn]] <- r

        # Print diagnostic page for this method
        if (!is.null(r$fitted)) {
          p <- tryCatch(
            plot_diagnostics(bench_data_complete, r, pheno = tr, method = mn,
                             row_col = row_col, col_col = col_col,
                             bench_label = bench_label),
            error = function(e) { message("  plot_diagnostics error: ", e$message); NULL }
          )
          if (!is.null(p)) print(p)
        }
      }  # end method loop

      # ---- Comparison page when running multiple methods --------------
      if (method %in% c("all", "all_mgcv") && length(results_this) > 1) {
        cp <- tryCatch(
          plot_comparison(bench_data_complete, results_this, pheno = tr,
                          row_col = row_col, col_col = col_col,
                          bench_label = bench_label),
          error = function(e) { message("  plot_comparison error: ", e$message); NULL }
        )
        if (!is.null(cp)) print(cp)
      }

      # ---- Accumulate BLUEs / BLUPs across traits and benches -------
      # Determine column name prefix (method tag when method = "all")
      make_cols <- function(df_pred, type, mn) {
        suffix <- if (method %in% c("all", "all_mgcv")) paste0("_", mn) else ""
        val_col <- paste0(tr, "_", type, suffix)
        se_col  <- paste0(tr, "_", type, "_SE", suffix)
        val_name <- if (type == "BLUE") "BLUE" else "BLUP"
        se_name  <- if (type == "BLUE") "BLUE_SE" else "BLUP_SE"
        setNames(df_pred[, c("Genotype", val_name, se_name)],
                 c("Genotype", val_col, se_col))
      }

      for (mn in methods_to_run) {
        r <- results_this[[mn]]

        if (!is.null(r$blues) && output_type %in% c("BLUEs", "both")) {
          add <- make_cols(r$blues, "BLUE", mn)
          names(add)[names(add) == "Genotype"] <- gt_col   # rename key column
          if (!is.null(bench_col) && bench_col != ".bench_internal")
            add[[bench_col]] <- bn
          blues_out <- if (is.null(blues_out)) add else
            merge(blues_out, add,
                  by = intersect(names(blues_out), names(add)), all = TRUE)
        }

        if (!is.null(r$blups) && output_type %in% c("BLUPs", "both")) {
          add <- make_cols(r$blups, "BLUP", mn)
          names(add)[names(add) == "Genotype"] <- gt_col   # rename key column
          if (!is.null(bench_col) && bench_col != ".bench_internal")
            add[[bench_col]] <- bn
          blups_out <- if (is.null(blups_out)) add else
            merge(blups_out, add,
                  by = intersect(names(blups_out), names(add)), all = TRUE)
        }
      }

      # ---- Append fitted/residual/spatial to fitted_out ------------
      # Use the first successful method's values
      primary <- Filter(function(r) !is.null(r$fitted), results_this)
      if (length(primary) > 0) {
        mn_primary <- names(primary)[1]
        r_primary  <- primary[[mn_primary]]
        row_idx    <- which(dfr[[bench_col_use]] == bn)
        # Model outputs exclude NA rows; write back only to complete rows
        complete_idx <- row_idx[complete_mask]

        multi <- method %in% c("all", "all_mgcv")
        suffix <- if (multi) paste0("_", mn_primary) else ""
        fitted_out[row_idx,    paste0(tr, "_Observed")]               <- bench_data_full[[tr]]
        fitted_out[complete_idx, paste0(tr, "_Fitted",       suffix)] <- r_primary$fitted
        fitted_out[complete_idx, paste0(tr, "_Residual",     suffix)] <- r_primary$residuals
        fitted_out[complete_idx, paste0(tr, "_SpatialTrend", suffix)] <- r_primary$spatial

        # If running multiple methods, also add other methods' fitted values
        if (multi) {
          for (mn in methods_to_run[-1]) {
            r2 <- results_this[[mn]]
            if (!is.null(r2$fitted)) {
              fitted_out[complete_idx, paste0(tr, "_Fitted_",       mn)] <- r2$fitted
              fitted_out[complete_idx, paste0(tr, "_Residual_",     mn)] <- r2$residuals
              fitted_out[complete_idx, paste0(tr, "_SpatialTrend_", mn)] <- r2$spatial
            }
          }
        }
      }

    }  # end trait loop
  }  # end bench loop

  # ------------------------------------------------------------------
  # 6. Remove internal bench column if added
  # ------------------------------------------------------------------
  if (bench_col_use == ".bench_internal") {
    fitted_out$.bench_internal <- NULL
    if (!is.null(blues_out) && ".bench_internal" %in% names(blues_out))
      blues_out$.bench_internal <- NULL
    if (!is.null(blups_out) && ".bench_internal" %in% names(blups_out))
      blups_out$.bench_internal <- NULL
  }

  # ------------------------------------------------------------------
  # 7. Write CSV outputs
  # ------------------------------------------------------------------
  write.csv(fitted_out, file.path(output_dir, "fitted_values.csv"),
            row.names = FALSE, na = "")

  if (!is.null(blues_out)) {
    # Reorder: gt_col first, then bench_col if present, then trait columns
    key_cols   <- c(gt_col, bench_col)
    key_cols   <- key_cols[key_cols %in% names(blues_out)]
    trait_cols_out <- setdiff(names(blues_out), key_cols)
    blues_out  <- blues_out[, c(key_cols, trait_cols_out), drop = FALSE]
    write.csv(blues_out, file.path(output_dir, "BLUEs.csv"),
              row.names = FALSE, na = "")
  }

  if (!is.null(blups_out)) {
    key_cols   <- c(gt_col, bench_col)
    key_cols   <- key_cols[key_cols %in% names(blups_out)]
    trait_cols_out <- setdiff(names(blups_out), key_cols)
    blups_out  <- blups_out[, c(key_cols, trait_cols_out), drop = FALSE]
    write.csv(blups_out, file.path(output_dir, "BLUPs.csv"),
              row.names = FALSE, na = "")
  }

  cat(sprintf("\nDone. Outputs written to: %s\n", output_dir))
  cat(sprintf("  spatial_diagnostics.pdf\n"))
  if (!is.null(blues_out)) cat("  BLUEs.csv\n")
  if (!is.null(blups_out)) cat("  BLUPs.csv\n")
  cat("  fitted_values.csv\n")
  cat("  outlier_report.txt\n")

  invisible(all_results)
}
