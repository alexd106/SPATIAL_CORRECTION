# SPATIAL_CORRECTION Project Memory

## Project Overview
Field trial spatial correction framework for R. Working directory: `/home/alexd106/Repos/SPATIAL_CORRECTION`.

## Key Files
- `scripts/fit_spatial_models.R` — main script: per-bench + joint multi-bench functions (renamed from SpatialCorrectionFlexible.R 2026-03-14)
- `scripts/simulate_spatial_data.R` — simulation framework (spatial types 1-7, per-bench params)
- `scripts/generate_section1_figures.R` — generates figures/data for report Part 1 (wheatdata)
- `scripts/generate_section2_figures.R` — generates figures/data for report Part 2 (joint models)
- `docs/spatial_correction_report.Rmd` — main report (two-part structure)
- `docs/simulate_spatial_data_guide.Rmd` — usage guide for simulate_spatial_data.R
- `data/wheatdata.rda` — real test data: 330 obs, 107 geno, 3 rep, 22 rows x 15 cols

## Report Structure (as of 2026-03-13)
- **Part 1:** Per-bench correction on replicated wheatdata. Compares SpATS/mgcv/sommer; mgcv_ps_re closes the gap (r=0.9996 with SpATS). Summary table replaces old 7-method plot.
- **Part 2:** Joint multi-bench models for unreplicated designs. Type 6 edge effects (centred), varying intensity. Highlights SpATS BLUP shrinkage limitation. Joint mgcv/sommer BLUEs recover true effects (RMSE ~2-3 vs 6 for BLUPs).

## Output/Figure Structure
- `docs/figures/s1_*.png` — Part 1 figures
- `docs/figures/s2_*.png` — Part 2 figures
- `docs/figures/sa_*.png` — Appendix figures
- `docs/figures/heatmap_*.png` — simulate_spatial_data_guide.Rmd figures (moved from docs/data/ 2026-03-14)
- `output/section1/` — wheatdata BLUEs, mgcv variant summary CSV
- `output/section2/` — comparison_summary.csv, residual_summary.csv

## Rmd Conventions (as of 2026-03-14)
Both Rmd files share the same YAML/setup structure: toc_depth=4, number_sections=true, fig dimensions in YAML, code_folding=hide, echo=FALSE, fig.align="center", out.width="100%". Both have a `packages` chunk after setup that uses `requireNamespace()` + `tryCatch(install.packages(...))` to silently check/install dependencies without halting on failure.

## Spatial Type 6 (Edge Effects)
- Added 2026-03-13: elevated at borders with diagonal asymmetry (top-left high, bottom-right low)
- Centred to mean zero per bench (avoids sommer BLUE offset from non-zero spatial mean)
- Controlled by spatial_scale (1=linear falloff, >1=sharper edge band)

## fit_spatial_models.R — Architecture
Per-bench: SpATS, mgcv (gam+te()), sommer (mmes+spl2Dc)
Joint multi-bench: `run_mgcv_joint()`, `run_sommer_joint()`
Main orchestrator: `run_spatial_correction()` (per-bench only; joint functions called directly)

## Critical sommer API Notes
- Use `mmes()` NOT `mmer()` — spl2Dc() only works with mmes()
- BLUEs: `b` = treatment contrasts → recover via `X_uniq %*% m$b`
- BLUPs: deviations in `m$uList`; add `m$b[1]` (intercept)
- BLUEs fail for unreplicated data with "Mat::min() has no elements"
- rand_base must NOT start with "~" when constructing BLUPs formula
- Non-zero spatial mean leaks into fixed effects (sommer spatial is random/mean-zero)

## Package Versions
- SpATS 1.0.19, mgcv 1.9.4, sommer 4.4.4

## Phase 2 Status
Joint models implemented. Report refactored. Remaining: partial genotype overlap, GxE, integration into run_spatial_correction().

## Removed Files
- `scripts/SpatialCorrectionSpATS.R` — original buggy script
- `scripts/test4_patch_comparison.R`, `scripts/test5_joint_model.R` — replaced by generate_section*.R
- `output/test1..5/` — replaced by output/section1, output/section2
