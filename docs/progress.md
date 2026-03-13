# Progress Summary: Flexible Multi-Bench Spatial Correction Framework

_Last updated: 2026-03-13_

---

## Status: Phase 2 — Report and guide complete; all spatial types centred; SpATS scale fix applied

---

## What Was Built

### `scripts/SpatialCorrectionFlexible.R` (~1010 lines)

A single self-contained R script implementing a flexible, multi-method spatial correction framework for field trial data. Written entirely fresh — does not source or depend on the original `SpatialCorrectionSpATS.R`.

---

## Functions Implemented

| Function | Purpose |
|---|---|
| `read_input(fn, rda_object)` | Read `.csv`, `.rda`, or `.RData` input |
| `replace_outliers_with_na(x)` | IQR outlier detection and replacement |
| `adaptive_nseg(n_unique)` | Compute nseg as `min(max(5, floor(n/2)), 20)` |
| `empirical_semivariogram(...)` | Manual pairwise semivariance, 15 distance bins |
| `run_SpATS_bench(...)` | SpATS correction; adaptive nseg; BLUEs/BLUPs |
| `run_mgcv_bench(..., smoother_type)` | `gam()` with configurable spatial smoother; BLUEs/BLUPs |
| `run_sommer_bench(...)` | `mmes()` + `spl2Dc()`; spatial as random; BLUEs/BLUPs |
| `run_mgcv_joint(...)` | Joint multi-bench GAM; `te(..., by=bench_f)` + nested row/col RE; BLUEs |
| `run_sommer_joint(...)` | Joint multi-bench `mmes()`; column-offset spl2Dc; BLUEs |
| `plot_heatmap(...)` | ggplot2 viridis tile heatmap |
| `plot_variogram(...)` | Empirical semivariogram line plot |
| `plot_diagnostics(...)` | 6-panel patchwork diagnostic page per method |
| `plot_comparison(...)` | Multi-method comparison page (spatial trends + scatter) |
| `run_spatial_correction(...)` | Main orchestration function (per-bench only) |

---

## Report Structure (`docs/spatial_correction_report.Rmd`)

Two-part HTML report with appendix. Generated via `generate_section1_figures.R`, `generate_section2_figures.R`, and `generate_appendix_figures.R`.

### Part 1 — Per-bench correction with replicated data
- Problem definition and model equations (SpATS, mgcv, sommer)
- Validation on `wheatdata` (330 plots, 107 genotypes, 3 reps)
- SpATS vs mgcv (default) vs sommer comparison
- mgcv variant analysis: gap traced to missing row/col random effects
- `mgcv_ps_re` closes gap entirely (r = 0.9996 with SpATS)

### Part 2 — Joint multi-bench models for unreplicated designs
- Why per-bench BLUEs fail with unreplicated data (model saturation)
- Per-bench SpATS BLUPs: shrinkage problem demonstrated
- Joint model specifications: `run_mgcv_joint()` and `run_sommer_joint()`
- Validation: 4 benches, type 6 edge effects, intensities 15/25/40/60
- Joint mgcv BLUEs: r = 0.986, RMSE = 2.1 (vs SpATS BLUPs RMSE = 6.9)
- Sommer BLUE upward shift explained: column-offset boundary artefact

### Appendix — Additional spatial pattern tests
- Row gradient (type 2), column gradient (type 3), localised hotspot (type 7)
- Same 4-bench unreplicated design, intensities 15/25/40/60
- Joint mgcv BLUEs RMSE: ~2.1 (gradients), ~2.3 (hotspot)
- Confirms joint mgcv generalises across spatial pattern types

---

## `simulate_spatial_data.R` — Spatial Types

| Type | Pattern | Centred |
|---|---|---|
| 1 | No spatial effect | — |
| 2 | Row gradient (half-sine) | Yes (as of 2026-03-13) |
| 3 | Column gradient (half-sine) | Yes (as of 2026-03-13) |
| 4 | Row + column gradient (additive) | Yes (as of 2026-03-13) |
| 5 | Central Gaussian patch | Yes (as of 2026-03-13) |
| 6 | Edge effects with diagonal asymmetry | Yes |
| 7 | Localised hotspot (off-centre Gaussian, ~20% of pots) | Yes |

**Important:** All spatial types used with joint models must be centred to mean zero. Non-centred patterns cause the spatial mean of the reference bench to leak into genotype BLUEs (constant upward shift). All types 2–7 are now centred. Type 1 returns all zeros and needs no centring.

---

## Key Technical Findings

### 1. mgcv performance gap (Part 1)
Default mgcv `te(row, col)` residual SD ~63 vs SpATS ~34. Gap driven entirely by missing row/col random effects. Adding `s(row_f, bs="re") + s(col_f, bs="re")` to a P-spline tensor product (`mgcv_ps_re`) closes the gap: r = 0.9996 with SpATS BLUEs.

### 2. Per-bench BLUEs fail for unreplicated data (Part 2)
With 1 obs per genotype per bench, treating genotype as fixed saturates the model. SpATS/sommer fail outright; mgcv fits but spatial smooth is suppressed. Only BLUPs (random genotype) are viable per-bench — but shrinkage compresses the genotype range.

### 3. Joint models recover proper BLUEs (Part 2)
`run_mgcv_joint()` and `run_sommer_joint()` fit a single model across all benches. Each genotype gets 4 observations (one per bench) — enough to estimate fixed effects. RMSE ~2.1 for mgcv BLUEs vs ~6.9 for SpATS BLUP averages.

### 4. Sommer column-offset boundary artefact (Part 2)
Even with centred spatial patterns, sommer joint BLUEs are systematically shifted upward. Root cause: B-spline basis functions at bench column boundaries are penalised toward zero in the large empty gap regions between benches. Type 6 edge effects are highest at column borders where this penalty is strongest → systematic underestimation → spatial mean leaks into bench fixed effects → all BLUEs inflated. Rankings preserved (r ≈ 0.985) but RMSE inflated (~5.8 vs mgcv ~2.1). A native `vsm(dsm(bench), spl2Dc(...))` would fix this; not yet supported in sommer 4.4.4.

### 5. Spatial pattern centring requirement
If a spatial pattern has a non-zero mean within a bench, and the joint mgcv model uses a zero-mean-constrained smooth (`te(..., by=bench_f)`), the mean of the reference bench's spatial pattern leaks into genotype BLUEs as a constant upward offset. **Fix:** ensure all spatial patterns are centred to mean zero per bench (`raw - mean(raw)`). Applied to all types 2–7 in `simulate_spatial_data.R`. Type 5 (central Gaussian patch) initially appeared near-zero-mean but was not — for default `spatial_scale=1` the tails at grid edges retain ~20% of peak value, giving mean ~30–40% of `spatial_intensity`. Now centred explicitly.

### 6. SpATS PSANOVA spatial surface non-zero mean
SpATS's PSANOVA decomposition constrains each component (row, column, interaction) independently, but the combined surface is not zero-mean over the observed data. `predict(model, which = gt_col)` evaluates the spatial surface at zero, not its mean, so BLUEs/BLUPs are systematically below the phenotype scale by ~`mean(spatial)`. On `wheatdata` this was ~64 units. **Fix:** in `run_SpATS_bench()`, centre the spatial surface and add `mean(spatial)` back to all BLUEs/BLUPs. This ensures all methods are on the same absolute scale. After fix: SpATS BLUPs RMSE improved from 6.85 → 6.15 (edge effects simulation).

---

## Bugs Fixed

| Bug | Fix |
|---|---|
| sommer `m$b` = treatment contrasts, not absolute BLUEs | Recover via `model.matrix() %*% m$b` |
| sommer `rand_base` started with `~` → double-tilde formula | Remove leading `~` from `rand_base` |
| sommer `m$fitted` = Xb only (not Xb+Zu) | Use `fitted(m)` for full predictions |
| sommer 4.4.4 lacks `vsm(dsm(), spl2Dc())` | Column-offset trick for per-bench surfaces |
| `gt_col == "Genotype"` deleted key column in output | Use `names(add)[...] <- gt_col` |
| Non-centred spatial types 2, 3, 4, 5 caused mgcv BLUE upward shift | Centre all spatial types to mean zero (`raw - mean(raw)`) — now applied to types 2–7 |
| SpATS PSANOVA surface has non-zero mean → BLUEs/BLUPs ~64 units below phenotype scale on wheatdata | Centre spatial surface; add mean(spatial) back to BLUEs/BLUPs in `run_SpATS_bench()` |

---

## Key Results Summary

| Scenario | Approach | r | RMSE |
|---|---|---|---|
| wheatdata (replicated, BLUEs) | SpATS | — | res SD 34.1 |
| wheatdata (replicated, BLUEs) | mgcv default | — | res SD 63.0 |
| wheatdata (replicated, BLUEs) | mgcv_ps_re | — | res SD 34.3 |
| 4-bench unreplicated, edge effects | Uncorrected | 0.950 | 3.96 |
| 4-bench unreplicated, edge effects | SpATS BLUPs (avg) | 0.982 | 6.15 |
| 4-bench unreplicated, edge effects | Joint mgcv BLUEs | 0.986 | 2.11 |
| 4-bench unreplicated, edge effects | Joint sommer BLUEs | 0.985 | 5.76 |
| 4-bench unreplicated, row gradient | SpATS BLUPs (avg) | 0.985 | 6.04 |
| 4-bench unreplicated, row gradient | Joint mgcv BLUEs | 0.986 | 2.07 |
| 4-bench unreplicated, col gradient | SpATS BLUPs (avg) | 0.985 | 6.03 |
| 4-bench unreplicated, col gradient | Joint mgcv BLUEs | 0.987 | 2.04 |
| 4-bench unreplicated, hotspot | SpATS BLUPs (avg) | 0.981 | 6.26 |
| 4-bench unreplicated, hotspot | Joint mgcv BLUEs | 0.984 | 2.26 |

---

## File Structure

```
scripts/
  SpatialCorrectionFlexible.R       — main framework
  simulate_spatial_data.R           — simulation (types 1-7, all centred)
  generate_section1_figures.R       — Part 1 figures/data
  generate_section2_figures.R       — Part 2 figures/data
  generate_appendix_figures.R       — Appendix figures/data

docs/
  spatial_correction_report.Rmd     — main report
  spatial_correction_report.html
  simulate_spatial_data_guide.Rmd   — usage guide for simulate_spatial_data.R
  simulate_spatial_data_guide.html
  figures/
    s1_*.png                        — Part 1 figures
    s2_*.png                        — Part 2 figures
    sa_*.png                        — Appendix figures

output/
  section1/                         — wheatdata BLUEs, mgcv variant summary
  section2/                         — comparison_summary.csv, residual_summary.csv
  appendix/                         — appendix_summary.csv
```

---

## Remaining Work

| Item | Description |
|---|---|
| Integration into `run_spatial_correction()` | Joint functions currently called directly; need orchestration with PDF/CSV output |
| Partial genotype overlap | Test information-sharing when genotype sets differ across benches |
| GxE modelling | Allow genotype-by-environment interaction in joint models |
| Real data validation | Validate joint models on actual multi-bench field trial data |
