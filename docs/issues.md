# Issues

_Last reviewed: 2026-03-15_

---

This is the canonical issue/backlog document for current bugs, risks,
documentation drift, and unresolved design questions.

It supersedes [live_issues.md](docs/legacy/live_issues.md).

---

## High Priority

### ~~1. `run_spatial_correction()` does not consistently handle missing phenotypes after outlier removal~~

**Status: Fixed**

NA rows from outlier removal are now filtered into `bench_data_complete` before
passing to model functions. Write-back uses `complete_idx` (the subset of
`row_idx` with valid phenotypes) so fitted/residual/spatial columns align
correctly. `_Observed` still uses the full `row_idx` to preserve original values
including NAs.

### ~~2. Exported `spatial` trend omits row and column random effects~~

**Status: Fixed**

Both fitters (`fit_mgcv_bench`, `fit_mgcv_joint`) now return `spatial_smooth`
(the `te()` surface only) and `spatial_total` (`te()` + row RE + column RE),
plus `row_re` and `col_re` vectors. `spatial_trends.csv` exports both columns.
The diagnostic heatmap in Panel 2 uses `spatial_smooth` (smooth and
interpretable); a separate `row_col_re_{trait}_{bench}.png` strip plot shows the
row and column RE magnitudes. The guide panel descriptions and output-file table
have been updated to reflect this split. BLUEs are unaffected — they were
computed from the full model throughout.

---

## Medium Priority

### 3. Multi-bench API advertises BLUP support that is not implemented

Summary:

The public API for `run_spatial_gam()` advertises `"BLUEs"`, `"BLUPs"`, and
`"both"` for `estimate_type`, but the joint multi-bench `mgcv` path supports
only BLUEs.

References:

- [spatial_correct_gam.R#L415](scripts/spatial_correct_gam.R#L415)
- [spatial_correct_gam.R#L500](scripts/spatial_correct_gam.R#L500)
- [spatial_correct_gam_guide.Rmd#L282](docs/spatial_correct_gam_guide.Rmd#L282)
- [spatial_correct_gam_guide.Rmd#L718](docs/spatial_correct_gam_guide.Rmd#L718)

### ~~4. Outlier detection is global by trait rather than bench-specific~~

**Status: Fixed**

Outlier detection (1.5 × IQR) is now applied per-bench within the main
bench × trait loop, so fences are computed from each bench's own distribution.
The outlier report now includes bench labels. Single-bench runs are unaffected
(trivially equivalent).

### 5. Main report is coherent as a results document but not self-contained as a reproducible workflow

Summary:

The main report expects pre-generated figures and CSVs but does not make those
prerequisites explicit enough.

References:

- [spatial_correction_report.Rmd#L180](docs/spatial_correction_report.Rmd#L180)
- [spatial_correction_report.Rmd#L291](docs/spatial_correction_report.Rmd#L291)
- [spatial_correction_report.Rmd#L328](docs/spatial_correction_report.Rmd#L328)
- [spatial_correction_report.Rmd#L404](docs/spatial_correction_report.Rmd#L404)

---

## Low Priority

### 6. Guide workflow is logically coherent, but path assumptions are fragile

Summary:

The guide assumes knitting from `docs/`, which is fine for rendering but
slightly fragile for readers running chunks interactively from the repo root.

References:

- [spatial_correct_gam_guide.Rmd#L89](docs/spatial_correct_gam_guide.Rmd#L89)
- [spatial_correct_gam_guide.Rmd#L146](docs/spatial_correct_gam_guide.Rmd#L146)
- [spatial_correct_gam_guide.Rmd#L159](docs/spatial_correct_gam_guide.Rmd#L159)

---

## Open Design Questions

- ~~Should exported “spatial trend” mean only the smooth `te(...)` surface, or the
  full location-driven correction including row and column random effects?~~
  **Resolved:** export both as `spatial_smooth` and `spatial_total` (see issue 2).
- Is `run_spatial_correction()` intended to be a production entry point or an
  internal analysis script?
