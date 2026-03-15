# Issues

_Last reviewed: 2026-03-15_

---

This is the canonical issue/backlog document for current bugs, risks,
documentation drift, and unresolved design questions.

It supersedes [live_issues.md](docs/legacy/live_issues.md).

---

## High Priority

### 1. `run_spatial_correction()` does not consistently handle missing phenotypes after outlier removal

Summary:

`run_spatial_correction()` does not safely subset or realign rows after
phenotype values are replaced with `NA`. Downstream model outputs can therefore
be written back against the wrong rows.

Impact:

- `fitted_values.csv` may be misaligned
- diagnostics may be wrong
- method comparisons may be wrong

References:

- [fit_spatial_models.R#L911](scripts/fit_spatial_models.R#L911)
- [fit_spatial_models.R#L1027](scripts/fit_spatial_models.R#L1027)

### 2. Exported `spatial` trend omits row and column random effects

Summary:

The exported `spatial` component in the `mgcv` fitters currently reflects only
the `te(...)` smooth, not the full location-driven correction term described in
the docs.

Impact:

- `spatial_trends.csv` is incomplete
- `spatial_surfaces.png` is incomplete
- the guide overstates what the spatial-trend output represents

References:

- [spatial_correct_gam.R#L302](scripts/spatial_correct_gam.R#L302)
- [spatial_correct_gam.R#L329](scripts/spatial_correct_gam.R#L329)
- [spatial_correct_gam.R#L471](scripts/spatial_correct_gam.R#L471)
- [spatial_correct_gam_guide.Rmd#L68](docs/spatial_correct_gam_guide.Rmd#L68)
- [spatial_correct_gam_guide.Rmd#L517](docs/spatial_correct_gam_guide.Rmd#L517)

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

### 4. Outlier detection is global by trait rather than bench-specific

Summary:

The current IQR rule is applied across the pooled trait distribution, which is
not appropriate for heterogeneous multi-bench data.

References:

- [fit_spatial_models.R#L860](scripts/fit_spatial_models.R#L860)
- [fit_spatial_models.R#L863](scripts/fit_spatial_models.R#L863)

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

- Should exported “spatial trend” mean only the smooth `te(...)` surface, or the
  full location-driven correction including row and column random effects?
- Is `run_spatial_correction()` intended to be a production entry point or an
  internal analysis script?
