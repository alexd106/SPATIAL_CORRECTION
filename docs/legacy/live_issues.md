# Live Issues

> Legacy issue log. The canonical issue backlog is now
> [issues.md](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/issues.md).

_Last reviewed: 2026-03-14_

---

## 1. `run_spatial_correction()` does not consistently handle missing phenotypes after outlier removal

### Severity

High

### Summary

`run_spatial_correction()` states that it works on rows with non-missing
phenotypes, but it does not actually subset those rows before fitting or before
writing fitted values back into `fitted_values.csv`.

After IQR-based outlier removal replaces some phenotype values with `NA`, the
downstream model fitters may silently drop rows internally while the assignment
logic still writes the returned fitted values, residuals, and spatial trends
against the full bench slice.

This creates a real risk of row misalignment in:

- `fitted_values.csv`
- diagnostic plots
- method comparison outputs

### Why it matters

If the fitted values returned by a model correspond only to complete cases, but
they are written back to all rows for the bench, the output files can become
silently corrupted.

### References

- [fit_spatial_models.R#L911](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/fit_spatial_models.R#L911)
- [fit_spatial_models.R#L1027](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/fit_spatial_models.R#L1027)

---

## 2. Exported `spatial` trend omits row and column random effects

### Severity

High

### Summary

In the `mgcv` fitters, the exported `spatial` component is taken from the
`te(...)` term only. The row and column random effects are not included, even
though the preferred model explicitly includes them and the documentation
describes them as part of the location-driven correction.

This means that:

- `spatial_trends.csv` is incomplete,
- `spatial_surfaces.png` is incomplete,
- the guide’s interpretation of the estimated spatial landscape is not strictly
  aligned with the implementation.

### Why it matters

Users reading the guide will reasonably infer that the “spatial trend” is the
full environmental correction term. At present it is only the smooth surface,
not the total location effect.

### References

- [spatial_correct_gam.R#L302](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/spatial_correct_gam.R#L302)
- [spatial_correct_gam.R#L329](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/spatial_correct_gam.R#L329)
- [spatial_correct_gam.R#L471](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/spatial_correct_gam.R#L471)
- [spatial_correct_gam_guide.Rmd#L68](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L68)
- [spatial_correct_gam_guide.Rmd#L517](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L517)

---

## 3. Multi-bench API advertises BLUP support that is not implemented

### Severity

Medium

### Summary

The public API for `run_spatial_gam()` advertises:

- `"BLUEs"`
- `"BLUPs"`
- `"both"`

for `estimate_type`.

However, the multi-bench `mgcv` implementation currently supports only BLUEs.
Users can request BLUPs in a joint multi-bench run and will not get the
advertised output.

### Why it matters

This is an API/documentation mismatch. At best it is confusing; at worst it
encourages users to assume a model was fitted that was never actually run.

### References

- [spatial_correct_gam.R#L415](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/spatial_correct_gam.R#L415)
- [spatial_correct_gam.R#L500](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/spatial_correct_gam.R#L500)
- [spatial_correct_gam_guide.Rmd#L282](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L282)
- [spatial_correct_gam_guide.Rmd#L718](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L718)

---

## 4. Outlier detection is global by trait rather than bench-specific

### Severity

Medium

### Summary

`run_spatial_correction()` applies IQR outlier replacement across the full trait
distribution, not within bench.

In a multi-bench design with intentionally heterogeneous spatial severity, a
valid value from a high-intensity bench may be treated as an outlier simply
because it is extreme relative to the pooled distribution.

### Why it matters

This is not just a coding choice; it changes the estimand. It can remove real
bench-specific signal and interacts poorly with the project’s own emphasis on
bench-level heterogeneity.

### References

- [fit_spatial_models.R#L860](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/fit_spatial_models.R#L860)
- [fit_spatial_models.R#L863](/home/alexd106/Repos/SPATIAL_CORRECTION/scripts/fit_spatial_models.R#L863)

---

## 5. Main report is coherent as a results document but not self-contained as a reproducible workflow

### Severity

Medium

### Summary

The main report reads pre-generated CSVs and figures from `../output` and
`figures/`, but it does not run the generation scripts or state clearly that the
artifacts must already exist.

As written, the report assumes prior execution of:

- `scripts/generate_section1_figures.R`
- `scripts/generate_section2_figures.R`
- `scripts/generate_appendix_figures.R`

### Why it matters

The report is logically coherent, but a user starting from a clean checkout
cannot knit it successfully without additional undocumented setup.

### References

- [spatial_correction_report.Rmd#L180](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correction_report.Rmd#L180)
- [spatial_correction_report.Rmd#L291](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correction_report.Rmd#L291)
- [spatial_correction_report.Rmd#L328](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correction_report.Rmd#L328)
- [spatial_correction_report.Rmd#L404](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correction_report.Rmd#L404)

---

## 6. Guide workflow is logically coherent, but path assumptions are fragile

### Severity

Low

### Summary

The guide presents a sensible end-to-end workflow:

1. simulate data,
2. inspect raw spatial bias,
3. fit the joint model,
4. inspect diagnostics,
5. validate against truth.

However, the examples assume the document is knitted from `docs/`:

- scripts are sourced via `../scripts/...`
- generated files are written to `"data"` and `"figures"`
- the prose describes these as `docs/data/` and `docs/figures/`

That is coherent for knitted use, but fragile for readers running chunks
interactively from the repo root.

### Why it matters

The logic is sound, but the execution context is implicit rather than explicit.
This can create unnecessary confusion even when the code itself is correct.

### References

- [spatial_correct_gam_guide.Rmd#L89](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L89)
- [spatial_correct_gam_guide.Rmd#L146](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L146)
- [spatial_correct_gam_guide.Rmd#L159](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd#L159)

---

## 7. General workflow assessment

### Guide

[spatial_correct_gam_guide.Rmd](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correct_gam_guide.Rmd)
is logically coherent and presents a sensible user journey. The main issues are
path fragility and a small amount of drift between documented outputs and the
current implementation.

### Report

[spatial_correction_report.Rmd](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/spatial_correction_report.Rmd)
is coherent as a methods/results document. Its main weakness is reproducibility
from a clean checkout, not conceptual flow.

---

## 8. Open questions

- Should the published “spatial trend” mean only the smooth `te(...)` surface,
  or the full location-driven correction including row and column random effects?
- Is `run_spatial_correction()` intended to be a production entry point, or
  mainly an internal analysis script? This affects how aggressively NA handling,
  outlier behaviour, and output guarantees should be tightened.
