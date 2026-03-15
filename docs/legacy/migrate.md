# Migration Plan: From Analysis Scripts to a Small R Package

> Legacy migration note. The canonical migration document is now
> [package_migration.md](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/package_migration.md).

---

## 1. Goal

The long-term goal is to turn the current spatial-correction framework into a
small GitHub-hosted R package that is easy to install, test, document, and
extend.

The package should provide a clean user-facing interface for:

- fitting spatial correction models to trial data,
- comparing supported modelling engines,
- generating diagnostics and summaries,
- simulating trial data for validation,
- later extending the framework to genotype-by-environment (GxE) modelling.

This migration should avoid carrying the current "large sourced script"
architecture directly into the package. R packages work best when functions are
defined under `R/`, documented with `roxygen2`, tested independently, and
exposed through a small, stable API.

---

## 2. Recommendation on File Structure

It is better to split the current large scripts into package-style modules, but
not into one file per tiny helper function.

The recommended approach is:

- split code by responsibility,
- keep public functions few and stable,
- keep engine-specific or helper functions internal unless they are genuinely
  useful to users,
- keep manuscript/report generation scripts outside the package core.

This is preferable to preserving the current layout of a few large scripts such
as `spatial_correct_gam.R` and `fit_spatial_models.R`, because package
maintenance becomes much easier when code is grouped into coherent modules.

---

## 3. Proposed Package Layout

A sensible target structure is:

```text
SPATIAL_CORRECTION/
|
|-- DESCRIPTION
|-- NAMESPACE
|-- README.md
|-- R/
|-- man/
|-- tests/
|   `-- testthat/
|-- vignettes/
|-- data/              # optional packaged example data
|-- inst/
|   `-- extdata/       # optional raw example input files
|-- scripts/           # report generation / manuscript scripts, not package core
`-- docs/              # rendered docs, notes, migration plans
```

Notes:

- `R/` contains package functions only
- `man/` is generated from `roxygen2`
- `tests/testthat/` contains unit and regression tests
- `vignettes/` contains user-facing long-form guides
- `scripts/` should remain available for report generation and reproducibility,
  but should call package functions rather than define them

---

## 4. Proposed `R/` Module Breakdown

Functions should be grouped by concern.

Recommended module layout:

- `R/aaa-package.R`
  - package-level documentation
  - imports and global notes

- `R/input.R`
  - input readers
  - column validation
  - type coercion
  - missing-data checks

- `R/utils-spatial.R`
  - `adaptive_nseg()`
  - semivariogram helpers
  - coordinate and grid utilities

- `R/utils-plots.R`
  - `plot_heatmap()`
  - `plot_variogram()`
  - shared plotting helpers

- `R/model-mgcv-single.R`
  - single-bench `mgcv` fitting functions

- `R/model-mgcv-joint.R`
  - joint `mgcv` fitting functions
  - later GxE extensions

- `R/model-spats.R`
  - SpATS wrappers

- `R/model-sommer.R`
  - sommer wrappers

- `R/diagnostics.R`
  - diagnostic plot builders
  - residual summaries
  - model-check helpers

- `R/simulate.R`
  - simulation functions
  - truth-generation helpers

- `R/compare.R`
  - method comparison wrappers
  - performance summaries

- `R/run.R`
  - high-level user-facing entry points

- `R/tidy.R`
  - S3 methods such as `print()`, `summary()`, `plot()`, `coef()`,
    `fitted()`, and `residuals()`

This is a better compromise than either:

- one huge script per method, or
- one file per helper function

---

## 5. Recommended Public API

The exported API should be intentionally small.

Suggested exported functions:

- `spatial_correct()`
- `compare_spatial_methods()`
- `simulate_field_trial()`
- `plot_diagnostics()`
- `plot_spatial_surface()`

Possible lower-level engine-specific functions can remain internal initially,
for example:

- `fit_mgcv_bench()`
- `fit_mgcv_joint()`
- `fit_spats_bench()`
- `fit_sommer_bench()`
- `fit_sommer_joint()`

These can be exported later only if there is a clear user need.

The preferred user workflow is a stable front door such as:

```r
fit <- spatial_correct(
  data = dfr,
  pheno = "yield",
  genotype = "Genotype",
  row = "Row",
  col = "Col",
  bench = "Bench",
  method = "mgcv_joint"
)
```

rather than requiring users to call engine-specific internal fitters directly.

---

## 6. Classes and Methods

The package should use **S3 classes**, not S4 or R6.

Reasoning:

- S3 is idiomatic for lightweight statistical package objects in R
- the current use case does not need formal S4 slot machinery
- R6 adds unnecessary object-oriented complexity for this workflow

Recommended result classes:

- `spatial_fit`
  - one fitted model result

- `spatial_compare`
  - multi-method comparison result

- `spatial_simulation`
  - simulated dataset plus truth

Recommended S3 methods for `spatial_fit`:

- `print.spatial_fit()`
- `summary.spatial_fit()`
- `plot.spatial_fit()`
- `coef.spatial_fit()`
- `fitted.spatial_fit()`
- `residuals.spatial_fit()`

Recommended S3 methods for `spatial_compare`:

- `print.spatial_compare()`
- `summary.spatial_compare()`
- `plot.spatial_compare()`

Recommended S3 methods for `spatial_simulation`:

- `print.spatial_simulation()`
- `summary.spatial_simulation()`
- `plot.spatial_simulation()`

This gives users a standard R workflow:

```r
fit <- spatial_correct(...)
summary(fit)
plot(fit)
coef(fit)
residuals(fit)
```

---

## 7. Suggested Structure of a `spatial_fit` Object

A `spatial_fit` object should be a structured list with a class attribute.

Suggested components:

- `call`
- `method`
- `model_type`
- `input_spec`
- `data`
- `blues`
- `blups`
- `spatial_trends`
- `fitted_values`
- `residuals`
- `metrics`
- `diagnostics`
- `engine_fit`

Where:

- `engine_fit` is the underlying `gam`, `SpATS`, or `mmes` object for advanced users
- `metrics` stores summary statistics such as RMSE, residual SD, correlations,
  EDF summaries, or information criteria
- `diagnostics` stores objects needed by `summary()` and `plot()`

This is preferable to returning an unstructured list with ad hoc naming.

---

## 8. Method Dispatch Design

The package should not force users into different top-level functions for each
engine unless this becomes necessary.

Preferred interface:

```r
spatial_correct(..., method = "mgcv_joint")
```

or:

```r
spatial_correct(..., engine = "mgcv", mode = "joint")
```

Either is acceptable, but one front door is better than multiple unrelated
entry points because:

- the API stays stable if internal implementations change,
- users learn one interface,
- adding future modes such as GxE is easier,
- testing and documentation are simpler.

---

## 9. What Should Stay Outside the Package Core

Not everything in the repository belongs inside the package.

These should remain outside the package core:

- manuscript/report generation scripts
- one-off figure assembly scripts
- rendered HTML documents
- exploratory notebooks
- large generated output files

Those assets can remain in:

- `scripts/`
- `docs/`
- `output/`

but they should consume the package API rather than define core logic.

---

## 10. Migration Plan

### Phase A: Stabilise the intended API

Decide which functions users should call directly.

Recommended initial user-facing API:

- `spatial_correct()`
- `compare_spatial_methods()`
- `simulate_field_trial()`

Everything else should be treated as internal until a clear reason emerges to
export it.

### Phase B: Extract stable code into `R/`

Move reusable functions out of:

- `scripts/spatial_correct_gam.R`
- `scripts/fit_spatial_models.R`
- `scripts/simulate_spatial_data.R`

and into grouped files under `R/`.

During this phase:

- remove `source("scripts/...")` dependencies from core logic,
- replace script-level side effects with functions,
- keep report scripts working by calling the new package functions.

### Phase C: Introduce classes

Replace loose list returns with structured objects such as:

```r
structure(
  list(...),
  class = "spatial_fit"
)
```

Then implement `print()`, `summary()`, and `plot()` methods.

### Phase D: Add documentation

Use `roxygen2` to document:

- exported functions,
- object classes,
- key arguments,
- return structures,
- worked examples using small data

Generate `NAMESPACE` and `man/` from the docs.

### Phase E: Add tests

Create `testthat` coverage for:

- input validation
- simulation behaviour
- model return structure
- class methods
- regression tests for known benchmark results

Focus especially on deterministic tests for:

- output columns
- dimensions
- presence of diagnostics
- basic recovery of simulated truth within tolerance

### Phase F: Convert guides into vignettes

Likely vignette set:

- getting started with `spatial_correct()`
- comparing `SpATS`, `mgcv`, and `sommer`
- simulating field-trial spatial data
- later: adding GxE modelling

The current long-form materials in `docs/` can be adapted into vignettes while
keeping rendered reports for the manuscript or project website.

### Phase G: Separate package from analysis layer

By the end of the migration:

- package logic lives in `R/`
- tests live in `tests/testthat/`
- vignettes live in `vignettes/`
- report/manuscript scripts live in `scripts/`
- rendered outputs remain in `docs/` and `output/`

This preserves reproducibility without mixing package code with analysis code.

---

## 11. What Not To Do

Avoid the following:

- keeping a package architecture based on sourcing large scripts
- exporting every helper function
- splitting into one file per tiny helper without a module structure
- adopting S4 unless formal slot validation becomes essential
- adopting R6 for what is fundamentally a statistical fitting workflow

The package should stay lightweight, idiomatic, and easy to extend.

---

## 12. Immediate Next Steps

The most practical next steps are:

1. Create a package skeleton with `DESCRIPTION`, `R/`, `tests/testthat/`, and
   `vignettes/`
2. Move the stable `mgcv` helpers and fitters into `R/`
3. Define the `spatial_fit` S3 class
4. Wrap the current recommended workflow in `spatial_correct()`
5. Convert the current guide into the first vignette

Once that is stable, the `SpATS`, `sommer`, comparison, and GxE layers can be
migrated incrementally.
