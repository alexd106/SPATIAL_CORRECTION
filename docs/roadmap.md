# Roadmap: Spatial Correction Framework

_Last updated: 2026-03-15_

---

## Purpose

This is the main planning and status document for the project.

It combines:

- project direction,
- current status,
- completed work,
- active priorities,
- near-term implementation order.

Related documents:

- [issues.md](docs/issues.md) — active bugs, risks, and open design questions
- [package_migration.md](docs/package_migration.md) — migration path to a small R package

---

## Goal

Build a flexible spatial-correction framework for field and glasshouse trial
data in R that:

- corrects plot-level spatial bias,
- supports both replicated single-bench and unreplicated multi-bench designs,
- provides defensible genotype BLUEs for downstream analysis,
- compares `SpATS`, `mgcv`, and `sommer`,
- evolves into a maintainable GitHub-hosted R package.

---

## Current Status

### Overall

Phase 1 is complete and the main Phase 2 modelling work is substantially
complete.

### Working now

- standalone `mgcv` workflow in [spatial_correct_gam.R](scripts/spatial_correct_gam.R)
- multi-method comparison framework in [fit_spatial_models.R](scripts/fit_spatial_models.R)
- simulation framework in [simulate_spatial_data.R](scripts/simulate_spatial_data.R)
- report and guide Rmds
- generated figures and output tables for the current validation work

### Current best-supported modelling conclusion

For the scenarios tested so far, the strongest general `mgcv` specification is:

```text
genotype fixed effect
+ per-bench 2D P-spline surface
+ row random effects
+ column random effects
```

For unreplicated within-bench designs, joint multi-bench models are the correct
direction; per-bench BLUEs are not practically useful.

---

## Completed Work

### Phase 1: Per-bench comparison

Completed:

- input reading for `.csv`, `.rda`, `.RData`
- per-bench correction with `SpATS`, `mgcv`, and `sommer`
- multiple `mgcv` smoother variants
- BLUE and BLUP extraction
- visual diagnostics and comparison pages
- real-data validation on `wheatdata`

Key finding:

- the main `mgcv` performance gap was driven by missing row/column random
  effects, not by the spline family itself
- adding row/column random effects closes the gap with `SpATS` almost entirely

### Phase 2: Joint multi-bench models

Completed:

- joint `mgcv` model across benches
- joint `sommer` model using the current column-offset workaround
- simulation validation across multiple spatial patterns
- guide and report documentation for the joint-model workflow

Key finding:

- joint `mgcv` BLUEs strongly outperform averaged per-bench BLUPs on RMSE in the
  tested unreplicated multi-bench scenarios

### Documentation and structure work

Completed:

- functionalised `run_spatial_gam()` workflow
- updated guide for the standalone `mgcv` script
- annotated file structure document
- package migration sketch
- live issue tracking

---

## Active Priorities

### 1. Stabilise the current analysis code

Priority items:

- fix the known output-alignment and API issues listed in
  [issues.md](docs/issues.md)
- tighten the meaning of exported “spatial trend”
- make Rmd workflows more explicit about prerequisites and execution context

### 2. Improve outlier handling

Direction:

- move away from pooled global IQR filtering
- keep deterministic validation checks
- add bench-specific robust raw flagging
- add model-based residual outlier flagging
- default to reporting rather than automatic exclusion

### 3. Extend the joint models

Main next extension:

- genotype-by-environment interaction (GxE)

Additional extension targets:

- partial genotype overlap across benches
- broader real-data validation
- uncertainty reporting and calibration

### 4. Prepare for package migration

Direction:

- move stable functions into package-style modules
- define a narrow public API
- introduce S3 result classes
- separate package code from manuscript/report scripts

---

## Planned Workstreams

### Workstream A: Statistical robustness

Includes:

- outlier handling redesign
- GxE modelling
- partial genotype overlap
- real-data validation

### Workstream B: Code quality and output consistency

Includes:

- fixing current live issues
- clarifying APIs
- improving output semantics
- reducing documentation drift

### Workstream C: Package architecture

Includes:

- modularising code under package-style structure
- documenting exported functions
- adding tests
- converting guides into vignettes

---

## Near-Term Implementation Order

Recommended order:

1. Fix the highest-risk issues in the current scripts.
2. Replace the current outlier logic with a clearer flag/report workflow.
3. Add GxE support to the joint `mgcv` model.
4. Validate on additional real datasets.
5. Begin package-style extraction into `R/`, `tests/`, and `vignettes/`.

---

## Legacy Detail

These older documents still contain useful detail, but are no longer the main
entry points:

- [plan.md](docs/legacy/plan.md) — detailed modelling notes and implementation sketches
- [progress.md](docs/legacy/progress.md) — historical progress log
