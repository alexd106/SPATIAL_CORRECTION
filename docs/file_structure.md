# Project File Structure

_Last updated: 2026-03-14_

---

```
SPATIAL_CORRECTION/
|
|-- data/
|   `-- wheatdata.rda              Real replicated field trial: 330 plots,
|                                  107 genotypes, 3 reps, 22 rows x 15 cols.
|                                  Used in Part 1 of the main report.
|
|-- scripts/
|   |-- spatial_correct_gam.R      Standalone mgcv spatial correction script.
|   |                              Defines run_spatial_gam() plus all helper and
|   |                              fitting functions. Source this file to use the
|   |                              functions interactively, or call run_spatial_gam()
|   |                              directly to run the full pipeline.
|   |
|   |-- fit_spatial_models.R       Multi-method comparison framework. Implements
|   |                              SpATS, mgcv (multiple smoother variants), and
|   |                              sommer per-bench and joint multi-bench corrections.
|   |                              Sources spatial_correct_gam.R for shared helpers.
|   |                              Main entry point: run_spatial_correction().
|   |
|   |-- simulate_spatial_data.R    Simulation framework for synthetic field trials.
|   |                              Supports 7 spatial pattern types (none, row/col
|   |                              gradients, Gaussian patch, edge effects, hotspot).
|   |                              All patterns centred to mean zero per bench.
|   |                              Main function: simulate_field_trial().
|   |
|   |-- generate_section1_figures.R  Generates figures and CSV outputs for Part 1
|   |                                of the main report (wheatdata, per-bench
|   |                                correction, method comparison).
|   |
|   |-- generate_section2_figures.R  Generates figures and CSV outputs for Part 2
|   |                                of the main report (joint multi-bench models,
|   |                                unreplicated design, method validation).
|   |
|   `-- generate_appendix_figures.R  Generates figures for the report appendix
|                                    (additional spatial pattern types: row gradient,
|                                    col gradient, localised hotspot).
|
|-- docs/
|   |-- spatial_correction_report.Rmd   Main technical report (two parts + appendix).
|   |                                   Part 1: per-bench correction on replicated
|   |                                   wheatdata. Part 2: joint multi-bench models
|   |                                   for unreplicated designs.
|   |-- spatial_correction_report.html  Compiled HTML version of the main report.
|   |
|   |-- spatial_correct_gam_guide.Rmd   Step-by-step user guide for spatial_correct_gam.R.
|   |                                   Walks through simulate -> run_spatial_gam() ->
|   |                                   interpret diagnostics -> validate, using a
|   |                                   4-bench simulated BNI dataset.
|   |-- spatial_correct_gam_guide.html  Compiled HTML version of the guide.
|   |
|   |-- simulate_spatial_data_guide.Rmd  Usage guide for simulate_spatial_data.R.
|   |                                    Documents all spatial types, parameters,
|   |                                    and plotting functions.
|   |-- simulate_spatial_data_guide.html Compiled HTML version.
|   |
|   |-- file_structure.md           This file.
|   |-- plan.md                     Original project planning notes.
|   |-- progress.md                 Running log of development decisions, technical
|   |                               findings, bugs fixed, and key results.
|   |
|   |-- data/
|   |   |-- BNI_simulation.csv      Simulated 4-bench BNI trial (960 plots). Generated
|   |   |                           by spatial_correct_gam_guide.Rmd via
|   |   |                           simulate_field_trial(). Input to run_spatial_gam().
|   |   |-- BNI_true_effects.csv    True genotype effects from the simulation. Used in
|   |   |                           the guide to validate correction quality.
|   |   `-- gam_output/             CSVs written by run_spatial_gam() during the guide:
|   |       |-- BLUEs.csv           Spatially corrected genotype estimates.
|   |       |-- spatial_trends.csv  Per-plot estimated spatial contributions.
|   |       `-- model_summary.csv   Residual SD, EDF, convergence per trait/bench.
|   |
|   `-- figures/
|       |-- guide_raw_observed.png        Raw BNI heatmaps (spatial_correct_gam_guide).
|       |-- guide_true_spatial.png        True spatial patterns (spatial_correct_gam_guide).
|       |-- guide_estimated_spatial.png   Estimated spatial surfaces (spatial_correct_gam_guide).
|       |-- guide_diagnostics_BNI_*.png   6-panel diagnostics per bench (spatial_correct_gam_guide).
|       |-- heatmap_*.png                 Spatial type illustrations (simulate_spatial_data_guide).
|       |-- s1_*.png                      Part 1 figures (wheatdata method comparison).
|       |-- s2_*.png                      Part 2 figures (joint model validation).
|       `-- sa_*.png                      Appendix figures (additional pattern types).
|
`-- output/
    |-- gam/                        Outputs from run_spatial_gam() on wheatdata (single-bench).
    |   |-- BLUEs.csv
    |   |-- spatial_trends.csv
    |   |-- model_summary.csv
    |   |-- diagnostics_yield_single.png
    |   |-- blues_distribution.png
    |   `-- spatial_surfaces.png
    |
    |-- section1/                   Outputs from generate_section1_figures.R.
    |   |-- BLUEs.csv               Wheatdata BLUEs from all methods.
    |   `-- mgcv_variants_summary.csv  Residual SD comparison across mgcv smoother variants.
    |
    |-- section2/                   Outputs from generate_section2_figures.R.
    |   |-- comparison_summary.csv  r and RMSE for each method across scenarios.
    |   `-- residual_summary.csv    Per-bench residual statistics.
    |
    `-- appendix/                   Outputs from generate_appendix_figures.R.
        `-- appendix_summary.csv    r and RMSE for appendix spatial pattern types.
```

---

## Key relationships

- `spatial_correct_gam.R` is the primary user-facing script. Source it to get
  `run_spatial_gam()` and its helper functions.
- `fit_spatial_models.R` sources `spatial_correct_gam.R` to inherit the shared
  helper functions (`read_input()`, `adaptive_nseg()`, `empirical_semivariogram()`,
  `plot_heatmap()`, `plot_variogram()`, `make_diagnostic_plots()`, `fit_mgcv_bench()`,
  `fit_mgcv_joint()`), then adds SpATS and sommer methods on top.
- The `generate_section*.R` and `generate_appendix_figures.R` scripts source
  `fit_spatial_models.R` (which in turn sources `spatial_correct_gam.R`) and
  write figures and CSVs consumed by the Rmd reports.
- The `spatial_correct_gam_guide.Rmd` sources only `spatial_correct_gam.R` and
  `simulate_spatial_data.R`; it does not depend on SpATS or sommer.
