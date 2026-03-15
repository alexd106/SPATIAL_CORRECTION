# Plan: Flexible Multi-Bench Spatial Correction Framework

> Legacy detailed planning document. The canonical high-level planning/status
> document is now [roadmap.md](/home/alexd106/Repos/SPATIAL_CORRECTION/docs/roadmap.md).
> Keep this file for detailed modelling notes and implementation sketches.

---

## 1. Problem Statement

`SpatialCorrectionSpATS.R` applies spatial correction to field trial data using the SpATS package. As implemented, it has several limitations from a statistical and practical standpoint:

1. **Hardcoded smooth complexity**: `nseg = c(5,10)` fixes the number of B-spline segments regardless of field dimensions. For a small bench (e.g. 10 rows × 30 cols), 10 column segments may over-fit; for a large bench (60 × 60), 5 row segments will grossly under-smooth.
2. **Single bench only**: Must be re-run manually for each field/bench, with no systematic multi-bench support.
3. **BLUPs only**: `genotype.as.random = TRUE` always. BLUPs are shrunken estimates — appropriate for within-trial ranking but biased and not recommended when combining across environments.
4. **Bugs**: Function name mismatch (`SpatialCorrectionSpATS_BLUPs` called but not defined); undefined variables `study` and `year`.

---

## 2. Statistical Justification for a More Flexible Framework

### Why spatial correction matters

Field trial observations follow:

$$y_{ij} = g_i + \xi(r_j, c_j) + \varepsilon_{ij}$$

where $g_i$ is the true genotype effect, $\xi(r_j, c_j)$ is a spatially structured environmental effect (nutrient gradients, soil moisture, drainage patterns, bench edge effects), and $\varepsilon_{ij} \sim N(0,\sigma^2_e)$ is plot-level noise. If $\xi$ is not modelled, genotype estimates are confounded with the spatial pattern, inflating or deflating apparent genetic differences.

### Why the current approach may be inadequate

- P-spline smoothness is controlled jointly by the number of knots (segments) and the penalty. Fixing `nseg = c(5,10)` means the model cannot adapt to fields of varying size or spatial complexity.
- The row and column random effects (`rnf`, `cnf`) partially overlap with the spatial smooth in what they capture — this double-counting can lead to poor identifiability of the spatial surface.
- Modern alternatives (mgcv, sommer) offer automatic smoothness selection (via REML) that adapts to the data without requiring user-specified knot counts.

---

## 3. Model Equations

### 3.1 Current SpATS model

For a single bench with $n$ plots, genotypes $i = 1, \ldots, G$, rows $r$ and columns $c$:

$$y_{rc} = \mu + g_i + f_{\text{SAP}}(c, r) + \rho_r + \kappa_c + \varepsilon_{rc}$$

| Term | Type | Description |
|------|------|-------------|
| $\mu$ | Fixed scalar | Grand mean |
| $g_i$ | Fixed (BLUEs) or random (BLUPs) | Genotype effect |
| $f_{\text{SAP}}(c, r)$ | Fixed smooth (2D P-spline) | Spatial surface via SAP; `nseg` controls basis dimension in each direction |
| $\rho_r$ | Random: $N(0, \sigma^2_r)$ | Residual row trend after spline |
| $\kappa_c$ | Random: $N(0, \sigma^2_c)$ | Residual column trend after spline |
| $\varepsilon_{rc}$ | $N(0, \sigma^2_e)$ | Plot-level residual |

The SAP smoother uses B-splines with a **second-difference penalty**. Smoothness is determined jointly by `nseg` and an estimated penalty parameter (variance ratio). The current script fixes `nseg = c(5,10)` — the penalty is estimated but the basis dimension is fixed.

**BLUPs**: `genotype.as.random = TRUE` → $g_i \sim N(0, \sigma^2_g)$, shrunken toward the mean.

**BLUEs**: `genotype.as.random = FALSE` → $g_i$ treated as a fixed parameter, no shrinkage.

---

### 3.2 mgcv::gam model

$$y_{rc} = \mu + g_i + s(r, c) + \varepsilon_{rc}$$

| Term | Type | Description |
|------|------|-------------|
| $g_i$ | Fixed (BLUEs) or random via `s(geno, bs="re")` (BLUPs) | Genotype effect |
| $s(r, c)$ | Penalised smooth, basis and penalty estimated jointly | 2D spatial surface |

Two smooth specifications to compare:

- **Isotropic**: `s(row, col, bs="tp")` — thin-plate regression spline, assumes equal smoothness in row and column directions.
- **Anisotropic**: `te(row, col)` / `t2(row, col)` — tensor product smooth, allows **different smoothness in row vs column direction**. Preferred for field trials with unequal plot spacing or when gradient scales differ between dimensions.

**Key advantage**: Smoothness parameter $\lambda$ estimated by **REML** — fully automatic, no manual `nseg` required.

**BLUEs**: `gam(pheno ~ 0 + Genotype + te(row, col), method="REML")` — genotype coefficients estimated as fixed effects directly.

**BLUPs**: replace `0 + Genotype` with `s(Genotype, bs="re")` — equivalent to $g_i \sim N(0, \sigma^2_g)$.

**Note on `gamm()`**: Do **not** use `mgcv::gamm()`. It uses `nlme` as its backend, which struggles with high-dimensional crossed random effects (many genotypes) and is numerically less stable. `gam()` is the correct choice here; `bam()` for large datasets.

---

### 3.3 sommer model

$$y_{rc} = \mu + g_i + \mathbf{u}_{rc} + \varepsilon_{rc}$$

The spatial surface is modelled as a **random effect**:

$$\mathbf{u} \sim N(\mathbf{0},\, \sigma^2_s \mathbf{K})$$

where $\mathbf{K}$ is the **2D P-spline kernel matrix** constructed by `spl2Dc(col, row, nsegments=c(nc, nr))` as a Kronecker product of 1D B-spline bases:

$$\mathbf{K} = \mathbf{B}_c \otimes \mathbf{B}_r$$

with a separable second-difference penalty.

**Key conceptual difference from SpATS and mgcv**: In SpATS and mgcv, the spatial smooth is a **fixed** function (penalised but treated as a fixed contribution). In sommer, the spatial surface is a **random effect** — it is shrunk toward zero. This can be advantageous when the spatial pattern is uncertain or weak, but may under-correct when the spatial effect is strong and structured.

**BLUEs**: `mmes(pheno ~ Genotype, random = ~spl2Dc(...))` — genotype fixed, spatial random.

**BLUPs**: `mmes(pheno ~ 1, random = ~Genotype + spl2Dc(...))` — both genotype and spatial are random.

**API note**: `spl2Dc()` requires `sommer::mmes()`, not `mmer()`. BLUEs from `mmes` use treatment contrasts; absolute BLUEs must be recovered via `model.matrix(~ gt_col) %*% m$b`.

---

## 4. Pros and Cons Summary

| Criterion | SpATS (enhanced) | mgcv::gam | sommer |
|-----------|-----------------|-----------|--------|
| **Designed for field trials** | Yes | No (general GAM) | Yes |
| **Smoothness selection** | Semi-auto (fixed `nseg`, penalty estimated) | Fully automatic (REML) | Semi-auto (fixed `nsegments`) |
| **Anisotropic smoothing** | Yes (separate `nseg` per direction) | Yes (`te()`) | Yes (separate segments) |
| **BLUEs** | Yes | Yes | Yes |
| **BLUPs** | Yes | Yes (`bs="re"`) | Yes |
| **Spatial surface type** | Fixed smooth (P-spline) | Fixed smooth (thin-plate or tensor product) | **Random effect** (P-spline) |
| **SE of BLUEs/BLUPs** | Yes (from `predict()`) | Yes (from model `vcov`) | Yes (from `model$Ci`) |
| **Multi-bench joint model** | Per-bench only | Extendable via `by=` or multi-level smooth | Natively extendable |
| **Computational speed** | Fast | Fast | Slower (iterative REML) |
| **Stability / maturity** | High | Very high | Moderate (API evolving) |
| **Unreplicated data** | BLUPs only recommended | BLUEs saturate; BLUPs work | BLUEs fail; BLUPs work |
| **Main limitation** | Fixed `nseg`; no joint multi-bench | No field-trial-specific workflow | Spatial-as-random may under-correct; slower |

---

## 5. Comparison Methodology

### 5.1 With simulated data (known truth)

**Primary metric — genotype estimate accuracy:**

$$r(\hat{g}_i,\, g^*_i) \quad \text{Pearson correlation (higher is better)}$$

$$\text{RMSE}(\hat{g}_i, g^*_i) = \sqrt{\frac{1}{G}\sum_{i=1}^{G}(\hat{g}_i - g^*_i)^2} \quad \text{(lower is better)}$$

**Secondary metric — spatial trend recovery:**

$$r(\hat{\xi}_{rc},\, \xi^*_{rc}) \quad \text{correlation of estimated vs true spatial surface}$$

**Baseline**: Uncorrected per-genotype phenotypic mean. A valid spatial correction must improve $r(\hat{g}, g^*)$ substantially over this baseline.

### 5.2 With real data (unknown truth — wheatdata)

When true effects are unknown:
1. **Moran's I** on residuals: $I \approx 0$ indicates no residual spatial autocorrelation; $I > 0$ indicates under-correction.
2. **Empirical semivariogram** of residuals vs raw values: residual variogram should be flatter (no spatial structure remaining) after correction.
3. **SE of BLUEs**: a more adequate spatial model typically yields smaller standard errors — though very small SE can also signal over-fitting.

### 5.3 Comparison visualisation

For each trait, a method-comparison page in the diagnostic PDF:
- Side-by-side heatmaps of estimated spatial trend (SpATS | mgcv | sommer)
- Side-by-side heatmaps of residuals
- Pairwise scatter plots of BLUEs across methods (high correlation = robustness; discrepancies = investigate)
- (Simulated only) Bar chart or table of $r(\hat{g}, g^*)$ and RMSE per method, with uncorrected baseline

---

## 6. Assessment of `simulate_spatial_data.R` for Phase 2

### Currently fit for purpose
- Generates multi-bench data with configurable spatial patterns (types 1–7, all centred to mean zero)
- Per-bench spatial type and intensity configurable via `spatial_type_per_bench` and `spatial_intensity_per_bench`
- Returns true genotype effects (`$true_effects`) and per-plot spatial contributions (`True_Spatial` column in `$data`) for validation
- Same genotype set across all benches (complete overlap — useful baseline for Phase 2)
- Heatmap functions for visualisation
- Usage guide: `docs/simulate_spatial_data_guide.Rmd`

### Gaps — status

| Gap | Status |
|-----|--------|
| **Per-bench spatial type/intensity** | **RESOLVED** — `spatial_type_per_bench` and `spatial_intensity_per_bench` parameters added |
| **True spatial effect saved in output** | **RESOLVED** — `True_Spatial` column included in `$data` |
| **No partial genotype overlap option** | **OPEN** — genotype sets are still identical across all benches |
| **No GxE option** | **OPEN** — all genotype effects identical across benches |
| **Single trait only** | **OPEN** — multi-trait requires multiple calls and joining |

### Verdict
`simulate_spatial_data.R` now supports **Phase 2** requirements for spatial type/intensity heterogeneity and true spatial recovery. Partial genotype overlap and GxE remain as future extensions.

---

## 7. Phased Implementation Plan

### Phase 1: Per-bench comparison — **COMPLETE**

**File**: `scripts/fit_spatial_models.R`

Written entirely fresh. All method-specific correction functions and all visualisation and utility functions implemented de novo in a single file.

**Implemented:**
- Read CSV (`.csv`) or R binary (`.rda` / `.RData`) input
  - Column names configurable (`gt_col`, `row_col`, `col_col`, `bench_col`)
  - When no `bench_col` is present or `bench_col = NULL`, treat the entire dataset as one bench
- Outlier detection via IQR method (optional, `outlier_iqr = TRUE/FALSE`)
- Per-bench spatial correction using:
  - `method = "SpATS"` — enhanced SpATS with adaptive `nseg`
  - `method = "mgcv"` — `gam()` with `te(row, col)`
  - `method = "sommer"` — `mmes()` with `spl2Dc()`
  - `method = "all"` — run all three and produce comparison outputs
- `output_type = "BLUEs"`, `"BLUPs"`, or `"both"`
- Diagnostic PDF with 6 panels per method per trait
- Output CSVs preserving all original file columns

**Verification results:**

| Test | Dataset | Method | r(ĝ, g*) | Notes |
|------|---------|--------|----------|-------|
| Test 1 | wheatdata (3 rep) | SpATS BLUEs | — | SpATS↔sommer r=0.992; res SD 34.1 |
| Test 1 | wheatdata (3 rep) | mgcv default BLUEs | — | SpATS↔mgcv r=0.871; res SD 63.0 |
| Test 1 | wheatdata (3 rep) | mgcv_ps_re BLUEs | — | r=0.9996 with SpATS; res SD 34.3 |
| Test 1 | wheatdata (3 rep) | sommer BLUEs | — | r=0.992 with SpATS; res SD 33.2 |
| Test 2 | Simulated (1 rep, spatial type 4) | Uncorrected | 0.881 | Baseline |
| Test 2 | Simulated | SpATS BLUPs | **0.958** | +0.077 vs baseline |
| Test 2 | Simulated | mgcv BLUPs | **0.958** | +0.077 vs baseline |
| Test 2 | Simulated | sommer BLUPs | **0.931** | +0.050 vs baseline |

**Note:** After the SpATS spatial surface centring fix, all three methods in Test 1 produce BLUEs on the same absolute phenotype scale (mean ≈ grand mean of yield). Previously SpATS BLUEs were ~64 units below mgcv/sommer due to the PSANOVA non-zero surface mean.

### Phase 2: Multi-bench joint models — **SUBSTANTIALLY COMPLETE**

**Implemented:**
- `run_mgcv_joint()` — joint GAM across all benches: `0 + Genotype + bench_f + te(..., by=bench_f) + s(row_f, bench_f, bs="re") + s(col_f, bench_f, bs="re")`
- `run_sommer_joint()` — joint `mmes()` with column-offset trick for independent per-bench surfaces
- Validated across 4 spatial pattern types: edge effects (type 6), row gradient (type 2), column gradient (type 3), localised hotspot (type 7)
- Joint mgcv consistently achieves r ≈ 0.984–0.987, RMSE ≈ 2.0–2.3 (vs SpATS BLUPs RMSE ≈ 6.0–6.3)
- Documented in `docs/spatial_correction_report.Rmd` (Part 2 + Appendix)

**Remaining for Phase 2:**

| Item | Description |
|------|-------------|
| Integration into `run_spatial_correction()` | Joint functions currently called directly; need orchestration with PDF/CSV output |
| Partial genotype overlap | Test information-sharing when genotype sets differ across benches |
| GxE modelling | Allow genotype-by-environment interaction in joint models |
| Real data validation | Validate joint models on actual multi-bench field trial data |

---

## 8. Genotype-by-Environment (GxE) Extension Plan

### 8.1 Motivation

The current joint models assume that each genotype has a single common effect
across all benches:

$$y_{ijrc} = \mu + g_i + e_j + \xi_j(r, c) + \varepsilon_{ijrc}$$

where $g_i$ is constant across environments (benches). This is a strong
assumption. In real trials, genotypes often respond differently across benches
or environments because of microclimate differences, management variation, or
true biological interaction.

If GxE is present but omitted, the current joint BLUE model may:

- overstate certainty in genotype means,
- bias average genotype estimates when responses are not parallel across benches,
- attribute some genotype-specific environment response to the spatial term or
  residual noise,
- mis-rank unstable genotypes that perform inconsistently across benches.

The next methodological extension should therefore quantify and model
genotype-by-environment interaction explicitly.

### 8.2 Proposed model extension

The recommended first implementation is to add a **random GxE deviation term**
to the current joint `mgcv` model:

$$y_{ijrc} = \mu + g_i + e_j + (ge)_{ij} + \xi_j(r, c) + \rho_{r(j)} + \kappa_{c(j)} + \varepsilon_{ijrc}$$

with

$$ (ge)_{ij} \sim N(0, \sigma^2_{ge}) $$

where:

- $g_i$ is the genotype main effect,
- $e_j$ is the fixed bench/environment effect,
- $(ge)_{ij}$ is the genotype-specific deviation within bench $j$,
- $\xi_j(r, c)$ is the bench-specific 2D spatial surface,
- $\rho_{r(j)}$ and $\kappa_{c(j)}$ are row and column random effects nested
  within bench,
- $\varepsilon_{ijrc} \sim N(0, \sigma^2_e)$ is the residual.

This preserves the current interpretation of genotype BLUEs as average genotype
performance across benches, while allowing bench-specific deviations around that
average.

### 8.3 Recommended `mgcv` implementation

Extend the current `run_mgcv_joint()` model by adding:

```r
s(Genotype, bench_f, bs = "re")
```

to the existing formula:

```r
pheno ~ 0 + Genotype + bench_f +
  te(Row, Col, by = bench_f, bs = c("ps", "ps")) +
  s(row_f, bench_f, bs = "re") +
  s(col_f, bench_f, bs = "re") +
  s(Genotype, bench_f, bs = "re")
```

Interpretation:

- `0 + Genotype` retains fixed genotype effects, i.e. average BLUEs across benches
- `bench_f` estimates mean bench differences
- `te(..., by = bench_f)` keeps an independent spatial surface per bench
- `s(row_f, bench_f, bs = "re")` and `s(col_f, bench_f, bs = "re")` retain
  the row/column random effects already shown to improve fit
- `s(Genotype, bench_f, bs = "re")` captures GxE as a penalised random
  genotype-by-bench deviation

This should be the first implementation because it is simple, regularised, and
compatible with the current `mgcv` framework.

### 8.4 Quantities to estimate and report

The GxE extension should return and document the following quantities:

1. **GxE variance component**
   - Estimate the variance associated with `s(Genotype, bench_f, bs = "re")`
   - Compare it with genotype main-effect variance (when applicable) and
     residual variance

2. **Per-genotype instability**
   - Extract the bench-specific random deviations for each genotype
   - Summarise instability as the SD or RMS of genotype-specific deviations
     across benches

3. **Change in genotype ranking**
   - Compare genotype rankings from:
     - joint model without GxE
     - joint model with GxE
   - Flag genotypes whose ranks or estimated means change materially

4. **Improvement in model fit**
   - Compare REML score, AIC, residual SD, and residual spatial structure
     between the no-GxE and GxE models

5. **Bench-specific adjusted predictions**
   - Optionally return genotype-by-bench adjusted predictions in addition to the
     across-bench BLUEs

### 8.5 Simulation extension

The simulation framework should be extended so genotype effects are decomposed
into:

$$g_{ij} = g_i + (ge)_{ij}$$

where $(ge)_{ij} \sim N(0, \sigma^2_{ge})$.

Implementation sketch in `simulate_spatial_data.R`:

```r
ge_dev <- matrix(rnorm(n_geno * n_bench, 0, sd_ge),
                 nrow = n_geno, ncol = n_bench)
```

and for each bench $j$:

```r
phenotype = geno_true[i] + ge_dev[i, j] + spatial_effect + error
```

Add a user-facing argument such as:

- `sd_ge = 0` by default
- `sd_ge > 0` to introduce controllable GxE magnitude

Validation scenarios should include:

- `sd_ge = 0` to confirm the GxE model does not overfit badly when no
  interaction exists
- small GxE
- moderate GxE
- strong GxE

### 8.6 Development plan

#### Step 1: Extend `run_mgcv_joint()`

- Add argument `allow_gxe = FALSE/TRUE`
- When `TRUE`, append `s(Genotype, bench_f, bs = "re")` to the formula
- Return a `gxe_summary` object containing:
  - estimated smoothing parameter / variance proxy,
  - per-genotype instability summaries,
  - optional genotype-by-bench adjusted deviations

#### Step 2: Extend the simulation script

- Add argument `sd_ge`
- Simulate genotype-specific bench deviations
- Save true GxE values alongside the existing true genotype and true spatial
  outputs

#### Step 3: Add validation experiments

- Fit joint models with and without GxE across a grid of `sd_ge` values
- Evaluate:
  - recovery of average genotype effects,
  - recovery of genotype instability,
  - change in RMSE and correlation with truth,
  - residual spatial diagnostics

#### Step 4: Add report documentation

- Add a new report subsection describing:
  - when GxE is necessary,
  - how the model was extended,
  - how much variance is attributable to GxE,
  - which genotypes are stable vs unstable

#### Step 5: Optional `sommer` implementation

- After `mgcv` GxE is validated, extend `run_sommer_joint()` with an analogous
  genotype-by-bench random term
- Compare whether `sommer` and `mgcv` give consistent conclusions about GxE

### 8.7 Caveats

- With one observation per genotype per bench, the GxE term must remain
  regularised; an overly rich interaction structure will compete with the
  spatial surface and residual term for the same information.
- The first implementation should therefore keep GxE random, not fixed.
- If GxE is large, the interpretation of a single across-bench genotype BLUE
  becomes less biologically meaningful; bench-specific adjusted estimates should
  then be reported alongside the average effect.
- Real-data validation is essential: simulation alone is not enough to establish
  that the added term is identifying biological interaction rather than
  absorbing residual structure.

---

## 9. Outlier Handling Plan

### 9.1 Motivation

The current framework includes optional outlier replacement via a global
1.5 × IQR rule. This is easy to implement, but it is not the most defensible
approach for spatially heterogeneous multi-bench trials.

Outliers in this setting can arise from several qualitatively different sources:

- **Measurement or data-entry errors**: impossible values, unit errors,
  transcription mistakes
- **Local experimental failures**: empty pots, irrigation failures, disease,
  broken plants
- **Legitimate biological extremes**: genuinely high- or low-performing plots
- **Spatially structured extremes**: values that appear extreme only because the
  underlying spatial model is incomplete

These cases should not all be treated identically. In particular, a pooled
distribution-based rule can incorrectly flag valid observations from benches with
stronger spatial effects or different mean/variance structure.

### 9.2 Statistical principles

The outlier strategy should follow these principles:

1. **Do not detect outliers from raw pooled phenotypes alone**
   - pooled thresholds ignore bench effects and spatial structure
   - they are especially inappropriate when benches differ in intensity

2. **Prefer model-aware residual diagnostics over raw-value thresholds**
   - outliers should be defined relative to the fitted spatial model, not only
     relative to the marginal phenotype distribution

3. **Use robust rules for flagging, not automatic deletion by default**
   - plots should usually be flagged first, then optionally excluded

4. **Separate hard errors from statistical anomalies**
   - impossible values can be removed deterministically
   - statistical outliers should be reviewed or handled by a configurable rule

5. **Operate at the correct level**
   - detection should be at least bench-specific
   - ideally bench × trait-specific, after fitting the spatial model

### 9.3 Recommended detection hierarchy

The most practical and statistically defensible workflow is a three-tier system.

#### Tier 1: deterministic validation checks

Apply these before any model fitting:

- missing genotype, row, column, or bench identifiers
- duplicated plot coordinates within bench
- impossible phenotype values where domain bounds are known
  - e.g. negative yields, BNI outside 0-100 if that scale is fixed
- non-numeric or malformed coordinates

These are not “statistical outliers”; they are data-quality failures and should
be reported separately.

#### Tier 2: bench-specific robust raw screening

Before model fitting, optionally flag extreme values within each bench × trait
using a robust rule such as:

- median absolute deviation (MAD), preferred
- or IQR, if a simpler fallback is needed

Recommended rule:

$$ z^\ast_i = \frac{x_i - \text{median}(x)}{1.4826 \cdot \text{MAD}(x)} $$

and flag if:

$$ |z^\ast_i| > 3.5 $$

Why this is preferred:

- more robust than SD-based z-scores
- less distorted by a few extreme values
- easy to explain and implement
- bench-specific, so it respects between-bench heterogeneity

This stage should produce **flags**, not immediate removal, unless the user
explicitly requests pre-fit exclusion.

#### Tier 3: model-based residual outlier detection

After fitting the spatial model, compute residual diagnostics within bench:

- raw residuals
- standardized residuals, if available
- studentized residuals, if feasible

Recommended practical rule for the first implementation:

- fit the preferred spatial model
- compute residuals within each bench × trait
- flag observations with
  - `|residual| > 3 × residual_SD_within_bench`, or preferably
  - `|robust_residual_z| > 3.5` using MAD on residuals

This is the most statistically meaningful detection layer because it asks:

"Is this plot extreme after accounting for genotype, bench, and spatial trend?"

That is the right question for this project.

### 9.4 Recommended default policy

The recommended default behaviour is:

- **always** run deterministic validation checks
- **optionally** run bench-specific raw screening
- **prefer** model-based residual flagging as the main statistical outlier method
- **do not automatically replace flagged values with `NA` by default**

Default user-facing policy:

- `outlier_method = "none"` by default for exclusion
- `flag_outliers = TRUE` by default for reporting

This keeps the workflow conservative and reproducible. Users can then choose one
of:

- report only
- exclude pre-fit raw outliers
- refit after excluding model-based residual outliers

### 9.5 Methods to support

Recommended methods to implement, in order:

1. **`"none"`**
   - no statistical outlier detection
   - still run validation checks

2. **`"mad_raw"`**
   - bench × trait MAD rule on raw phenotypes
   - good pragmatic baseline

3. **`"iqr_raw"`**
   - bench × trait IQR rule on raw phenotypes
   - keep for backward compatibility, but not preferred

4. **`"mad_resid"`**
   - MAD rule on model residuals within bench × trait
   - recommended default statistical method once implemented

5. **`"sd_resid"`**
   - residual SD rule within bench × trait
   - less robust than MAD, but easy to communicate

The long-run preferred default should be:

- `outlier_method = "mad_resid"`

but only after the implementation supports safe two-stage fitting and clear
reporting.

### 9.6 Implementation sketch

Add an explicit outlier subsystem rather than embedding the logic inside one
small helper.

Recommended helper functions:

```r
validate_trial_data()
flag_outliers_raw()
flag_outliers_residual()
apply_outlier_policy()
summarise_outliers()
```

Suggested behaviour:

- `validate_trial_data()`
  - detect impossible values, duplicate coordinates, malformed inputs

- `flag_outliers_raw(data, method = "mad_raw", by = c("bench", "trait"))`
  - return the original data plus columns such as:
    - `outlier_raw_flag`
    - `outlier_raw_score`
    - `outlier_raw_method`

- `flag_outliers_residual(data, residual_col, method = "mad_resid", by = c("bench", "trait"))`
  - return:
    - `outlier_resid_flag`
    - `outlier_resid_score`
    - `outlier_resid_method`

- `apply_outlier_policy(data, policy = "report")`
  - `"report"`: keep all rows, only annotate
  - `"set_na"`: set flagged phenotype values to `NA`
  - `"drop"`: remove flagged rows before refit

- `summarise_outliers()`
  - generate per-trait and per-bench counts for the output report

### 9.7 Integration into the current workflow

#### `run_spatial_correction()`

Replace the current boolean `outlier_iqr = TRUE/FALSE` with a more explicit API,
for example:

```r
outlier_method = "none"
outlier_policy = "report"
outlier_threshold = 3.5
outlier_scope = "bench"
```

Recommended staged workflow:

1. validate input data
2. run optional raw outlier flagging
3. fit model on full data unless policy says otherwise
4. run residual outlier flagging
5. if requested, refit after exclusion
6. write both the final outputs and an outlier report

#### `run_spatial_gam()`

Add the same interface, but keep the first version simpler:

- support validation checks
- support `mad_raw` and `mad_resid`
- default to reporting only

This keeps the user-facing standalone script aligned with the broader framework.

### 9.8 Output requirements

Replace the current plain-text outlier report with a richer output set:

- `outlier_report.csv`
  - one row per flagged plot
  - includes bench, genotype, coordinates, trait, observed value, score, method,
    and action taken

- `outlier_summary.csv`
  - counts by trait, bench, and method

- optional appendix page in diagnostics PDF
  - heatmap of flagged plots
  - residual distribution before/after exclusion

This is preferable to a free-text report because it is auditable and easier to
reuse downstream.

### 9.9 Recommended first implementation

The most sensible first implementation is:

1. keep deterministic validation checks always on
2. deprecate global pooled IQR filtering
3. add bench-specific `mad_raw` flagging
4. add bench-specific `mad_resid` flagging after model fit
5. default to `policy = "report"` rather than automatic removal

This is a good balance between:

- statistical defensibility
- implementation complexity
- transparency for users

### 9.10 Caveats

- In unreplicated designs, an extreme residual may reflect a true genotype
  effect, a local environmental failure, or model misspecification; it should
  not automatically be treated as error.
- Residual-based outlier detection is only as good as the fitted model; if the
  spatial model is poor, the residual outlier flags may be misleading.
- Any exclusion rule should be logged explicitly so that downstream analyses are
  reproducible and auditable.

---

## 10. Output Specification

### `spatial_diagnostics.pdf`

One page per trait × method. Each page contains a 6-panel grid:
1. Raw observed trait heatmap
2. Corrected (fitted) trait heatmap
3. Estimated spatial trend heatmap
4. Residuals heatmap
5. Variogram of raw values
6. Variogram of residuals

If `method = "all"`: additional comparison page per trait:
- Side-by-side spatial trend heatmaps (SpATS | mgcv | sommer)
- Pairwise BLUE scatter plots with correlation coefficient

### `BLUEs.csv` and/or `BLUPs.csv`

- One row per genotype (× bench when multi-bench)
- **All original file columns preserved** (Genotype, Bench, all auxiliary columns)
- Appended: `<trait>_BLUE`, `<trait>_BLUE_SE` per trait
- If `method = "all"`: separate columns per method (`<trait>_BLUE_SpATS`, `<trait>_BLUE_mgcv`, `<trait>_BLUE_sommer`)
- Same structure for BLUPs

### `fitted_values.csv`

- One row per plot
- **All original file columns preserved** (including row, col, rep, all auxiliary variables)
- Appended per trait: `<trait>_Observed`, `<trait>_Fitted`, `<trait>_Residual`, `<trait>_SpatialTrend`
- Merge uses original data frame as left-hand side — no plots dropped

### `outlier_report.txt`

Per-trait: count and values of outliers replaced by NA.

---

## 11. Function Architecture (`fit_spatial_models.R`)

```
read_input(fn, rda_object)             # read CSV or .rda; return data frame

replace_outliers_with_na(x)            # IQR outlier replacement

adaptive_nseg(n_unique)                # min(max(5, floor(n/2)), 20)

empirical_semivariogram(data, ...)     # manual pairwise semivariance, 15 bins

# Per-bench functions (replicated data or BLUPs for unreplicated)
run_SpATS_bench(bench_data, ...)       # SpATS; adaptive nseg; BLUEs/BLUPs
                                       # Note: spatial surface centred; mean added back to BLUEs/BLUPs
run_mgcv_bench(bench_data, ...)        # gam() te(); configurable smoother_type; BLUEs/BLUPs
run_sommer_bench(bench_data, ...)      # mmes() + spl2Dc(); BLUEs/BLUPs

# Joint multi-bench functions (unreplicated data, BLUEs)
run_mgcv_joint(dfr, ...)              # joint GAM: 0+Genotype+bench_f+te(...,by=bench_f)+nested RE
run_sommer_joint(dfr, ...)            # joint mmes() with column-offset spl2Dc

plot_heatmap(data, row_col, col_col, value_col, ...)   # ggplot2 viridis tile
plot_variogram(vgram, title)                            # empirical semivariogram
plot_diagnostics(bench_data, result, pheno, method, .) # 6-panel patchwork
plot_comparison(bench_data, results_list, pheno, .)    # multi-method comparison

run_spatial_correction(               # main orchestration function (per-bench only)
  fn, trait_cols, bench_col,          # joint functions called directly — not yet integrated
  gt_col, row_col, col_col,
  method, output_type,
  outlier_iqr, output_dir, rda_object
)
```

---

## 12. Usage Examples

```r
source("scripts/fit_spatial_models.R")

# Test 1: Real data (wheatdata), all methods, BLUEs
run_spatial_correction(
  fn          = "data/wheatdata.rda",
  rda_object  = "wheatdata",
  trait_cols  = "yield",
  gt_col      = "geno",
  row_col     = "row",
  col_col     = "col",
  bench_col   = NULL,        # single bench
  method      = "all",
  output_type = "BLUEs",
  outlier_iqr = TRUE,
  output_dir  = "output/wheatdata"
)

# Test 2: Simulated multi-bench data, SpATS only, BLUPs
run_spatial_correction(
  fn          = "data/BNI_simulation.csv",
  trait_cols  = "BNI",
  bench_col   = "Bench",
  gt_col      = "Genotype",
  row_col     = "Row",
  col_col     = "Col",
  method      = "SpATS",
  output_type = "BLUPs",
  outlier_iqr = FALSE,
  output_dir  = "output/simulation"
)
```
