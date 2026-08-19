# BiometricalCLPM

**Biometrical Cross-Lagged Panel Models for longitudinal twin data.**

The Biometrical CLPM integrates the *biometrical (genetic) simplex* and *Direction-of-Causation*
(DoC) models into the twin cross-lagged panel model, so that stability is modelled as separate
genetic and environmental developmental processes, cross-lagged causation can be tested against
cross-lagged genetic/environmental confounding, and proximal (cross-sectional) causal effects can
be estimated alongside distal (lagged) ones.

This repository contains reusable OpenMx builders for the model (`R/`), a worked general
specification (`Model/`), and the analysis scripts for the application to smoking quantity and
depressive symptoms in TEDS young adults (`TEDS/`).

## Citation

> Singh, M., Hunter, M. D., Assary, E., Verhulst, B., Peterson, R. E., Maes, H. H. M.,
> Dolan, C. V., Eley, T. C., & Neale, M. C. (2025). *Causation Between Smoking Quantity and
> Depressive Symptoms in Young Adults: Evidence from Novel Cross-Lagged Twin Models*.
> medRxiv. https://doi.org/10.1101/2025.11.18.25340516

Preregistration: https://doi.org/10.17605/OSF.IO/K47BF

<details>
<summary><strong>Abstract</strong></summary>

**Background**: Smoking and depression often co-occur in young adults, potentially due to
reciprocal causal effects and/or shared underlying etiological factors. Here, we tested causal
hypotheses between smoking quantity (cigarettes per day; *CigDay*) and depressive symptoms
(*DepSx*) using novel Biometrical Cross-Lagged Panel Models (CLPMs), which integrate twin
developmental and direction-of-causation analyses for more robust causal inference in longitudinal
twin studies.

**Methods**: Study sample included 10,034 participants (61.6% females; 4,112 twin pairs) from the
Twins Early Development Study, with up to six repeated assessments spanning ages 21–29. The
Biometrical CLPM provided three key innovations over standard twin CLPM: distinct autoregressive
processes of latent genetic and environmental factors; cross-lagged genetic and environmental
liabilities, contrasting confounding and causation; and cross-sectional causal effects,
distinguishing between proximal (short-term) and distal (lagged) causation.

**Results**: In the two largest waves assessed approximately four years apart, genetic liabilities
for both variables remained highly stable, while environmental influences showed considerable
temporal variation. In the standard CLPM, *CigDay* predicted *DepSx* four years later. However, in
the Biometrical CLPM, the cross-lagged phenotypic associations were attenuated and non-significant
when accounting for the differential genetic and environmental developmental processes. Adding the
cross-sectional causal paths revealed significant bidirectional *proximal* effects, with a stronger
*CigDay→DepSx* effect. Across all six waves with intervals up to two years, analyses consistently
showed significant bidirectional cross-lagged associations.

**Conclusions**: The findings support a reinforcing feedback loop between *CigDay* and *DepSx*
during young adulthood, involving bidirectional effects that persist for up to two years but
dissipate over longer intervals.

</details>

---

## Requirements

| | |
|---|---|
| R | ≥ 4.1 (analyses in the paper: 4.4.0) |
| [OpenMx](https://openmx.ssri.psu.edu/) | ≥ 2.21.0 (analyses in the paper: 2.21.3) |

```r
install.packages("OpenMx")
```

`R/biometricalCLPM.R` needs **only OpenMx**. The TEDS data-preparation scripts additionally use
`dplyr`, `tidyr`, `haven`, `ggplot2`, `GGally` and `psych`; the univariate ACE scripts use
`miFunctions.R` (Hermine Maes' helper functions, bundled here).

> The univariate ACE scripts set the NPSOL optimizer via `miFunctions.R`. NPSOL is only available
> in the OpenMx build from the [University of Virginia repository](https://openmx.ssri.psu.edu/installing-openmx);
> with the CRAN build, use `mxOption(NULL, "Default optimizer", "CSOLNP")` (or SLSQP) instead.

---

## Quick start

```r
source("R/biometricalCLPM.R")

## Column names the model expects, for two phenotypes over 3 waves
vn <- clpmVarNames("X", "Y", waveLabels = 1:3)
vn$selVars
#> "X1_Tw1" "Y1_Tw1" "X2_Tw1" "Y2_Tw1" "X3_Tw1" "Y3_Tw1"
#> "X1_Tw2" "Y1_Tw2" "X2_Tw2" "Y2_Tw2" "X3_Tw2" "Y3_Tw2"

## 1. Standard twin CLPM (Fig. 1A): phenotypic autoregression + cross-lags
m0   <- biometricalCLPM("X", "Y", waveLabels = 1:3,
                        mzData = mzDat, dzData = dzDat,      # wide, one row per family
                        estimator = "ML", components = c("A", "E"))
fit0 <- mxRun(m0)

## 2. Biometrical CLPM (Fig. 1B): stability as separate A and E simplex processes
fit1 <- mxRun(useBiometricalAR(fit0))
mxCompare(fit1, fit0)          # is phenotypic autoregression too restrictive?

## 3. Sources of the cross-lagged association: causation vs. confounding (Fig. 1C)
fitHybA <- mxRun(useLiabilityCrossLags(fit1, "A"))        # phenotypic CL + genetic CL
fitHybE <- mxRun(useLiabilityCrossLags(fit1, "E"))        # phenotypic CL + environmental CL
fitConf <- mxRun(useLiabilityCrossLags(fit1, c("A", "E"), dropPhenotypic = TRUE))
fitNull <- mxRun(setPhenotypicCrossLags(fit1, free = FALSE))   # no cross-lagged association

mxCompare(fitHybA, fit1)                     # causation is nested in the hybrid model
mxCompare(fit1, fitNull)                     # any cross-lagged association at all?
mxCompare(fitConf, list(fitHybA, fitHybE))   # equal df — select on AIC

## 4. Proximal direction-of-causation at the last wave (Fig. 1D)
fitA <- mxRun(addProximalDoC(fitNull, wave = 3, direction = "x->y", dropCov = "A"))
fitB <- mxRun(addProximalDoC(fitNull, wave = 3, direction = "y->x", dropCov = "A"))
fitC <- mxRun(addProximalDoC(fitNull, wave = 3, direction = "both", dropCov = c("A", "E")))
mxCompare(fitNull, list(fitA, fitB, fitC))   # equal df — select on AIC

## Standardised estimates
clpmStdEstimates(fit1, "stdBetaA", se = TRUE)   # genetic autoregression / cross-lags
clpmStdEstimates(fit1, "stdBeta")               # phenotypic paths
```

Each modifier returns a **new, unfitted `MxModel`**; pass it to `mxRun()` (or `mxTryHard()`) yourself.
Modifiers accept fitted models as input, so the previous solution supplies the start values — which
is how the published scripts walk through the model space.

### The model space

| Paper | Model | Call |
|---|---|---|
| Fig. 1A | Standard twin CLPM (phenotypic AR) | `biometricalCLPM(...)` |
| Fig. 1B | Biometrical CLPM (A/C/E simplex AR) | `useBiometricalAR(fit)` |
| — | Hybrid: phenotypic cross-lags *and* liability cross-lags | `useLiabilityCrossLags(fit, "A")` |
| Fig. 1C | Cross-lagged genetic/environmental confounding, no phenotypic causation | `useLiabilityCrossLags(fit, c("A","E"), dropPhenotypic = TRUE)` |
| — | No cross-lagged association ("null") | `setPhenotypicCrossLags(fit, free = FALSE)` |
| Fig. 1D | Proximal (cross-sectional) DoC | `addProximalDoC(fit, wave, direction, dropCov)` |

The causation model (Fig. 1B) is nested within each hybrid model, and the null model is nested
within the causation model. The hybrid models and the pure confounding model have the *same* number
of free parameters, so they are compared on pseudo-AIC rather than by a difference test — likewise
the competing proximal DoC models. Rejecting confounding while causation fits as well is what
strengthens the causal interpretation of the cross-lagged paths.

Identification of the DoC submodels rests on the two phenotypes having **different ACE variance
decompositions**. At a wave with A, C and E variance there are five potential sources of within-wave
covariance (two causal paths, `covA`, `covC`, `covE`) and at most three may be estimated; with AE
variance, at most two of four. `addProximalDoC()` fixes the covariances you name in `dropCov` and
warns if you free more causal paths than you drop covariances. Confirm with
`mxCheckIdentification()`.

---

## Function reference (`R/biometricalCLPM.R`)

| Function | Purpose |
|---|---|
| `biometricalCLPM()` | Build the multi-group (MZ/DZ) twin CLPM. Continuous and/or ordinal phenotypes; ML or WLS; A/C/E or A/E. |
| `useBiometricalAR()` | Fix phenotypic autoregression at 0, free the latent A/C/E autoregressions. |
| `useLiabilityCrossLags()` | Free cross-lagged paths between latent A/C/E factors; optionally drop the phenotypic cross-lags. |
| `setPhenotypicCrossLags()` | Free or fix the phenotypic (distal) cross-lags, by direction and/or wave. |
| `addProximalDoC()` | Free proximal cross-sectional causal path(s) at a wave, fixing the corresponding within-wave innovation covariance(s). |
| `clpmVarNames()` | Manifest variable names / `selVars` in model order. |
| `clpmLabels()` | Every parameter label, grouped by matrix and role. |
| `clpmStart()` | Start values. |
| `clpmStdEstimates()` | Tidy a standardised algebra (`stdBeta`, `stdBetaA`, `psiA`, `rA`, …) into a data frame, optionally with `mxSE()` standard errors. |
| `clpmSumStats()` | Build WLS summary statistics for one zygosity group and strip the raw data. |

### Data your model needs

* **Wide (family-wise) format**, one row per twin pair, columns named as `clpmVarNames(...)$selVars`.
* **Separate MZ and DZ data frames.** Opposite-sex DZ pairs were pooled with same-sex DZ pairs in
  the paper; supply your own grouping if you want them modelled separately.
* **Ordinal phenotypes** must be `mxFactor()`s (`clpmSumStats()` does this for you) and are handled
  through the liability-threshold model.
* Incomplete pairs are retained under FIML (`estimator = "ML"`); under WLS they contribute to the
  pairwise summary statistics.

### Estimation

The paper used **weighted least squares** on summary statistics, because full-information ML with
several ordinal variables per person becomes intractable at six waves:

```r
sv  <- clpmVarNames("DepSxResRN_T", "CigDay3L_T", waveLabels = 1:6)$selVars
ov  <- grep("CigDay3L_T", sv, value = TRUE)
mzs <- clpmSumStats(mzWide, sv, ordVars = ov, nThresholds = 2)   # raw data stripped
dzs <- clpmSumStats(dzWide, sv, ordVars = ov, nThresholds = 2)

m <- biometricalCLPM("DepSxResRN_T", "CigDay3L_T", waveLabels = 1:6,
                     mzData = mzs, dzData = dzs,
                     estimator = "WLS", components = c("A", "E"),
                     ordinal = "y", nThresholds = 2)
```

Under WLS, goodness-of-fit is Browne's pseudo-χ² with a corresponding pseudo-AIC, and nested models
are compared with the Satorra–Bentler scaled-difference χ² (`mxCompare()` reports these). When two
models have the same number of parameters the SB χ² is necessarily zero; the paper then selected on
pseudo-AIC, treating ΔAIC > 2 as a meaningfully worse fit.

### Parameter labels

`s` and `t` are consecutive wave labels; `w` is a single wave label. Labels are stable across
analyses, which is why the two-wave TEDS models use `waveLabels = c(1, 6)` — the same paths carry
the same names as in the six-wave models.

| Label | Matrix | Meaning |
|---|---|---|
| `AR_bx{t}x{s}`, `AR_by{t}y{s}` | `BE` | phenotypic autoregression |
| `CL_by{t}x{s}` | `BE` | **distal** (lagged) effect x → y |
| `CL_bx{t}y{s}` | `BE` | **distal** (lagged) effect y → x |
| `Inst_by{w}x{w}` | `BE` | **proximal** (cross-sectional) effect x → y |
| `Inst_bx{w}y{w}` | `BE` | **proximal** (cross-sectional) effect y → x |
| `VAx{w}`, `VAy{w}`, `covAxy{w}` | `PSA` | innovation (co)variance of the A factors (likewise `VC*`/`VE*`) |
| `AR_Ax{t}x{s}`, `AR_Ay{t}y{s}` | `BEA` | genetic autoregression (likewise `AR_C*`, `AR_E*`) |
| `CL_Ay{t}x{s}`, `CL_Ax{t}y{s}` | `BEA` | cross-lagged genetic liability (likewise `CL_C*`, `CL_E*`) |
| `b0_{var}` | `b0` | intercept / mean |
| `th{k}_{var}` | `threG` | threshold *k* of an ordinal variable |

`clpmLabels(waveLabels)` returns all of them programmatically:

```r
L <- clpmLabels(c(1, 6))
L$pheno$clyx     # "CL_by6x1"   distal x -> y
L$A$ar           # "AR_Ax6x1" "AR_Ay6y1"
L$pheno$instxy   # "Inst_bx1y1" "Inst_bx6y6"
```

### Algebras available on a fitted model

`V` (expected phenotypic covariance), `VA`/`VC`/`VE`, `rA`/`rC`/`rE` (genetic and environmental
correlations), `psiA`/`psiC`/`psiE` (standardised innovation covariances), `stdBeta`,
`stdBetaA`/`stdBetaC`/`stdBetaE` (standardised regressions), and `US` (unstandardised and
standardised variance components). All are at the top level of the fitted object:

```r
round(fit1$stdBetaA$result, 2)
round(mxSE(stdBetaA, fit1), 3)
```

---

## Repository layout

```
R/
  biometricalCLPM.R        Reusable model builders (OpenMx only). Start here.

Model/
  Biometrical_CLPM_Model_MultipleTimePoints.R
                           Original flat specification of the model for two
                           continuous phenotypes over an arbitrary number of
                           waves, with the full nested model sequence written
                           out. Kept as the readable reference for what the
                           functions in R/ build.

TEDS/                      Analyses reported in the manuscript
  1_DataWrangle_Residualise_DepSx_CigDay.R
        Cleaning and exclusions; sMFQ sum scores; CigDay 3-level ordinal
        variable; regression of DepSx on age and sex followed by rank-based
        inverse-normal transformation.
        in : Data/teds_v2_Jun2024.sav
        out: Data/TEDS_DepSx_byFam_Jun2024.tsv, Data/TEDS_CigDay_byFam_Jun2024.tsv

  2a_UnivariateACE_DepSx.R / 2b_UnivariateACE_CigDay.R
        Univariate ACE models for each wave (continuous, and ordinal with
        age/sex as definition variables).

  3_Get_DepSx_CigDay_SumStats_for_WLS.R
        Merges the two phenotypes, splits by zygosity and builds the WLS
        summary statistics used by scripts 4 and 5.
        out: Out/{MZ,DZ}sumStatDepSmk.RData        (six waves)
             Out/{MZ,DZ}sumStatDepSmk_T16.RData    (TEDS21 + TEDS26)

  4_BiometricalCLPM_Dep_Smk_TEDS21_26.R
        Primary analysis. Two waves (TEDS21 = T1, TEDS26 = T6, ~4 years apart):
        standard twin CLPM, Biometrical CLPM, cross-lagged confounding models,
        and the proximal DoC models.

  5_BiometricalCLPM_Dep_Smk_TEDS21_COVIDall_26.R
        All six waves (TEDS21, COVID Phases 1–4, TEDS26), intervals up to
        two years.

miFunctions.R              Label/value helpers used by the univariate ACE
                           scripts (Hermine Maes).
```

The TEDS scripts read from `Data/` and write to `Out/`, `Plots/` and `Logs/` relative to the
working directory; create those directories before running.

## Data availability

TEDS data are available on request; see the
[TEDS data access policy](https://www.teds.ac.uk/researchers/teds-data-access-policy).
No individual-level data are included in this repository. The WLS summary statistics produced by
script 3 (and by `clpmSumStats()`, which strips the raw data from the returned object) contain no
individual-level information.

## Reproducing the published analyses

With TEDS data in `Data/`, run `TEDS/1` → `TEDS/3` → `TEDS/4` (two-wave, the primary analysis) and
`TEDS/5` (six-wave). Scripts `2a`/`2b` produce the univariate ACE estimates reported in the
supplement and are not required by `3`–`5`. The scripts are written to be stepped through
interactively; they print model summaries and write parameter tables to `Out/`.

`R/biometricalCLPM.R` builds matrices that are element-for-element identical (values, labels and
free/fixed status) to those constructed in `Model/` and `TEDS/4`–`5`, so it reproduces the published
models rather than approximating them.

## Applying the model to a new twin cohort

1. Get your data into wide family-wise format with columns named by `clpmVarNames()`.
2. Fit univariate ACE models per wave first, to check that A, C and E are estimable and that the two
   phenotypes differ in their variance decomposition — the DoC part of the model has no leverage if
   they do not.
3. Build with `biometricalCLPM()` and walk the model space in the order given in
   [The model space](#the-model-space), comparing each step with `mxCompare()`.
4. Check `mxCheckIdentification()` after adding proximal DoC paths.
5. If a model fails to converge, `mxTryHard()` and better start values (`clpmStart()`) usually
   resolve it; with WLS, confirm your summary statistics are positive definite
   (`eigen(cov2cor(mzs$observedStats$cov))$values`).

## Note on `Model/` vs `R/`

`Model/Biometrical_CLPM_Model_MultipleTimePoints.R` is the original specification and is kept
deliberately flat and explicit: every matrix is constructed and printed in place, which makes it the
best document for reading *what the model is*. It is a template to copy and edit, not a script to
run as-is — the data objects (`datmz_`, `datdz_`, `selVars_`) are placeholders you fill in.
`R/biometricalCLPM.R` is the same model wrapped in functions, and is what you want if you are
applying it to your own data.

## Licence

MIT — see [LICENSE](LICENSE).
