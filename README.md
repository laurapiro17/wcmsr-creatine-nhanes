# Dietary creatine and depression risk in NHANES 2013–March 2020

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/abstract-in%20press-brightgreen)](#status)
[![R](https://img.shields.io/badge/R-survey-blue)](#methods)

Secondary analysis of NHANES extending Bakian et al. 2020 to a post-2013 US cohort. The inverse association between dietary creatine intake and depression replicates — but only once seafood is counted as a source of creatine, a choice every study in this area makes and none states. A spline dose-response model and a formal interaction test with antidepressant medication use are reported alongside it.

## Research question

Does dietary creatine intake remain inversely associated with depression risk in US adults across NHANES cycles 2013–March 2020, and how much does the answer depend on how the exposure is defined?

## Background

Bakian et al. 2020 (*Translational Psychiatry*) reported an inverse association between dietary creatine intake and PHQ-9 depression risk in NHANES 2005–2012 (N = 22,692; AOR 0.68, Q4 vs Q1). Ostojic et al. 2025 (*Nutritional Neuroscience*) extended this to a Korean cohort (KNHANES 2022) using quartiles, and flagged the lack of a formal interaction test with antidepressant medication as a key limitation.

## Contribution

This analysis extends the prior work in three ways.

It draws on three NHANES cycles (2013–March 2020, roughly 15,000–20,000 adults),
fully independent of the 2005–2012 sample analysed by Bakian.

It re-runs the same survey-weighted model under four definitions of what counts as
dietary creatine, differing only in whether seafood is included and at what share.
Q4 vs Q1 moves from AOR 0.81 (95% CI 0.60–1.07) with meat, poultry and eggs alone to
0.72 (0.54–0.95) once seafood enters, close to Bakian's 0.68. The exposure definition
is doing much of the work, and it is a choice studies in this area make silently.

And it adds two analyses a quartile contrast cannot give: a spline dose-response
model, so threshold, plateau or U-shaped patterns become visible, and a test of the
multiplicative interaction between dietary creatine and prescription antidepressant
or anxiolytic use, which is the gap Ostojic et al. 2025 flagged.

## Data

NHANES public-release files, three cycles:

| Cycle | File suffix |
|---|---|
| 2013–2014 | `_H` |
| 2015–2016 | `_I` |
| 2017–March 2020 (pre-pandemic) | `P_` |

Per-cycle files: `DEMO` (demographics + survey weights), `DR1TOT` and `DR2TOT` (24-hour dietary recall, days 1 and 2), `DPQ` (PHQ-9), `RXQ_RX` (prescription medications), `BMX` (BMI), `SMQ` (smoking), `PAQ` (physical activity), `HUQ` (healthcare access). Pulled via the `nhanesA` R package.

## Methods

- Creatine intake is derived from FPED animal-protein categories (`PF_MPE`, Meat + Poultry + Eggs, as the primary proxy; sensitivity definitions add seafood from `PF_SSNS` at 30% and 50% shares) and converted to grams of creatine using Bakian's average of 0.11 g per oz-equivalent, i.e. 3.88 mg per g of animal protein. Details in [`refs/creatine_methodology.md`](refs/creatine_methodology.md).
- Depression is defined as PHQ-9 ≥ 10 (moderate-to-severe).
- The primary result is a survey-weighted logistic regression (`survey::svyglm`) on creatine quartiles, re-run under each of the four exposure definitions. The adjustment set follows Bakian 2020: age, sex, race/ethnicity, income-to-poverty ratio, education, BMI, smoking, physical activity, healthcare access.
- The dose-response model uses a restricted cubic spline basis (`splines::ns`, interior knot at the median with boundary knots at the 10th and 90th percentiles of intake), with a design-corrected Wald test for non-linearity.
- The interaction test uses a multiplicative term between creatine and antidepressant or anxiolytic use, assessed with a design-corrected Wald test.
- Sex stratification was pre-specified, motivated by the larger effect estimate Bakian reported in females.
- Combined-cycle weights follow NCHS guidance (`weight × cycle_years / total_years`), with strata `SDMVSTRA` and clusters `SDMVPSU`.
- The four exposure definitions and alternative covariate sets are in [`scripts/06_sensitivity.R`](scripts/06_sensitivity.R).

## Repository structure

```
scripts/
  01_setup.R               Package install + first-cycle smoke test
  02_download_cycles.R     NHANES 2013-14, 2015-16, 2017-Mar 2020 + covariates
  03_creatine.R            FPED → creatine intake derivation
  04_adjusted.R            Survey-weighted logistic regression (Bakian replication)
  05_splines_interaction.R Restricted cubic splines + medication interaction
  06_sensitivity.R         4 creatine definitions, education-adjusted models

refs/
  creatine_methodology.md  Per-food creatine content, Bakian 2020 derivation

POSIT_CLOUD_SETUP.md       Reproducibility on Posit Cloud Free (no local R needed)
CITATION.cff               Machine-readable citation metadata
```

## Status

The abstract was accepted for oral presentation (Original Research) at the 5th IJMS
World Conference of Medical Student Research, held virtually on 11–12 July 2026, and
is in press for the conference supplement of the *International Journal of Medical
Students*.

## Reproducibility

End-to-end re-run from raw NHANES files:

1. Clone the repository.
2. Open `scripts/01_setup.R` in R (≥ 4.3) or Posit Cloud — installs `nhanesA`, `survey`, `rms`, `dplyr`, `tidyr`, `ggplot2`, `splines`, `mgcv`, `broom`, `haven`.
3. Run scripts in numbered order. Expected sample sizes after `01_setup.R` are documented in [`POSIT_CLOUD_SETUP.md`](POSIT_CLOUD_SETUP.md).

## Citation

Use the *Cite this repository* button on the GitHub sidebar (auto-generated from `CITATION.cff`), or:

> Piñero Roig, L. (2026). *Dietary Creatine Intake and Depression in U.S. Adults: a Multi-Cycle NHANES Replication and Sensitivity Analysis of Bakian 2020.* GitHub repository, https://github.com/laurapiro17/wcmsr-creatine-nhanes

## License

MIT — see [LICENSE](LICENSE).
