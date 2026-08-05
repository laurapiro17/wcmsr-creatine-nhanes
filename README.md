# Dietary creatine and depression risk in NHANES 2013–March 2020

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-under%20peer%20review-orange)](#status)
[![R](https://img.shields.io/badge/R-survey%20%C2%B7%20rms-blue)](#methods)

Secondary analysis of NHANES extending Bakian et al. 2020 with a post-2013 US cohort, a dose-response model based on restricted cubic splines, and a formal interaction test with antidepressant medication use.

## Research question

Does dietary creatine intake remain inversely associated with depression risk in US adults across NHANES cycles 2013–March 2020, and how is this association modified by antidepressant medication use?

## Background

Bakian et al. 2020 (*Translational Psychiatry*) reported an inverse association between dietary creatine intake and PHQ-9 depression risk in NHANES 2005–2012 (N = 22,692; AOR 0.68, Q4 vs Q1). Ostojic et al. 2025 (*Nutritional Neuroscience*) extended this to a Korean cohort (KNHANES 2022) using quartiles, and flagged the lack of a formal interaction test with antidepressant medication as a key limitation.

## Contribution

This analysis extends the prior work in three ways.

It draws on three NHANES cycles (2013–March 2020, roughly 15,000–20,000 adults),
fully independent of the 2005–2012 sample analysed by Bakian. It models the
dose-response with restricted cubic splines (`rms::rcs`, 4 knots) rather than
quartile cut-offs, so threshold, plateau or U-shaped patterns that a Q4 vs Q1
contrast cannot show become visible. And it tests the multiplicative interaction
between dietary creatine and prescription antidepressant or anxiolytic use, which
is the gap Ostojic et al. 2025 flagged.

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
- The main model is a survey-weighted logistic regression (`survey::svyglm`) with restricted cubic splines on creatine intake. The adjustment set follows Bakian 2020: age, sex, race/ethnicity, income-to-poverty ratio, education, BMI, smoking, physical activity, healthcare access.
- The interaction test uses a multiplicative term between creatine and antidepressant or anxiolytic use, assessed with a design-corrected Wald test.
- Sex stratification was pre-specified, motivated by the larger effect estimate Bakian reported in females.
- Combined-cycle weights follow NCHS guidance (`weight × cycle_years / total_years`), with strata `SDMVSTRA` and clusters `SDMVPSU`.
- Sensitivity analyses cover four creatine-intake definitions and alternative covariate sets (`scripts/06_sensitivity.R`).

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

Abstract submitted to the **5th IJMS World Conference of Medical Student Research** (virtual, 11–12 July 2026). Currently under peer review.

## Reproducibility

End-to-end re-run from raw NHANES files:

1. Clone the repository.
2. Open `scripts/01_setup.R` in R (≥ 4.3) or Posit Cloud — installs `nhanesA`, `survey`, `rms`, `dplyr`, `tidyr`, `ggplot2`, `splines`, `mgcv`, `broom`, `haven`.
3. Run scripts in numbered order. Expected sample sizes after `01_setup.R` are documented in [`POSIT_CLOUD_SETUP.md`](POSIT_CLOUD_SETUP.md).

## Citation

Use the *Cite this repository* button on the GitHub sidebar (auto-generated from `CITATION.cff`), or:

> Piñero Roig, L. (2026). *Dietary Creatine Intake and Depression Risk in NHANES 2013–March 2020: a Dose-Response Re-analysis.* GitHub repository, https://github.com/laurapiro17/wcmsr-creatine-nhanes

## License

MIT — see [LICENSE](LICENSE).
