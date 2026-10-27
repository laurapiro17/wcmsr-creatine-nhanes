# Dietary creatine and depression risk in NHANES 2013–March 2020

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/abstract-in%20press-brightgreen)](#status)
[![R](https://img.shields.io/badge/R-survey-blue)](#methods)

Secondary analysis of NHANES extending Bakian et al. 2020 to a post-2013 US cohort. The inverse association between dietary creatine intake and depression replicates — but only once seafood is counted as a source of creatine, a choice every study in this area makes and none states. A spline dose-response model and a formal interaction test with antidepressant medication use are reported alongside it.

## Erratum (October 2026)

The exposure in the published abstract was built on `PF_MPE` from the `heiscore.data` package. Despite its name, that variable does not contain poultry: checked against the USDA Food Patterns Equivalents Database, it equals meat + cured meat + organ meat + eggs. Poultry is the largest protein-food component in this population (about 1.47 oz-eq/day), so roughly a third of the intended meat + poultry + eggs exposure was missing, in all four definitions. `scripts/03_creatine.R` and `scripts/06_sensitivity.R` now rebuild meat + poultry + eggs as `PFALLPROTLEG − PFSEAPLANTLEG`. The code as it stood for the abstract is tagged `abstract-as-published`.

Re-running the identical analysis (same covariates, survey design and N = 14,664; PHQ-9 ≥ 10):

| Definition | Q4 vs Q1, published | Q4 vs Q1, corrected | Per 0.1 g/day, published | Per 0.1 g/day, corrected |
|---|---|---|---|---|
| A (no seafood) | 0.81 (0.60–1.07) | 0.79 (0.59–1.06) | 0.97 (0.94–1.01), p = 0.11 | 0.98 (0.95–1.01), p = 0.20 |
| B (closest to Bakian) | 0.72 (0.54–0.95) | **0.75 (0.56–1.01)** | 0.96 (0.93–1.00), p = 0.043 | **0.97 (0.94–1.00), p = 0.099** |
| C | 0.65 (0.46–0.91) | 0.69 (0.50–0.94) | 0.96 (0.93–0.99), p = 0.014 | 0.97 (0.94–1.00), p = 0.031 |
| D | 0.67 (0.50–0.91) | 0.70 (0.52–0.94) | 0.96 (0.92–0.99), p = 0.027 | **0.97 (0.94–1.00), p = 0.065** |

Three published estimates are no longer statistically significant (in bold). The 20–39 subgroup estimates move from 0.60–0.72 to 0.64–0.80, non-significant throughout. The direction of the finding, a stronger association when seafood is counted, is unchanged, but the abstract's statement that the association replicates under the definition closest to Bakian's no longer holds. A correction is being requested from the journal.

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
Q4 vs Q1 moves from AOR 0.79 (95% CI 0.59–1.06) with meat, poultry and eggs alone to
0.75 (0.56–1.01) once seafood enters, against Bakian's 0.68 (corrected values; see Erratum). The exposure definition
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

- Creatine intake is derived from FPED animal-protein categories (meat + poultry + eggs, rebuilt as `PFALLPROTLEG − PFSEAPLANTLEG` because `PF_MPE` omits poultry, as the primary proxy; sensitivity definitions add seafood from `PF_SSNS` at 30% and 50% shares) and converted to grams of creatine using Bakian's average of 0.11 g per oz-equivalent, i.e. 3.88 mg per g of animal protein. Details in [`refs/creatine_methodology.md`](refs/creatine_methodology.md).
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
