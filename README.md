# WCMSR 2026 — Creatine & Depression in NHANES

**Conference**: 5th IJMS World Conference of Medical Student Research (virtual, 11-12 Jul 2026)
**Abstract deadline**: 10 May 2026
**Fee**: $40 USD (PayPal @editorijms)

## Question
Does dietary creatine intake remain inversely associated with depression risk in U.S. adults using post-Bakian NHANES cycles (2013–March 2020), and what is the dose-response shape and interaction with antidepressant medication?

## Precedent
- Bakian et al. 2020 *Translational Psychiatry* — NHANES 2005-2012, N=22,692, AOR 0.68 (Q4 vs Q1).
- Ostojic et al. 2025 *Nutritional Neuroscience* — KNHANES 2022 Korean cohort, quartiles only, no medication interaction.

## Novelty (3 angles)
1. Post-Bakian US cohort (2013–March 2020), ~15-20K participants combining 3 cycles.
2. Restricted cubic splines for dose-response (vs Bakian's quartiles).
3. Formal interaction test with antidepressant/anxiolytic use (Ostojic flagged this as a limitation).

## Cycles
- 2013-2014 (suffix `_H`)
- 2015-2016 (suffix `_I`)
- 2017-March 2020 pre-pandemic (prefix `P_`)

## Files needed per cycle
| File | Purpose |
|---|---|
| `DEMO` | age, sex, race, income, weights |
| `DR1TOT` + `DR2TOT` | dietary recall day 1 & 2 |
| `DPQ` | PHQ-9 depression screener |
| `RXQ_RX` | prescription medications (antidepressant flag) |
| `BMX` | BMI |
| `SMQ` | smoking |
| `PAQ` | physical activity |
| `HUQ` | healthcare access |

## Creatine calculation (Bakian)
Creatine intake = sum across food items of (food gram weight × creatine content per gram). Creatine content per food code from FNDDS + ad-hoc lookup table for meat/fish.

## Schedule (28 Apr → 10 May)
| Day | Task |
|---|---|
| 1-2 | Setup R env + download all NHANES files |
| 3 | Combine cycles + construct survey weights `WTSAF` |
| 4-5 | Replicate Bakian quartile analysis as smoke test |
| 6-7 | Restricted cubic splines dose-response |
| 8-9 | Interaction test medication × creatine |
| 10-11 | Write abstract (250-300 words) + sensitivity |
| 12 | AI-audit + cross-check + submit ($40 PayPal) |
