# Creatine Calculation Methodology

## Bakian 2020 approach (replicate exactly)

**Average creatine content**: 3.88 g per kg of animal protein (≈ 3.88 mg/g).

**Range**: 0.06 g/oz (frankfurters/sausage/luncheon meats, ≈ 2.12 mg/g) → 0.16 g/oz (fish high in n-3 fatty acids, ≈ 5.65 mg/g). Average across sources = 0.11 g/oz.

**Source**: extensive literature review, key reference Balsom et al. 1994 *Sports Medicine*.

**Conversion**: 1 oz = 28.3495 g, so 0.11 g/oz = 3.88 mg/g.

## Foods INCLUDED (MyPyramid / FPED categories)

- Meat: beef, pork, veal, lamb, game
- Organ meats
- Frankfurters / sausage / luncheon meats
- Poultry: chicken, turkey
- Fish/shellfish high in n-3 fatty acids
- Fish/shellfish low in n-3 fatty acids

## Foods EXCLUDED

- Eggs (low/no creatine)
- Dairy (low/no creatine)
- Cranberries (Bakian explicitly excluded — minor source)
- **Supplemental creatine** (analyzed separately)

## Per-food creatine content (from literature, raw)

| Food | mg/g (raw) | Source |
|---|---|---|
| Pork | 5.0 | Various reviews |
| Beef | 4.5 | Various reviews |
| Salmon | 4.5 | Various reviews |
| Tuna | 4.0 | Various reviews |
| Chicken | 3.4 | Various reviews |
| Turkey | ~3.0 | Various reviews |
| Frankfurters/sausage | 2.12 | Bakian Table S1 lower bound |
| Fish high n-3 (e.g., herring) | 5.65 | Bakian Table S1 upper bound |
| **AVERAGE (Bakian)** | **3.88** | Bakian 2020 |

Cooked values are lower than raw. NHANES dietary recall captures foods as consumed,
i.e. usually cooked, but Bakian applied the 3.88 g/kg average regardless. Replicated
here for comparability.

## Computation per participant

```
creatine_intake_g_per_day = 
  (animal_protein_grams_day1 × 3.88e-3) average with
  (animal_protein_grams_day2 × 3.88e-3)
```

Using **MPED** (NHANES ≤ 2015-2016) or **FPED** (NHANES 2015-2016 onward) summed equivalents for meat/poultry/fish categories.

## Bakian's quartile cutoffs (for reference / sanity check)

- Q1: 0 – 0.26 g/day
- Q4: 0.70 – 3.16 g/day

If our quartile cutoffs differ wildly from these in the post-Bakian cohort, that's a flag for methodology error.

## Validation paper (NHANES 2017-2018 with same 3.88 g/kg)

Smith-Ryan et al. 2021 (PMC8498075) — used the 3.88 g/kg Bakian coefficient in NHANES 2017-2018 for older adults (≥65, N=1,221) for cardiovascular/metabolic outcomes (NOT depression). Confirms methodology transferable to modern cycles.

## Decisions taken

- **Blanket coefficient instead of per-food values.** Bakian's Supplementary Tables
  S1/S2 would allow per-food precision, but Smith-Ryan 2021 applied the blanket
  3.88 g/kg to modern cycles, so the blanket coefficient is used here.
- **Supplement users.** Excluded in a sensitivity analysis, as in Bakian
  (NHANES dietary supplements file, DSQ).
- **FPED vs MPED in 2013-2014.** The 2013-2014 cycle has its own FPED release,
  which is the one used.
