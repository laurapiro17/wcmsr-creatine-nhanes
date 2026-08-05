# Running the analysis on Posit Cloud

No local R installation is needed. Posit Cloud Free (25 h/month) covers the whole
pipeline; peak memory stays under 1 GB.

## Setup

1. Create a project at https://posit.cloud/ (**New Project** → **New RStudio Project**).
2. Upload `scripts/01_setup.R` from the **Files** pane.
3. Create the folder `data/raw/`.
4. Open `01_setup.R` and click **Source** (`Cmd+Shift+S`). Package installation plus
   the first NHANES download takes about 3 minutes.

Then run the remaining scripts in numbered order.

## Expected output of `01_setup.R`

```
Sample sizes:
  DEMO   : ~15500
  DPQ    : ~8965
  DR1TOT : ~14300
  DR2TOT : ~12000

Adults >=20 with valid PHQ-9: ~5000
Depression cases (PHQ-9 >= 10): ~430
Crude prevalence: ~8.6%
```

A crude prevalence of 8-9% is consistent with Vahratian et al. 2020 and indicates
the PHQ-9 recode is working. A markedly different figure points to a coding error
in the DPQ block rather than to a real finding.

## Troubleshooting

- If `nhanesA` fails to install:
  `install.packages("nhanesA", repos = "https://cloud.r-project.org")`
- NHANES servers occasionally time out. Re-running the download block is safe:
  files are cached under `data/raw/`.
