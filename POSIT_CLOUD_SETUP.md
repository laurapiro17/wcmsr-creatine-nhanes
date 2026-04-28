# Posit Cloud Setup — 5 minuts

## Pas 1 — Crear compte (1 min)
1. Vés a **https://posit.cloud/**
2. Sign Up amb el teu compte de GitHub (laurapiro17)
3. Tria **Cloud Free** (25h/mes, suficient per 12 dies)

## Pas 2 — Crear projecte (1 min)
1. Botó **New Project** → **New RStudio Project**
2. Nom: `wcmsr-creatine-nhanes`
3. Espera ~30s a que arrenqui

## Pas 3 — Pujar fitxers (1 min)
A la pestanya **Files** (panel inferior dret):
1. Click **Upload**
2. Puja: `~/Projects/wcmsr-creatine-nhanes/scripts/01_setup.R`
3. Crea carpeta `data/raw/` (botó **New Folder** dos cops)

## Pas 4 — Executar setup (~3 min)
1. Obre `01_setup.R` al panel superior esquerre
2. Click **Source** (botó dalt a la dreta del panel) o `Cmd+Shift+S`
3. Espera la instal·lació de paquets + download

## Què hauries de veure

```
Sample sizes:
  DEMO   : ~15500
  DPQ    : ~8965
  DR1TOT : ~14300
  DR2TOT : ~12000

Adults ≥20 with valid PHQ-9: ~5000
Depression cases (PHQ-9 ≥ 10): ~430
Crude prevalence: ~8.6%
```

Si surt **~8-9% prevalence** → metodologia OK. Si no, em dius el número i debuguem.

## Si tens problemes
- `nhanesA` falla: prova `install.packages("nhanesA", repos="https://cloud.r-project.org")`
- Memory issue Posit Cloud free: tanca altres projectes; aquest no hauria de passar de 1GB RAM
