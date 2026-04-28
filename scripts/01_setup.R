# ============================================================
# WCMSR Creatine-Depression NHANES — Setup
# Day 1: install packages + download first cycle as smoke test
# ============================================================

# --- 1. Install packages (one-time, ~2-3 min on Posit Cloud) ---
pkgs <- c("nhanesA", "survey", "dplyr", "tidyr", "ggplot2",
          "haven", "broom", "splines", "mgcv", "rms")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# --- 2. Load ---
library(nhanesA)
library(survey)
library(dplyr)
library(tidyr)
library(haven)

# --- 3. Smoke test: download NHANES 2017-March 2020 pre-pandemic files ---
# Pre-pandemic uses prefix P_
# DEMO = demographics, DPQ = PHQ-9, DR1TOT = day 1 dietary, DR2TOT = day 2 dietary

cat("Downloading NHANES 2017-March 2020 pre-pandemic core files...\n")

demo_p   <- nhanes("P_DEMO")                      # demographics + survey weights
dpq_p    <- nhanes("P_DPQ", translated = FALSE)   # PHQ-9 raw codes 0-3
dr1tot_p <- nhanes("P_DR1TOT")                    # day 1 dietary totals
dr2tot_p <- nhanes("P_DR2TOT")                    # day 2 dietary totals

cat("Sample sizes:\n")
cat("  DEMO   :", nrow(demo_p), "\n")
cat("  DPQ    :", nrow(dpq_p), "\n")
cat("  DR1TOT :", nrow(dr1tot_p), "\n")
cat("  DR2TOT :", nrow(dr2tot_p), "\n")

# --- 4. Quick PHQ-9 sanity check (depression cases adults ≥20) ---
adults <- demo_p %>%
  filter(RIDAGEYR >= 20) %>%
  inner_join(dpq_p, by = "SEQN")

# PHQ-9 = sum of DPQ010..DPQ090, score ≥10 = depression
phq_items <- paste0("DPQ0", c("10","20","30","40","50","60","70","80","90"))
adults$phq9 <- rowSums(adults[, phq_items], na.rm = FALSE)
adults$depressed <- adults$phq9 >= 10

cat("\nAdults ≥20 with valid PHQ-9:", sum(!is.na(adults$phq9)), "\n")
cat("Depression cases (PHQ-9 ≥ 10):", sum(adults$depressed, na.rm = TRUE), "\n")
cat("Crude prevalence:",
    round(mean(adults$depressed, na.rm = TRUE) * 100, 2), "%\n")

# Expected: ~8.6% per Vahratian 2020 — sanity check.

# --- 5. Save raw cache ---
saveRDS(list(demo = demo_p, dpq = dpq_p,
             dr1 = dr1tot_p, dr2 = dr2tot_p),
        "data/raw/nhanes_p_2017_mar2020.rds")

cat("\nDone. If prevalence ~8-9%, methodology is on track.\n")
