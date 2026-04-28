# ============================================================
# WCMSR Creatine-Depression NHANES — Day 3
# Compute dietary creatine intake from FPED via heiscore.data
# Apply Bakian: 0.11 g creatine per oz_eq of animal protein
# ============================================================

# --- 1. Install + load heiscore.data ---
if (!"heiscore.data" %in% rownames(installed.packages())) {
  install.packages("heiscore.data")
}
library(heiscore.data)
library(dplyr)
library(tidyr)

# --- 2. Inspect column structure of one cycle to confirm names ---
cat("=== FPED 2017-March 2020 column names ===\n")
data("fped_1720")
print(head(names(fped_1720), 60))
cat("\nDimensions:", dim(fped_1720), "\n")
cat("First row example:\n")
print(head(fped_1720, 2))

# --- 3. Define animal protein columns (Bakian's MyPyramid subgroups) ---
# Standard FPED column names. We grep them defensively in case of variant naming.
animal_pf_patterns <- c("PF_MEAT", "PF_CUREDMEAT", "PF_ORGAN",
                        "PF_POULT", "PF_SEAFD_HI", "PF_SEAFD_LOW")

find_animal_cols <- function(df) {
  hits <- character()
  for (p in animal_pf_patterns) {
    m <- grep(paste0("^", p, "$"), names(df), value = TRUE, ignore.case = TRUE)
    hits <- c(hits, m)
  }
  unique(hits)
}

# --- 4. Compute creatine per-person per-day for each cycle ---
# Bakian: total animal protein (oz_eq) × 0.11 g creatine / oz_eq
KCREAT_PER_OZ <- 0.11   # grams of creatine per oz equivalent of animal protein

compute_creatine <- function(fped, cycle_label, day_label) {
  ac <- find_animal_cols(fped)
  if (length(ac) == 0) {
    cat("WARNING: no animal protein columns matched in", cycle_label, day_label, "\n")
    return(NULL)
  }
  cat("Cycle", cycle_label, day_label, "- columns matched:", paste(ac, collapse=", "), "\n")

  # SEQN column (some FPED tables use SEQN, some seqn)
  seqn_col <- intersect(c("SEQN", "seqn"), names(fped))[1]
  if (is.na(seqn_col)) stop("No SEQN column in ", cycle_label)

  # Day indicator column might be DRDDAY (1 or 2), or files might be split per day
  if ("DRDDAY" %in% names(fped)) {
    fped <- fped %>% rename(day = DRDDAY)
  } else if ("drdday" %in% names(fped)) {
    fped <- fped %>% rename(day = drdday)
  } else {
    # If no day column, assume single day
    fped$day <- as.integer(day_label)
  }

  fped %>%
    rename(SEQN = !!seqn_col) %>%
    mutate(across(all_of(ac), ~ as.numeric(.))) %>%
    mutate(animal_oz = rowSums(across(all_of(ac)), na.rm = TRUE),
           creatine_g = animal_oz * KCREAT_PER_OZ) %>%
    select(SEQN, day, animal_oz, creatine_g)
}

# Load all 3 cycles
data("fped_1314"); data("fped_1516"); data("fped_1720")

cat("\n=== Computing creatine per cycle ===\n")
cre_1314 <- compute_creatine(fped_1314, "2013-14", "both")
cre_1516 <- compute_creatine(fped_1516, "2015-16", "both")
cre_1720 <- compute_creatine(fped_1720, "2017-Mar2020", "both")

# --- 5. Average day 1 and day 2 creatine per person ---
avg_creatine <- function(df) {
  if (is.null(df)) return(NULL)
  df %>%
    group_by(SEQN) %>%
    summarise(
      n_days       = n(),
      animal_oz    = mean(animal_oz, na.rm = TRUE),
      creatine_g   = mean(creatine_g, na.rm = TRUE),
      .groups = "drop"
    )
}

cre_1314_avg <- avg_creatine(cre_1314)
cre_1516_avg <- avg_creatine(cre_1516)
cre_1720_avg <- avg_creatine(cre_1720)

cat("\n=== Per-person creatine summary (mean across days) ===\n")
for (lab_df in list(list("2013-14", cre_1314_avg),
                    list("2015-16", cre_1516_avg),
                    list("2017-Mar2020", cre_1720_avg))) {
  if (!is.null(lab_df[[2]])) {
    s <- summary(lab_df[[2]]$creatine_g)
    cat(sprintf("  %s : N=%d  median=%.3f g/d  mean=%.3f g/d  Q3=%.3f g/d\n",
                lab_df[[1]], nrow(lab_df[[2]]),
                s["Median"], s["Mean"], s["3rd Qu."]))
  }
}

# --- 6. Merge with adults_all_cycles from Day 2 ---
adults <- readRDS("data/processed/adults_all_cycles.rds")

cre_all <- bind_rows(
  cre_1314_avg %>% mutate(cycle = "2013-2014"),
  cre_1516_avg %>% mutate(cycle = "2015-2016"),
  cre_1720_avg %>% mutate(cycle = "2017-Mar2020")
)

adults_cre <- adults %>% left_join(cre_all, by = c("SEQN","cycle"))

cat("\n=== Merge result ===\n")
cat("Adults total:", nrow(adults_cre), "\n")
cat("With creatine value:", sum(!is.na(adults_cre$creatine_g)), "\n")
cat("With creatine + valid PHQ-9:",
    sum(!is.na(adults_cre$creatine_g) & !is.na(adults_cre$phq9)), "\n")

# --- 7. Smoke test: replicate Bakian quartile pattern (univariate) ---
analytic <- adults_cre %>%
  filter(!is.na(creatine_g), !is.na(phq9))
analytic$creatine_q <- cut(analytic$creatine_g,
  breaks = quantile(analytic$creatine_g, probs = seq(0,1,0.25), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))

cat("\n=== Quartile cutoffs (g/day creatine) ===\n")
print(quantile(analytic$creatine_g, probs = seq(0,1,0.25), na.rm = TRUE))

cat("\n=== Crude depression prevalence by creatine quartile ===\n")
print(analytic %>% group_by(creatine_q) %>%
        summarise(N = n(),
                  cases = sum(depressed, na.rm = TRUE),
                  prev_pct = round(mean(depressed, na.rm = TRUE)*100, 2)))

# Univariate OR Q4 vs Q1 (no covariates yet, no survey weights)
m0 <- glm(depressed ~ creatine_q, data = analytic, family = binomial)
cat("\n=== Univariate OR vs Q1 (NO weights, NO covariates) ===\n")
print(round(exp(cbind(OR = coef(m0), confint.default(m0))), 3))

cat("\nBakian Q4-vs-Q1 AOR ref: 0.68 (CI 0.52-0.88) — full adjustment + weights\n")
cat("If our crude OR is in the 0.5-0.85 range, we're on the right track.\n")

# --- 8. Save analytic dataset ---
saveRDS(adults_cre, "data/processed/adults_with_creatine.rds")
cat("\nSaved data/processed/adults_with_creatine.rds\n")
cat("\nDay 3 done. Next (Day 4-5): survey weights + adjusted logistic + splines.\n")
