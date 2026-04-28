# ============================================================
# WCMSR Creatine-Depression NHANES — Day 3 (v2)
# heiscore.data only provides HEI-2020 aggregates (PF_MPE).
# Use PF_MPE = Meat + Poultry + Eggs as animal protein proxy.
# Bakian: 0.11 g creatine per oz_eq of meat/poultry.
# Limitation: excludes seafood (~15% of US creatine intake).
# ============================================================

if (!"heiscore.data" %in% rownames(installed.packages())) {
  install.packages("heiscore.data")
}
library(heiscore.data)
library(dplyr)
library(tidyr)

KCREAT_PER_OZ <- 0.11   # Bakian average

# --- 1. Inspect 3 cycles for column-name variation ---
data("fped_1314"); data("fped_1516"); data("fped_1720")

inspect <- function(df, name) {
  mpe   <- grep("PF_MPE",  names(df), value = TRUE, ignore.case = TRUE)
  ssns  <- grep("PF_SSNS", names(df), value = TRUE, ignore.case = TRUE)
  legs  <- grep("PF_LEGUMES", names(df), value = TRUE, ignore.case = TRUE)
  seapl <- grep("PFSEAPLANTLEG", names(df), value = TRUE, ignore.case = TRUE)
  cat(sprintf("%s: MPE=%s | SSNS=%s | LEG=%s | SEAPL=%s\n",
              name,
              paste(mpe, collapse=","),
              paste(ssns, collapse=","),
              paste(legs, collapse=","),
              paste(seapl, collapse=",")))
}
cat("=== Column inventory ===\n")
inspect(fped_1314, "2013-14")
inspect(fped_1516, "2015-16")
inspect(fped_1720, "2017-Mar2020")

# --- 2. Compute creatine per person (mean of day 1 + day 2 PF_MPE) ---
compute_creatine <- function(df, label) {
  mpe_cols <- grep("PF_MPE", names(df), value = TRUE, ignore.case = TRUE)
  if (length(mpe_cols) < 1) stop("No PF_MPE columns in ", label)

  # Find SEQN col
  seqn_col <- intersect(c("SEQN","seqn"), names(df))[1]

  # Subset and coerce to numeric
  df2 <- df %>% select(all_of(c(seqn_col, mpe_cols))) %>%
    rename(SEQN = !!seqn_col) %>%
    mutate(across(all_of(mpe_cols), ~ as.numeric(.)))

  # Mean across day 1 and day 2 (NA-aware)
  df2$animal_oz <- rowMeans(df2[, mpe_cols, drop = FALSE], na.rm = TRUE)
  df2$animal_oz[is.nan(df2$animal_oz)] <- NA
  df2$creatine_g <- df2$animal_oz * KCREAT_PER_OZ

  df2 %>% select(SEQN, animal_oz, creatine_g)
}

cre_1314 <- compute_creatine(fped_1314, "2013-14")
cre_1516 <- compute_creatine(fped_1516, "2015-16")
cre_1720 <- compute_creatine(fped_1720, "2017-Mar2020")

cat("\n=== Per-person animal protein + creatine summary ===\n")
for (lab_df in list(list("2013-14", cre_1314),
                    list("2015-16", cre_1516),
                    list("2017-Mar2020", cre_1720))) {
  s <- summary(lab_df[[2]]$creatine_g)
  cat(sprintf("  %s : N=%d  median=%.3f g/d  mean=%.3f g/d  Q3=%.3f g/d  max=%.3f g/d\n",
              lab_df[[1]], nrow(lab_df[[2]]),
              s["Median"], s["Mean"], s["3rd Qu."], s["Max."]))
}

# --- 3. Merge with adults_all_cycles from Day 2 ---
adults <- readRDS("data/processed/adults_all_cycles.rds")
cre_all <- bind_rows(
  cre_1314 %>% mutate(cycle = "2013-2014"),
  cre_1516 %>% mutate(cycle = "2015-2016"),
  cre_1720 %>% mutate(cycle = "2017-Mar2020")
)

adults_cre <- adults %>% left_join(cre_all, by = c("SEQN","cycle"))

cat("\n=== Merge result ===\n")
cat("Adults total:", nrow(adults_cre), "\n")
cat("With creatine value:", sum(!is.na(adults_cre$creatine_g)), "\n")
cat("With creatine + valid PHQ-9:",
    sum(!is.na(adults_cre$creatine_g) & !is.na(adults_cre$phq9)), "\n")

# --- 4. Smoke test: replicate Bakian quartile pattern (univariate) ---
analytic <- adults_cre %>% filter(!is.na(creatine_g), !is.na(phq9))

q <- quantile(analytic$creatine_g, probs = seq(0,1,0.25), na.rm = TRUE)
analytic$creatine_q <- cut(analytic$creatine_g, breaks = q,
                           include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))

cat("\n=== Quartile cutoffs (g/day creatine) ===\n")
print(round(q, 3))
cat("Bakian reference: Q1 0-0.26, Q4 0.70-3.16\n")

cat("\n=== Crude depression prevalence by creatine quartile ===\n")
print(analytic %>% group_by(creatine_q) %>%
        summarise(N = n(),
                  cases = sum(depressed, na.rm = TRUE),
                  prev_pct = round(mean(depressed, na.rm = TRUE)*100, 2),
                  .groups="drop"))

m0 <- glm(depressed ~ creatine_q, data = analytic, family = binomial)
cat("\n=== Univariate OR vs Q1 (NO weights, NO covariates) ===\n")
print(round(exp(cbind(OR = coef(m0), confint.default(m0))), 3))

cat("\nBakian Q4-vs-Q1 AOR ref: 0.68 (CI 0.52-0.88) — adjusted + weighted\n")
cat("Crude OR Q4-vs-Q1 expected to be similar or stronger (no adjustment yet).\n")

# --- 5. Save analytic dataset ---
saveRDS(adults_cre, "data/processed/adults_with_creatine.rds")
cat("\nSaved data/processed/adults_with_creatine.rds\n")
cat("\nDay 3 done. Next (Day 4-5): survey weights + adjusted logistic + splines.\n")
cat("Limitation logged: PF_MPE proxy excludes seafood (~15% of US creatine).\n")
