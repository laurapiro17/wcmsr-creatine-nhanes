# ============================================================
# WCMSR Creatine-Depression NHANES — Day 5b (sensitivity)
# 4 definitions of dietary animal protein → creatine
#   A. PF_MPE only (current main analysis)
#   B. PF_MPE + 0.30 × PF_SSNS (US-avg ~30% of SSNS is seafood)
#   C. PF_MPE + (PFSEAPLANTLEG - PF_LEGUMES) (animal+plant)
#   D. PF_MPE + 0.50 × PF_SSNS (high-seafood-share scenario)
# Add education to covariates (Bakian had it; we missed it).
# Run for each: adjusted Q4 vs Q1, continuous, 20-39 stratum.
# ============================================================

library(heiscore.data)
library(dplyr)
library(tidyr)
library(survey)
options(survey.lonely.psu = "adjust")

KCREAT_PER_OZ <- 0.11

data("fped_1314"); data("fped_1516"); data("fped_1720")

# ---- helper to compute mean across day 1/day 2 for any column set ----
mean_oz <- function(df, cols) {
  m <- as.matrix(df[, cols, drop = FALSE])
  rowMeans(matrix(as.numeric(m), nrow = nrow(df)), na.rm = TRUE)
}

build_animal <- function(df, defn) {
  d1 <- intersect(c("DR1_PF_MPE","DR1_PF_SSNS","DR1_PFSEAPLANTLEG","DR1T_PF_LEGUMES"), names(df))
  d2 <- intersect(c("DR2_PF_MPE","DR2_PF_SSNS","DR2_PFSEAPLANTLEG","DR2T_PF_LEGUMES"), names(df))
  for (c in c(d1,d2)) df[[c]] <- as.numeric(df[[c]])

  mpe   <- (df$DR1_PF_MPE + df$DR2_PF_MPE) / 2
  ssns  <- (df$DR1_PF_SSNS + df$DR2_PF_SSNS) / 2
  seapl <- (df$DR1_PFSEAPLANTLEG + df$DR2_PFSEAPLANTLEG) / 2
  legs  <- (df$DR1T_PF_LEGUMES + df$DR2T_PF_LEGUMES) / 2

  switch(defn,
         A = mpe,
         B = mpe + 0.30 * ssns,
         C = mpe + (seapl - legs),
         D = mpe + 0.50 * ssns)
}

run_for_definition <- function(defn) {
  cat(sprintf("\n========== Definition %s ==========\n", defn))

  cre <- bind_rows(
    data.frame(SEQN = fped_1314$SEQN, animal_oz = build_animal(fped_1314, defn)),
    data.frame(SEQN = fped_1516$SEQN, animal_oz = build_animal(fped_1516, defn)),
    data.frame(SEQN = fped_1720$SEQN, animal_oz = build_animal(fped_1720, defn))
  )
  cre$creatine_g <- cre$animal_oz * KCREAT_PER_OZ
  cre$SEQN <- as.integer(cre$SEQN)

  # Load the Day 4 prepared dataset (has weights, covariates already recoded)
  m4 <- readRDS("data/processed/day4_models.rds")
  a <- m4$analytic
  a$SEQN <- as.integer(a$SEQN)

  # Replace creatine with new definition
  a <- a %>% select(-creatine_g, -animal_oz, -creatine_q) %>%
       left_join(cre %>% select(SEQN, animal_oz, creatine_g), by = "SEQN")

  a <- a %>% filter(!is.na(creatine_g), !is.na(phq9), wt_combined > 0)

  q <- quantile(a$creatine_g, probs = seq(0,1,0.25), na.rm = TRUE)
  a$creatine_q <- cut(a$creatine_g, breaks = q,
                     include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))

  cat(sprintf("  N=%d | median=%.3f g/d | Q4 cutoff=%.3f\n",
              nrow(a), median(a$creatine_g, na.rm=TRUE),
              q[4]))

  # Add education if available
  if ("DMDEDUC2" %in% names(a)) {
    a$educ <- factor(as.character(a$DMDEDUC2))
  } else {
    a$educ <- factor("Unknown")
  }

  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                   weights = ~wt_combined, data = a, nest = TRUE)

  cov_form <- "+ age + sex + race + income + bmi + smoke + pa_active + healthcare + antidep + cycle"
  if (length(unique(a$educ)) > 1) cov_form <- paste(cov_form, "+ educ")

  # Q4 vs Q1
  fit_q <- tryCatch(svyglm(as.formula(paste("depressed ~ creatine_q", cov_form)),
                           design = des, family = quasibinomial()), error=function(e) NULL)
  if (!is.null(fit_q)) {
    co <- coef(fit_q); ci <- confint(fit_q)
    q4 <- grep("creatine_qQ4", names(co))
    OR <- exp(co[q4]); CI <- exp(ci[q4,])
    cat(sprintf("  Q4 vs Q1 adj : AOR = %.3f (%.3f - %.3f)\n", OR, CI[1], CI[2]))
  }

  # Continuous (per 0.1 g/d)
  fit_c <- tryCatch(svyglm(as.formula(paste("depressed ~ I(creatine_g*10)", cov_form)),
                           design = des, family = quasibinomial()), error=function(e) NULL)
  if (!is.null(fit_c)) {
    co <- summary(fit_c)$coefficients
    rn <- "I(creatine_g * 10)"
    OR <- exp(co[rn,"Estimate"])
    CI <- exp(co[rn,"Estimate"] + c(-1.96,1.96)*co[rn,"Std. Error"])
    p  <- co[rn,"Pr(>|t|)"]
    cat(sprintf("  Continuous   : AOR/0.1g = %.3f (%.3f - %.3f) p=%.3f\n",
                OR, CI[1], CI[2], p))
  }

  # 20-39 stratum
  a$age_grp <- ifelse(a$age <= 39, "20-39", "40+")
  des2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                    weights = ~wt_combined, data = a, nest = TRUE)
  sub <- subset(des2, age_grp == "20-39")
  fit_y <- tryCatch(svyglm(as.formula(paste("depressed ~ creatine_q", cov_form)),
                           design = sub, family = quasibinomial()), error=function(e) NULL)
  if (!is.null(fit_y)) {
    co <- coef(fit_y); ci <- confint(fit_y)
    q4 <- grep("creatine_qQ4", names(co))
    if (length(q4)) {
      OR <- exp(co[q4]); CI <- exp(ci[q4,])
      cat(sprintf("  20-39 stratum: AOR = %.3f (%.3f - %.3f) | n=%d\n",
                  OR, CI[1], CI[2], nrow(sub$variables)))
    }
  }
}

for (defn in c("A","B","C","D")) run_for_definition(defn)

cat("\n========== Sensitivity summary ==========\n")
cat("If A-D give similar AOR, finding is robust to seafood handling.\n")
cat("If C (animal+plant) gives much weaker effect, plant proteins dilute signal.\n")
