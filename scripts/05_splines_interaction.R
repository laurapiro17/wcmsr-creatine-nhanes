# ============================================================
# WCMSR Creatine-Depression NHANES — Day 5
# 1) Continuous + restricted cubic splines (dose-response)
# 2) Effect modification by antidepressant/anxiolytic use
# 3) Sex stratification (Bakian found stronger effect in females)
# ============================================================

library(dplyr)
library(survey)
library(splines)

m <- readRDS("data/processed/day4_models.rds")
analytic <- m$analytic
des <- m$design

# Drop zero-weight again
analytic2 <- analytic[analytic$wt_combined > 0 & !is.na(analytic$wt_combined), ]
des2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                  weights = ~wt_combined, data = analytic2, nest = TRUE)
options(survey.lonely.psu = "adjust")

cat("Analytic N (with weights):", nrow(analytic2), "\n\n")

# ============================================================
# 1. CONTINUOUS CREATINE (per 0.1 g/day increase)
# ============================================================
cat("=== 1. Continuous creatine (per 0.1 g/day) ===\n")
m_cont <- svyglm(depressed ~ I(creatine_g*10) + age + sex + race + income +
                   bmi + smoke + pa_active + healthcare + antidep + cycle,
                 design = des2, family = quasibinomial())
co <- summary(m_cont)$coefficients
or_cont <- exp(co["I(creatine_g * 10)", "Estimate"])
ci_cont <- exp(co["I(creatine_g * 10)", "Estimate"] +
                 c(-1.96, 1.96) * co["I(creatine_g * 10)", "Std. Error"])
cat(sprintf("  AOR per 0.1 g/day creatine = %.3f (%.3f - %.3f) | p=%.4f\n",
            or_cont, ci_cont[1], ci_cont[2],
            co["I(creatine_g * 10)", "Pr(>|t|)"]))

# ============================================================
# 2. RESTRICTED CUBIC SPLINES (3 knots) — pre-compute basis columns
# ============================================================
cat("\n=== 2. Restricted cubic splines (3 knots) ===\n")
knots <- quantile(analytic2$creatine_g,
                  probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
cat("Knots at:", round(knots, 3), "\n")

# Pre-compute spline basis as columns in analytic2, then re-create design
spline_mat <- splines::ns(analytic2$creatine_g,
                          knots = knots[2],
                          Boundary.knots = knots[c(1, 3)])
colnames(spline_mat) <- c("cre_ns1","cre_ns2")
analytic2 <- cbind(analytic2, as.data.frame(spline_mat))

des2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                  weights = ~wt_combined, data = analytic2, nest = TRUE)

m_rcs <- svyglm(depressed ~ cre_ns1 + cre_ns2 +
                  age + sex + race + income + bmi + smoke + pa_active +
                  healthcare + antidep + cycle,
                design = des2, family = quasibinomial())

# Test for non-linearity: compare linear vs spline (Wald on cre_ns2 != 0)
nl_test <- tryCatch(regTermTest(m_rcs, ~ cre_ns2),
                    error = function(e) NULL)
if (!is.null(nl_test)) {
  cat("Non-linearity test (cre_ns2 = 0):\n"); print(nl_test)
}

# Print spline coefficients
cat("Spline term coefficients (link scale):\n")
co <- summary(m_rcs)$coefficients
print(round(co[c("cre_ns1","cre_ns2"), ], 4))

# ============================================================
# 3. INTERACTION: creatine × antidepressant/anxiolytic use
# ============================================================
cat("\n=== 3. Interaction creatine × antidepressant ===\n")
m_int <- svyglm(depressed ~ creatine_q * antidep + age + sex + race + income +
                  bmi + smoke + pa_active + healthcare + cycle,
                design = des2, family = quasibinomial())

# Wald test for the interaction term(s)
int_idx <- grep(":antidep", names(coef(m_int)))
int_test <- tryCatch(regTermTest(m_int, ~ creatine_q:antidep),
                     error = function(e) NULL)
if (!is.null(int_test)) {
  cat("Interaction Wald test:\n"); print(int_test)
}

# Stratified estimates
cat("\n--- Stratified Q4-vs-Q1 AOR by antidepressant use ---\n")
for (ad in c(0, 1)) {
  sub <- subset(des2, antidep == ad)
  fit <- tryCatch(svyglm(depressed ~ creatine_q + age + sex + race + income +
                           bmi + smoke + pa_active + healthcare + cycle,
                         design = sub, family = quasibinomial()),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    co <- coef(fit); ci <- confint(fit)
    q4 <- grep("creatine_qQ4", names(co))
    if (length(q4)) {
      OR <- exp(co[q4]); CI <- exp(ci[q4, ])
      label <- ifelse(ad == 0, "No antidep/anxio", "On antidep/anxio")
      cat(sprintf("  %-18s : AOR = %.3f (%.3f - %.3f) | n=%d\n",
                  label, OR, CI[1], CI[2], nrow(sub$variables)))
    }
  }
}

# ============================================================
# 4. SEX STRATIFICATION (Bakian: stronger in females)
# ============================================================
cat("\n=== 4. Sex stratification ===\n")
for (s in levels(analytic2$sex)) {
  sub <- subset(des2, sex == s)
  fit <- tryCatch(svyglm(depressed ~ creatine_q + age + race + income +
                           bmi + smoke + pa_active + healthcare + antidep + cycle,
                         design = sub, family = quasibinomial()),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    co <- coef(fit); ci <- confint(fit)
    q4 <- grep("creatine_qQ4", names(co))
    if (length(q4)) {
      OR <- exp(co[q4]); CI <- exp(ci[q4, ])
      cat(sprintf("  %-8s : AOR = %.3f (%.3f - %.3f) | n=%d\n",
                  s, OR, CI[1], CI[2], nrow(sub$variables)))
    }
  }
}

# ============================================================
# 5. AGE STRATIFICATION (Bakian: stronger in 20-39)
# ============================================================
cat("\n=== 5. Age stratification (20-39 vs 40+) ===\n")
analytic2$age_grp <- factor(ifelse(analytic2$age <= 39, "20-39", "40+"))
des3 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                  weights = ~wt_combined, data = analytic2, nest = TRUE)
for (g in c("20-39","40+")) {
  sub <- subset(des3, age_grp == g)
  fit <- tryCatch(svyglm(depressed ~ creatine_q + sex + race + income +
                           bmi + smoke + pa_active + healthcare + antidep + cycle,
                         design = sub, family = quasibinomial()),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    co <- coef(fit); ci <- confint(fit)
    q4 <- grep("creatine_qQ4", names(co))
    if (length(q4)) {
      OR <- exp(co[q4]); CI <- exp(ci[q4, ])
      cat(sprintf("  Age %-6s : AOR = %.3f (%.3f - %.3f) | n=%d\n",
                  g, OR, CI[1], CI[2], nrow(sub$variables)))
    }
  }
}

cat("\n=== DAY 5 SUMMARY ===\n")
cat(sprintf("Continuous: AOR per 0.1 g/d = %.3f (%.3f-%.3f)\n",
            or_cont, ci_cont[1], ci_cont[2]))
cat("See above for splines, interaction, sex+age stratification.\n")
cat("Bakian comparators: females (0.62), 20-39 (0.52), no-meds (stronger)\n")
