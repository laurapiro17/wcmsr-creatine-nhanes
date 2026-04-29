# ============================================================
# WCMSR Creatine-Depression NHANES — Day 4
# Survey-weighted adjusted logistic regression
# Covariates following Bakian 2020:
#   age, sex, race/ethnicity, income, education, BMI,
#   smoking, physical activity, healthcare access,
#   antidepressant/anxiolytic use
# Combined weight per NCHS guidance: WT × cycle_years / total_years
# ============================================================

library(dplyr)
library(survey)
library(tidyr)

a <- readRDS("data/processed/adults_with_creatine.rds")

# --- 1. Diagnostic: what columns are available? ---
cat("=== Column inventory in adults_with_creatine ===\n")
nm <- names(a)

inv <- list(
  weights_2day = grep("WTDR2D",   nm, value = TRUE),
  weights_pp   = grep("WTDR2DPP", nm, value = TRUE),
  psu          = grep("SDMVPSU",  nm, value = TRUE),
  stratum      = grep("SDMVSTRA", nm, value = TRUE),
  age          = grep("RIDAGEYR", nm, value = TRUE),
  sex          = grep("RIAGENDR", nm, value = TRUE),
  race         = grep("RIDRETH",  nm, value = TRUE),
  income       = grep("INDFMPIR", nm, value = TRUE),
  educ         = grep("DMDEDUC",  nm, value = TRUE),
  bmi          = grep("BMXBMI",   nm, value = TRUE),
  smoke        = grep("SMQ020",   nm, value = TRUE),
  paq          = grep("PAQ",      nm, value = TRUE)[1:3],
  huq          = grep("HUQ030",   nm, value = TRUE)
)
for (k in names(inv)) cat(sprintf("  %-12s : %s\n", k, paste(inv[[k]], collapse=", ")))

# --- 2. Fetch dietary survey weights (not in DEMO; live in DR1TOT) ---
library(nhanesA)
cat("\nFetching dietary weights from DR1TOT files...\n")

get_weights <- function(table, wt_col) {
  d <- nhanes(table)
  if (!"SEQN" %in% names(d)) return(NULL)
  if (!wt_col %in% names(d)) {
    cat("  WARNING:", wt_col, "not found in", table, "\n")
    return(NULL)
  }
  d[, c("SEQN", wt_col)] %>% setNames(c("SEQN", "wt_orig"))
}

w_H <- get_weights("DR1TOT_H", "WTDR2D")
w_I <- get_weights("DR1TOT_I", "WTDR2D")
# Pre-pandemic uses WTDRD1PP and WTDR2DPP — try both names
w_P <- tryCatch(get_weights("P_DR1TOT", "WTDR2DPP"),
                error = function(e) get_weights("P_DR1TOT", "WTDRD1PP"))

# Fall back: scan column names if not found
if (is.null(w_P)) {
  d <- nhanes("P_DR1TOT")
  cat("P_DR1TOT weight columns available:",
      paste(grep("^WT", names(d), value = TRUE), collapse=", "), "\n")
  # Use whichever exists (prefer 2-day, fallback 1-day)
  for (cand in c("WTDR2DPP","WTDRD1PP","WTDR2D","WTDRD1")) {
    if (cand %in% names(d)) {
      w_P <- d[, c("SEQN", cand)] %>% setNames(c("SEQN","wt_orig"))
      cat("Using", cand, "as P-cycle weight.\n"); break
    }
  }
}

w_all <- bind_rows(
  w_H %>% mutate(cycle_letter = "H"),
  w_I %>% mutate(cycle_letter = "I"),
  w_P %>% mutate(cycle_letter = "P")
)
w_all$SEQN <- as.integer(w_all$SEQN)
a$SEQN <- as.integer(a$SEQN)

a <- a %>% left_join(w_all %>% select(SEQN, wt_orig), by = "SEQN")

# --- 2b. Construct combined weight ---
# Cycle years: 2013-14 = 2y, 2015-16 = 2y, 2017-March 2020 = 3.2y
total_years <- 2 + 2 + 3.2
a$cycle_years <- c(H = 2, I = 2, P = 3.2)[a$cycle]
a$wt_combined <- a$wt_orig * a$cycle_years / total_years

cat("\n=== Combined weight summary ===\n")
print(summary(a$wt_combined))
cat("Sum of weights (should ≈ US adult pop):", round(sum(a$wt_combined, na.rm=TRUE)/1e6, 1), "M\n")

# --- 3. Prepare analytic dataset (complete-case for adjusted model) ---
analytic <- a %>%
  filter(!is.na(creatine_g), !is.na(phq9), !is.na(wt_combined),
         !is.na(SDMVPSU), !is.na(SDMVSTRA))

# Quartiles using SURVEY-WEIGHTED quantiles (standard practice)
# But here we just use unweighted quartiles for now to be comparable to crude
q <- quantile(analytic$creatine_g, probs = seq(0,1,0.25), na.rm = TRUE)
analytic$creatine_q <- cut(analytic$creatine_g, breaks = q,
                           include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))

cat("\n=== Analytic sample size (complete-case, with weights+PSU+strata) ===\n")
cat(nrow(analytic), "\n")

# --- 4. Recode covariates ---
analytic <- analytic %>%
  mutate(
    age       = as.numeric(RIDAGEYR),
    sex       = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    race      = factor(RIDRETH1),
    income    = as.numeric(INDFMPIR),
    bmi       = as.numeric(BMXBMI),
    smoke     = ifelse(SMQ020 == 1, 1, 0),       # 1 = smoked ≥100 cig
    pa_active = ifelse(!is.na(PAQ605) & PAQ605 == 1, 1, 0),
    healthcare= ifelse(!is.na(HUQ030) & HUQ030 == 1, 1, 0),
    antidep   = as.integer(has_antidep_or_anxio | has_antidep | has_anxio),
    cycle     = factor(cycle)
  )

# --- 5. Survey design ---
des <- svydesign(
  ids     = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~wt_combined,
  data    = analytic,
  nest    = TRUE
)
options(survey.lonely.psu = "adjust")

# --- 6. Crude weighted Q4 vs Q1 ---
m_crude_w <- svyglm(depressed ~ creatine_q,
                    design = des, family = quasibinomial())
cat("\n=== Crude WEIGHTED OR (no covariates) ===\n")
print(round(exp(cbind(OR = coef(m_crude_w), confint(m_crude_w))), 3))

# --- 7. Adjusted weighted logistic ---
m_adj <- svyglm(depressed ~ creatine_q + age + sex + race + income +
                  bmi + smoke + pa_active + healthcare + antidep + cycle,
                design = des, family = quasibinomial())

cat("\n=== ADJUSTED WEIGHTED OR (Bakian-style covariates) ===\n")
adj_tab <- round(exp(cbind(OR = coef(m_adj), confint(m_adj))), 3)
print(adj_tab[grep("creatine_q", rownames(adj_tab)), ])

cat("\nBakian Q4-vs-Q1 AOR: 0.68 (CI 0.52-0.88)\n")

# --- 8. Save model objects for Day 5 (splines) ---
saveRDS(list(analytic = analytic, design = des,
             m_crude_w = m_crude_w, m_adj = m_adj,
             adj_tab = adj_tab),
        "data/processed/day4_models.rds")

cat("\nDay 4 done. Saved data/processed/day4_models.rds\n")
cat("Next (Day 5): restricted cubic splines + interaction antidep × creatine.\n")
