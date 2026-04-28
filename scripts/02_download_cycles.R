# ============================================================
# WCMSR Day 2 — Download cycles _H (2013-14) + _I (2015-16)
# + covariates for all 3 cycles + harmonize into adults table
# ============================================================

library(nhanesA)
library(dplyr)
library(tidyr)

# --- 1. Cycles to download ---
# _H = 2013-2014, _I = 2015-2016, P_ = 2017-March 2020 pre-pandemic
cycles <- list(
  H  = list(suffix = "_H",  label = "2013-2014"),
  I  = list(suffix = "_I",  label = "2015-2016"),
  P  = list(suffix = "P_",  label = "2017-Mar2020")
)

# --- 2. File set per cycle ---
# Core: DEMO, DPQ, DR1TOT, DR2TOT
# Covariates: BMX (BMI), SMQ (smoking), PAQ (physical activity),
#             HUQ (healthcare access), RXQ_RX (prescription meds)
file_bases <- c("DEMO", "DPQ", "DR1TOT", "DR2TOT",
                "BMX", "SMQ", "PAQ", "HUQ", "RXQ_RX")

# --- 3. Helper to build NHANES file name ---
build_name <- function(base, suffix) {
  if (startsWith(suffix, "_")) paste0(base, suffix)
  else paste0(suffix, base)
}

# --- 4. Download everything ---
all_data <- list()
for (cyc_id in names(cycles)) {
  cyc <- cycles[[cyc_id]]
  cat("\n=== Cycle", cyc$label, "(", cyc$suffix, ") ===\n")
  cyc_data <- list()
  for (base in file_bases) {
    fname <- build_name(base, cyc$suffix)
    cat("  ", fname, "...")
    df <- tryCatch(
      nhanes(fname, translated = (base != "DPQ")),
      error = function(e) {
        cat(" FAILED:", conditionMessage(e), "\n"); NULL
      }
    )
    if (!is.null(df)) {
      cat(" n=", nrow(df), "\n")
      cyc_data[[base]] <- df
    }
  }
  all_data[[cyc_id]] <- cyc_data
}

# --- 5. Build harmonized adults table per cycle ---
phq_items <- paste0("DPQ0", c("10","20","30","40","50","60","70","80","90"))

# Common antidepressant generic names (lowercase) for string match against RXDDRUG
antidep_names <- c(
  "fluoxetine","sertraline","citalopram","escitalopram","paroxetine","fluvoxamine",
  "venlafaxine","duloxetine","desvenlafaxine","levomilnacipran",
  "bupropion","mirtazapine","trazodone","vilazodone","vortioxetine","nefazodone",
  "amitriptyline","nortriptyline","imipramine","desipramine","clomipramine","doxepin","trimipramine","protriptyline",
  "phenelzine","tranylcypromine","selegiline","isocarboxazid"
)
anxiolytic_names <- c(
  "alprazolam","lorazepam","diazepam","clonazepam","oxazepam","temazepam",
  "buspirone","hydroxyzine"
)

# Collapse RXQ_RX to one row per SEQN
summarize_rx <- function(rx) {
  if (is.null(rx)) return(NULL)
  rx %>%
    mutate(drug = tolower(as.character(RXDDRUG))) %>%
    group_by(SEQN) %>%
    summarise(
      n_meds      = sum(!is.na(drug) & drug != ""),
      has_antidep = any(drug %in% antidep_names, na.rm = TRUE),
      has_anxio   = any(drug %in% anxiolytic_names, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(has_antidep_or_anxio = has_antidep | has_anxio)
}

build_adults <- function(cyc_id, cyc_data) {
  if (is.null(cyc_data$DEMO) || is.null(cyc_data$DPQ)) return(NULL)

  a <- cyc_data$DEMO %>%
    filter(RIDAGEYR >= 20) %>%
    inner_join(cyc_data$DPQ, by = "SEQN")

  # Coerce PHQ-9 numeric, set 7/9 to NA
  for (col in phq_items) {
    if (col %in% names(a)) {
      v <- a[[col]]
      if (is.factor(v)) v <- as.character(v)
      v <- suppressWarnings(as.numeric(v))
      v[v %in% c(7, 9)] <- NA
      a[[col]] <- v
    }
  }
  a$phq9 <- rowSums(a[, phq_items], na.rm = FALSE)
  a$depressed <- a$phq9 >= 10
  a$cycle <- cyc_id

  # Single-row covariates first
  for (cov in c("BMX", "SMQ", "PAQ", "HUQ")) {
    if (!is.null(cyc_data[[cov]])) {
      a <- left_join(a, cyc_data[[cov]], by = "SEQN", suffix = c("", paste0(".", cov)))
    }
  }

  # Multi-row RXQ_RX: collapse first
  rx_summary <- summarize_rx(cyc_data$RXQ_RX)
  if (!is.null(rx_summary)) {
    a <- left_join(a, rx_summary, by = "SEQN")
    a$has_antidep[is.na(a$has_antidep)] <- FALSE
    a$has_anxio[is.na(a$has_anxio)] <- FALSE
    a$has_antidep_or_anxio[is.na(a$has_antidep_or_anxio)] <- FALSE
    a$n_meds[is.na(a$n_meds)] <- 0
  }
  a
}

adults_by_cycle <- lapply(names(all_data), function(cid) build_adults(cid, all_data[[cid]]))
names(adults_by_cycle) <- names(all_data)

# --- 6. Summary table ---
cat("\n========================================\n")
cat(" DAY 2 SUMMARY\n")
cat("========================================\n")
for (cid in names(adults_by_cycle)) {
  a <- adults_by_cycle[[cid]]
  if (is.null(a)) next
  n_total <- nrow(a)
  n_phq <- sum(!is.na(a$phq9))
  n_dep <- sum(a$depressed, na.rm = TRUE)
  prev <- mean(a$depressed, na.rm = TRUE) * 100
  cat(sprintf("  %s : adults=%d, valid PHQ-9=%d, dep cases=%d, prev=%.2f%%\n",
              cycles[[cid]]$label, n_total, n_phq, n_dep, prev))
}

# --- 7. Combine all cycles (only common columns) ---
common_cols <- Reduce(intersect, lapply(adults_by_cycle, names))
cat("\nCommon columns across cycles:", length(common_cols), "\n")

adults_all <- bind_rows(lapply(adults_by_cycle, function(a) a[, common_cols, drop = FALSE]))
cat("Combined adults: n =", nrow(adults_all), "\n")
cat("  Total valid PHQ-9:", sum(!is.na(adults_all$phq9)), "\n")
cat("  Total depression cases:", sum(adults_all$depressed, na.rm = TRUE), "\n")
cat("  Pooled crude prevalence:",
    round(mean(adults_all$depressed, na.rm = TRUE) * 100, 2), "%\n")

# --- 8. Save ---
tryCatch({
  if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
  saveRDS(adults_all, "data/processed/adults_all_cycles.rds")
  saveRDS(all_data, "data/raw/all_cycles_raw.rds")
  cat("\nSaved data/processed/adults_all_cycles.rds\n")
}, error = function(e) {
  cat("Save skipped (non-critical):", conditionMessage(e), "\n")
})

cat("\nDay 2 done. Next: FPED + creatine calculation (Day 3).\n")
