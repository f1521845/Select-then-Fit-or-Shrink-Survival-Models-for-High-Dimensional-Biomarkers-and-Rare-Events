# 02_SUPPORT2.R  (generated from notebooks/02_SUPPORT2.ipynb; same analysis)
# Run in RStudio: open this file and click 'Source' (or Ctrl/Cmd + Shift + S).
# Tables print to the Console, figures appear in the Plots pane, and everything
# is also saved under results/realdata/<dataset>/.

##############################################################################
# SUPPORT2: large clinical cohort (real-data application 1)
##############################################################################
#
# Role in the study: large n, low dimension, common events, missing data typical of EHR.
# Tests H1 (with large n and low p, the Cox-type models perform about equally well) and
# H4 (flexible models gain only if there are non-linear effects). SUPPORT2 is also the
# base for the plasmode simulation.
#
# Data: 9,105 seriously ill hospitalised adults (SUPPORT study); outcome = time to death.
# Loaded from SurvSet (data/prepare_data.py).
#
# Analysis: repeated 5-fold cross-validation (10 repeats). Imputation, standardisation
# and tuning happen inside each training fold. Because of their cost, the Bayesian arms
# run on the first 10 outer folds, mirroring the reduced replication count used in the
# simulation.

# ---- settings ------------------------------------------------------------------
# QUICK = FALSE runs the full protocol analysis (5-fold x 10 repeats; hours to days).
# Results in QUICK mode are NOT for interpretation.
QUICK <- TRUE   # TRUE = fast check (~5-10 min). Set to FALSE for the full protocol analysis.

# Working directory must be the survsim folder (RStudio: Session > Set Working Directory >
# To Source File Location). The scripts live next to 01_simulation.R.
if (!file.exists("R/realdata.R")) stop("Set the working directory to the survsim folder first ",
                                       "(Session > Set Working Directory > To Source File Location).")
source("R/realdata.R")
suppressPackageStartupMessages(library(ggplot2))
options(width = 160)

DS <- "support2"
out_dir <- file.path("results", "realdata", DS, if (QUICK) "quick" else "full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- if (QUICK) list(K = 3, repeats = 1, m_imp = 1, xgb = 3, bayes_folds = 3,
                       stan = list(chains = 2, iter = 600, warmup = 300, cores = 2),
                       bart = list(ndpost = 200, nskip = 100), subsample = 1500) else
               list(K = 5, repeats = 10, m_imp = 5, xgb = 50, bayes_folds = 10,
                    stan = list(chains = 4, iter = 2000, warmup = 1000, cores = 4),
                    bart = list())

MODELS <- c("S1_cox", "S2_ridge", "S3_bayes_normal", "S4_lasso_enet", "S5_bayes_horseshoe",
            "S6_rsf", "S7_xgb", "S8_select_fit", "S8b_select_bayes", "S8leak_select_fit", "S9_bart")
cat("Mode:", if (QUICK) "QUICK check" else "FULL protocol analysis", "\n")


# ---- 1. Data and preprocessing decisions -----------------------------------
#
# - Excluded (leakage / not baseline): sfdm2 (functional status measured at 2 months),
#   surv2m and surv6m (survival predictions from the original SUPPORT model), and dnr
#   (DNR order, which can be written after study entry). fac_num_co duplicates num_co,
#   and dzclass is a coarser grouping of dzgroup (perfectly collinear), so both are
#   dropped.
# - Missing data: kept as a realistic feature of EHR data. Multiple imputation by chained
#   equations (mice, m = 5) is fitted within each training fold and applied to its test
#   fold. The outcome is not used in the imputation model, so test outcomes cannot leak.
# - Time scale: days converted to years.

if (!file.exists("data/support2.csv")) system("python3 data/prepare_data.py")
if (!file.exists("data/support2.csv")) stop("data/support2.csv not found. Unzip survsim_data.zip into the survsim folder ",
                                 "(it contains data/support2.csv), or run: pip install SurvSet; python data/prepare_data.py")
raw <- read.csv("data/support2.csv")
drop <- c("pid", "event", "time", "num_surv2m", "num_surv6m", "fac_sfdm2", "fac_dnr", "fac_num_co", "fac_dzclass")
df <- raw[, setdiff(names(raw), drop)]
for (v in grep("^fac_", names(df), value = TRUE)) df[[v]] <- factor(df[[v]])
names(df) <- sub("^(num|fac)_", "", names(df))
y <- data.frame(time = raw$time / 365.25, status = raw$event)
if (!is.null(cfg$subsample)) {                          # QUICK mode only
  set.seed(1); keep <- sort(sample(nrow(df), cfg$subsample)); df <- df[keep, ]; y <- y[keep, ]
}
levels_ref <- lapply(Filter(is.factor, df), levels)
miss <- sort(colMeans(is.na(df)), decreasing = TRUE)
print(data.frame(variable = names(miss), percent_missing = round(100 * miss, 1))[miss > 0, ])


# ---- 2. Descriptives -------------------------------------------------------

n <- nrow(y); events <- sum(y$status); p <- ncol(encode(df[complete.cases(df), , drop = FALSE], levels_ref))
tstar <- unname(quantile(y$time[y$status == 1], c(0.25, 0.5, 0.75)))
cat(sprintf("n = %d, events = %d (%.1f%%), candidate predictors p = %d, EPV = %.2f\n",
            n, events, 100 * events / n, p, events / p))
cat(sprintf("Median follow-up (reverse KM): %.2f years\n",
            summary(survfit(Surv(y$time, 1 - y$status) ~ 1))$table[["median"]]))
cat(sprintf("Evaluation times t* (25/50/75th percentiles of event times): %s years\n",
            paste(round(tstar, 2), collapse = ", ")))
km <- survfit(Surv(time, status) ~ 1, data = y)
plot(km, xlab = "Time (years)", ylab = "Survival", main = paste(DS, ": Kaplan-Meier estimate"),
     conf.int = TRUE, col = "#2F5D8A")
abline(v = tstar, lty = 3, col = "grey40")


# ---- What does the practitioner guide predict? -----------------------------
#
# The simulation (01_simulation.R) produces a practitioner guide that recommends a model
# for each combination of sample size, event rate and number of predictors. Before
# looking at the results, we record what the guide predicts for a dataset like this one.
# Section 7 checks whether the real data agree (secondary objective 5: generalisability).

guide <- read_guide()
if (!is.null(guide)) print(guide_lookup(guide, n, events / n, p)[, c("n", "event_frac", "p", "epv", "best_mean", "recommended")], row.names = FALSE) else
  cat("No guide yet: run 01_simulation.R (demo or collect mode) first.\n")


# ---- 3. Out-of-fold performance (repeated cross-validation) ----------------

spec <- list(m_imp = cfg$m_imp, screen_K = NULL, protect_regex = NULL,
             tstar = tstar, grid = seq(tstar[3] / 50, tstar[3], length.out = 50),
             stan_args = cfg$stan, xgb_configs = cfg$xgb, bart_args = cfg$bart, levels = levels_ref)

# Each outer fold is cached in out_dir, so an interrupted run resumes where it stopped.
t0 <- Sys.time()
cv <- run_cv(df, y, MODELS, spec, K = cfg$K, repeats = cfg$repeats,
             bayes_folds = cfg$bayes_folds, cache_dir = file.path(out_dir, "folds"))
cat("Elapsed:", format(round(Sys.time() - t0, 1)), "\n")
if (length(cv$failures)) { cat("Failures / not estimable:\n"); print(table(cv$failures)) }
saveRDS(cv, file.path(out_dir, "cv_results.rds"))


# ---- 4. Results ------------------------------------------------------------

w <- summarise_cv(cv$metrics)
tab <- format_cv_table(w)
write.csv(tab, file.path(out_dir, "table_performance.csv"), row.names = FALSE)
print(tab)

g <- plot_cv_results(cv$metrics, title = paste(DS, ": out-of-fold performance"))
ggsave(file.path(out_dir, "fig_performance.png"), g, width = 11, height = 4.5, dpi = 150)
print(g)


# ---- 5. Optimism of apparent performance -----------------------------------

# Optimism: apparent (training) C minus out-of-fold C, per model
opt <- data.frame(Model = w$label,
                  apparent_C = round(w$mean.apparent_c, 3),
                  out_of_fold_C = round(w$mean.harrell_c, 3))
opt$optimism <- round(opt$apparent_C - opt$out_of_fold_C, 3)
print(opt[!is.na(opt$apparent_C), ])


# ---- 6. Associations: do the methods agree on the hazard ratios? -----------

# Associations: hazard ratios per SD from the Cox-type models fitted to all data
cf <- fit_full_coefs(df, y, spec)
write.csv(cf, file.path(out_dir, "table_hazard_ratios.csv"), row.names = FALSE)
g2 <- plot_coefs(cf, title = paste(DS, ": hazard ratios per SD by method"))
ggsave(file.path(out_dir, "fig_hazard_ratios.png"), g2, width = 10, height = 7, dpi = 150)
print(g2)


# ---- 7. Does the practitioner guide hold? ----------------------------------

chk <- guide_check(guide, w, n, events / n, p)


# ---- 8. Interpretation (write after the FULL run) --------------------------
#
# Prompts, in protocol order. Do not interpret QUICK-mode output.
#
# - Prediction: which models are within one fold-SD of the best IBS? Do calibration
#   slopes show overfitting (< 1) or over-shrinkage (> 1)?
# - Bayesian vs frequentist pairs: ridge vs Bayesian-normal, lasso vs horseshoe. Are they
#   distinguishable here?
# - Select-then-fit (H6): how large is the optimism? How much does leakage inflate the
#   out-of-fold C?
# - Guide: is the result consistent with the practitioner guide? If not, which feature of
#   this dataset (correlation, non-linearity, censoring) might explain it?

