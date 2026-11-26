# 04_DLBCL.R
# Run in RStudio: open this file and click 'Source' (or Ctrl/Cmd + Shift + S).
# Tables print to the Console, figures appear in the Plots pane, and everything
# is also saved under results/realdata/<dataset>/.

##############################################################################
# DLBCL: high-dimensional gene expression (real-data application 3)
##############################################################################
#
# Role in the study: p (7,399 genes) far exceeds n (240) and events (138). Tests H3
# (horseshoe vs lasso in sparse high-dimensional data) and H6 (select-then-fit optimism).
#
# Data: diffuse large B-cell lymphoma microarray data (SurvSet DLBCL); outcome = death.
#
# Analysis: repeated 5-fold cross-validation (10 repeats). Inside each training fold, the
# top 1,000 genes by univariate Cox score are kept, and the same screened set is used by
# every model (MCMC on 7,399 genes is impractical). Classical Cox is not estimable (p >
# n).

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

DS <- "DLBCL"
out_dir <- file.path("results", "realdata", DS, if (QUICK) "quick" else "full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- if (QUICK) list(K = 3, repeats = 1, m_imp = 1, xgb = 3, bayes_folds = 3,
                       stan = list(chains = 2, iter = 600, warmup = 300, cores = 2),
                       bart = list(ndpost = 200, nskip = 100), screen_K = 200) else
               list(K = 5, repeats = 10, m_imp = 5, xgb = 50, bayes_folds = Inf,
                    stan = list(chains = 4, iter = 2000, warmup = 1000, cores = 4),
                    bart = list(), screen_K = 1000)

MODELS <- c("S1_cox", "S2_ridge", "S3_bayes_normal", "S4_lasso_enet", "S5_bayes_horseshoe",
            "S6_rsf", "S7_xgb", "S8_select_fit", "S8b_select_bayes", "S8leak_select_fit", "S9_bart")
cat("Mode:", if (QUICK) "QUICK check" else "FULL protocol analysis", "\n")


# ---- 1. Data ---------------------------------------------------------------

if (!file.exists("data/DLBCL.csv")) system("python3 data/prepare_data.py")
if (!file.exists("data/DLBCL.csv")) stop("data/DLBCL.csv not found. The data/ folder of the repository is missing. Download it from GitHub ",
                                 "(it contains data/DLBCL.csv), or run: pip install SurvSet; python data/prepare_data.py")
raw <- read.csv("data/DLBCL.csv")
df <- raw[, grep("^num_", names(raw))]; names(df) <- sub("^num_", "", names(df))
y <- data.frame(time = raw$time, status = raw$event)
levels_ref <- list()
cat(sprintf("%d patients x %d genes; missing values: %d\n", nrow(df), ncol(df), sum(is.na(df))))


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


# ---- 3. Out-of-fold performance (repeated cross-validation, screening inside folds) ---

spec <- list(m_imp = cfg$m_imp, screen_K = cfg$screen_K, protect_regex = NULL,
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


# ---- 5. Select-then-fit and optimism (H6) ----------------------------------

# Compare the select-then-fit model's apparent C with its out-of-fold C, and the honest
# version with the leaky one. In the leak variant, gene selection saw the test fold.

# Optimism: apparent (training) C minus out-of-fold C, per model
opt <- data.frame(Model = w$label,
                  apparent_C = round(w$mean.apparent_c, 3),
                  out_of_fold_C = round(w$mean.harrell_c, 3))
opt$optimism <- round(opt$apparent_C - opt$out_of_fold_C, 3)
print(opt[!is.na(opt$apparent_C), ])


# ---- 6. Which genes are selected, and how stable is the selection? ---------

# Selection stability across outer folds (mean pairwise Jaccard; 1 = identical sets)
stab <- lapply(c(S4_lasso_enet = "S4_lasso_enet", S5_bayes_horseshoe = "S5_bayes_horseshoe",
                 S8_select_fit = "S8_select_fit"), function(m) selection_stability(cv$selected, m))
stab <- stab[!vapply(stab, is.null, TRUE)]
print(data.frame(Model = MODEL_LABELS_RD[names(stab)],
           mean_selected = sapply(stab, function(s) round(s$mean_size, 1)),
           jaccard = sapply(stab, function(s) round(s$jaccard, 2)),
           top_features = sapply(stab, function(s) {
             k <- seq_len(min(5, length(s$freq)))
             if (!length(k)) "none selected" else
               paste(sprintf("%s (%.0f%%)", names(s$freq)[k], 100 * s$freq[k]), collapse = ", ") }),
           row.names = NULL))


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

