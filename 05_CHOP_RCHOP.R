# 05_CHOP_RCHOP.R
# Run in RStudio: open this file and click 'Source' (or Ctrl/Cmd + Shift + S).
# Tables print to the Console, figures appear in the Plots pane, and everything
# is also saved under results/realdata/<dataset>/.

##############################################################################
# CHOP <U+2192> R-CHOP: external validation (real-data application 4)
##############################################################################
#
# Role in the study: a genuine external validation across treatment eras. All tuning,
# screening and fitting use the 181 CHOP-treated patients; the frozen models are then
# evaluated once on the 233 R-CHOP-treated patients. Compares internal (cross-validated)
# estimates with external performance, the real-data check of H6's optimism gap. Also
# tests H3.
#
# Data: Lenz et al. (2008), GEO GSE10846, as distributed in the CRAN package bujar (chop:
# 181 patients, 105 deaths; rchop: 233 patients, 60 deaths; 3,833 probe sets each).

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

DS <- "chop_rchop"
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

get_bujar <- function(name) {
  f <- file.path("data", paste0(name, ".rda"))
  if (!file.exists(f)) {
    ok <- tryCatch({ download.file(sprintf("https://github.com/cran/bujar/raw/master/data/%s.rda", name), f, mode = "wb"); TRUE },
                   error = function(e) FALSE)
    if (!ok) stop("Download failed: install the bujar package or place ", name, ".rda in data/")
  }
  e <- new.env(); load(f, envir = e); get(name, envir = e)
}
chop <- get_bujar("chop"); rchop <- get_bujar("rchop")
stopifnot(identical(colnames(chop), colnames(rchop)))
to_df <- function(m) { d <- as.data.frame(m[, -(1:2)]); names(d) <- make.names(colnames(m)[-(1:2)]); d }
df <- to_df(chop);  y <- data.frame(time = chop[, 1], status = chop[, 2])
df_ext <- to_df(rchop); y_ext <- data.frame(time = rchop[, 1], status = rchop[, 2])
levels_ref <- list()
cat(sprintf("CHOP: %d patients, %d events | R-CHOP: %d patients, %d events | %d probe sets\n",
            nrow(df), sum(y$status), nrow(df_ext), sum(y_ext$status), ncol(df)))
cat(sprintf("Max follow-up: CHOP %.1f, R-CHOP %.1f years\n", max(y$time), max(y_ext$time)))
cat(sprintf("Survival times of exactly 0: CHOP %d, R-CHOP %d (moved to half the smallest positive time inside glmnet and XGBoost, which reject them)\n",
            sum(y$time == 0), sum(y_ext$time == 0)))


# ---- 2. Descriptives (development data: CHOP) ------------------------------

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


# ---- 3. Internal validation: repeated cross-validation within CHOP ---------

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


# ---- 4. Internal results ---------------------------------------------------

w <- summarise_cv(cv$metrics)
tab <- format_cv_table(w)
write.csv(tab, file.path(out_dir, "table_performance.csv"), row.names = FALSE)
print(tab)

g <- plot_cv_results(cv$metrics, title = paste(DS, ": out-of-fold performance"))
ggsave(file.path(out_dir, "fig_performance.png"), g, width = 11, height = 4.5, dpi = 150)
print(g)


# ---- 5. External validation: train on all of CHOP, test once on R-CHOP -----
#
# R-CHOP follow-up is shorter, so the evaluation times t\* for the external test are the
# 25th, 50th and 75th percentiles of R-CHOP event times. The leak variant is excluded: it
# would let gene selection see the validation cohort.

tstar_ext <- unname(quantile(y_ext$time[y_ext$status == 1], c(0.25, 0.5, 0.75)))
spec_ext <- modifyList(spec, list(tstar = tstar_ext, grid = seq(tstar_ext[3] / 50, tstar_ext[3], length.out = 50)))
f_ext <- file.path(out_dir, "external.rds")
if (file.exists(f_ext)) ext <- readRDS(f_ext) else {
  set.seed(2026)
  ext <- eval_split(df, y, df_ext, y_ext, setdiff(MODELS, "S8leak_select_fit"), spec_ext)
  ext$preds <- NULL; saveRDS(ext, f_ext)
}
we <- reshape(ext$metrics, idvar = "model", timevar = "metric", direction = "wide")
names(we) <- sub("^value\\.", "", names(we))
we$Model <- MODEL_LABELS_RD[we$model]
print(we[, c("Model", "ibs", "harrell_c", "uno_c", "td_auc", "calib_slope", "ici")])


# ---- 6. Optimism: apparent vs internal (CV) vs external discrimination -----

cmp <- merge(data.frame(model = as.character(w$model), Model = w$label,
                        apparent_C = w$mean.apparent_c, internal_CV_C = w$mean.harrell_c),
             data.frame(model = we$model, external_C = we$harrell_c), by = "model", all.x = TRUE)
cmp$internal_minus_external <- cmp$internal_CV_C - cmp$external_C
write.csv(cmp, file.path(out_dir, "table_optimism.csv"), row.names = FALSE)
print(cmp[order(-cmp$external_C), -1])


# ---- 7. Selection stability within CHOP ------------------------------------

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


# ---- 8. Does the practitioner guide hold? ----------------------------------

chk <- guide_check(guide, w, n, events / n, p)


# ---- 9. Interpretation (write after the FULL run) --------------------------
#
# - External vs internal: which models lose the most discrimination and calibration on
#   R-CHOP? The shift from CHOP to R-CHOP (rituximab era) is a real change in population
#   and treatment, so some loss is expected for every model.
# - H6: does select-then-fit show the largest gap between apparent and external
#   performance?
# - Calibration: the baseline hazard differs between eras. Separate the calibration slope
#   (relative risks transport well) from calibration in the large.

