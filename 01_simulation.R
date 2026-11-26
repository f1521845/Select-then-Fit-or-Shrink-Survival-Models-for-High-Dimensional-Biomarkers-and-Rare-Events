# =============================================================================
#  01_simulation.R
#  Select-then-fit or shrink? Simulation study (protocol Sections 4, 6, 7)
#  and the practitioner guide (which model for which n, event rate and p).
# =============================================================================
#
#  One script, three modes (set MODE below or with the SURVSIM_MODE variable):
#
#    "demo"    Small local run (6 scenarios x 3 replications, short MCMC) that
#              exercises every step, from simulation to the guide. Minutes to
#              about an hour on a laptop. Results are NOT for the paper.
#    "local"   The full design (83 scenarios) run sequentially on this machine.
#              Only realistic for a subset: set SCENARIOS below.
#    "collect" Summarise the full design after the HPC run
#              (hpc/run_task.R + hpc/submit.slurm write results/main/).
#
#  Outputs (results/<mode>/paper/):
#    summary_metrics.csv     mean and Monte Carlo SE per scenario x model x measure
#    summary_coefficients.csv bias, RMSE, coverage, exaggeration of log-HRs
#    diagnostics.csv         MCMC diagnostics and failures
#    hypotheses.csv          H1-H6: pre-specified criterion, scenarios supporting
#    fig_*.png               figures for the paper
#    practitioner_guide.csv / .md, fig_guide_heatmap.png
#
#  Run from the project root:  Rscript 01_simulation.R
# =============================================================================

MODE <- Sys.getenv("SURVSIM_MODE", "demo")
SCENARIOS <- NULL          # "local" mode: optional vector of scenario ids to run

source("R/design.R")
source("R/run.R")
source("R/guide.R")
suppressPackageStartupMessages(library(ggplot2))

out_root <- file.path("results", MODE)
rep_dir <- if (MODE == "collect") "results/main" else file.path(out_root, "reps")
sc_dir <- if (MODE == "collect") "results/scenarios" else file.path(out_root, "scenarios")
paper_dir <- file.path(out_root, "paper")
for (d in c(rep_dir, sc_dir, paper_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

# -----------------------------------------------------------------------------
# 1. Design
# -----------------------------------------------------------------------------
design <- design_table()

if (MODE == "demo") {
  # A small cross-section of the design: low/high EPV, linear/nonlinear, high-dim,
  # and one prior-grid scenario, with reduced settings.
  keep <- with(design,
    (block == "A" & event_frac == 0.3 & dgm == "sparse" & n %in% c(200, 1000)) |
    (block == "A" & event_frac == 0.3 & dgm == "nonlinear" & n == 1000 & p == 10) |
    (block == "B" & n == 200 & p == 500 & event_frac == 0.4 & dgm == "sparse"))
  design <- design[keep, ]
  design$reps_main <- 3; design$reps_bayes <- 3; design$n_test <- 5000
  stan_args <- list(chains = 2, iter = 1000, warmup = 500, cores = 2)
  xgb_configs <- 5
  bart_args <- list(ndpost = 300, nskip = 150)
} else {
  if (MODE == "local" && !is.null(SCENARIOS)) design <- design[design$id %in% SCENARIOS, ]
  stan_args <- list(chains = 4, iter = 2000, warmup = 1000,
                    cores = min(4, parallel::detectCores()))
  xgb_configs <- 50
  bart_args <- list()
}
write.csv(design, file.path(paper_dir, "design.csv"), row.names = FALSE)
log_msg(MODE, " mode: ", nrow(design), " scenarios")

# -----------------------------------------------------------------------------
# 2. Calibrate each scenario (event fraction, evaluation times t*)
# -----------------------------------------------------------------------------
scenarios <- lapply(seq_len(nrow(design)), function(i) {
  f <- file.path(sc_dir, sprintf("sc_%03d.rds", design$id[i]))
  if (file.exists(f)) return(readRDS(f))
  if (MODE == "collect") stop("Scenario file missing: ", f, " (run hpc/run_task.R calibrate)")
  sc <- calibrate(scenario_from_row(design[i, ]))
  saveRDS(sc, f); sc
})
names(scenarios) <- design$id

# -----------------------------------------------------------------------------
# 3. Run replications (demo / local). Existing replications are skipped.
# -----------------------------------------------------------------------------
if (MODE != "collect") {
  for (i in seq_len(nrow(design))) {
    row <- design[i, ]; sc <- scenarios[[as.character(row$id)]]
    for (r in seq_len(row$reps_main)) {
      f <- file.path(rep_dir, sprintf("s%03d_r%03d.rds", row$id, r))
      if (file.exists(f)) next
      models <- ALL_MODELS
      if (isTRUE(row$variants)) models <- c(models, PRIOR_VARIANTS)
      if (r > row$reps_bayes) models <- setdiff(models, BAYES_MODELS)
      log_msg(sprintf("scenario %d (%s: n=%d, p=%d, events=%.0f%%, %s) rep %d",
                      row$id, row$block, row$n, row$p, 100 * row$event_frac, row$dgm, r))
      res <- run_rep(sc, r, models = models, stan_args = stan_args,
                     xgb_configs = xgb_configs, bart_args = bart_args, verbose = FALSE)
      saveRDS(res, f)
    }
  }
}

# -----------------------------------------------------------------------------
# 4. Collect
# -----------------------------------------------------------------------------
files <- list.files(rep_dir, "^s\\d+_r\\d+\\.rds$", full.names = TRUE)
files <- files[as.integer(sub("^s(\\d+)_r.*", "\\1", basename(files))) %in% design$id]
R <- lapply(files, readRDS)
M <- do.call(base::rbind, lapply(R, `[[`, "metrics"))
C <- do.call(base::rbind, lapply(R, `[[`, "coefs"))
D <- do.call(base::rbind, lapply(R, `[[`, "diags"))
I <- do.call(base::rbind, lapply(R, `[[`, "info"))
log_msg(length(files), " replication files collected")

big <- function(x) formatC(x, format = "d", big.mark = ",")
lab <- with(design, setNames(sprintf("%s%d: n=%s, p=%s, %g%% events, %s", block, id,
                                     big(n), big(p), 100 * event_frac, dgm), id))
model_lab <- function(m) ifelse(m %in% names(MODEL_LABELS), MODEL_LABELS[m], m)

# -----------------------------------------------------------------------------
# 5. Summaries with Monte Carlo standard errors
# -----------------------------------------------------------------------------
mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
S <- aggregate(value ~ scenario + model + metric + t_index,
               transform(M, t_index = ifelse(is.na(t_index), 0L, t_index)),
               function(x) c(mean = mean(x), mcse = mcse(x), n = length(x)))
S <- cbind(S[, 1:4], as.data.frame(S$value))
S <- merge(S, design[, c("id", "block", "n", "p", "event_frac", "dgm", "rho", "censoring")],
           by.x = "scenario", by.y = "id")
write.csv(S, file.path(paper_dir, "summary_metrics.csv"), row.names = FALSE)

coef_summary <- NULL
if (!is.null(C) && nrow(C)) {
  C$type <- ifelse(C$true != 0, "true", "null")
  C$cover <- ifelse(is.na(C$lo), NA, C$lo <= C$true & C$true <= C$hi)
  C$exag <- ifelse(C$true != 0, abs(C$est) / abs(C$true), NA)
  coef_summary <- do.call(base::rbind, lapply(split(C, list(C$scenario, C$model, C$type), drop = TRUE), function(d)
    data.frame(scenario = d$scenario[1], model = d$model[1], type = d$type[1],
               bias = mean(d$est - d$true), rmse = sqrt(mean((d$est - d$true)^2)),
               coverage = mean(d$cover), coverage_mcse = mcse(as.numeric(d$cover)),
               exaggeration = mean(d$exag), n = nrow(d))))
  write.csv(coef_summary, file.path(paper_dir, "summary_coefficients.csv"), row.names = FALSE)
}
if (!is.null(D) && nrow(D)) write.csv(D, file.path(paper_dir, "diagnostics.csv"), row.names = FALSE)

# convenience accessor: mean / mcse of one measure
get_m <- function(metric, t = 0) {
  s <- S[S$metric == metric & S$t_index == t, ]
  s
}

# -----------------------------------------------------------------------------
# 6. Hypotheses H1-H6 (pre-specified criteria; protocol Section 2)
# -----------------------------------------------------------------------------
verdict <- function(k, n) if (n == 0) "not testable in this run" else
  sprintf("%s (%d of %d scenarios)", if (k / n > 0.5) "supported" else "not supported", k, n)
coxtype <- c("S1_cox", "S2_ridge", "S3_bayes_normal", "S4_lasso_enet", "S5_bayes_horseshoe")
H <- list()

# H1: large n, low dimension, linear PH -> Cox-type models within 2 MCSE on IBS
ibs <- get_m("ibs")
sc1 <- design$id[design$block == "A" & design$n >= 1000 & design$dgm == "sparse"]
k <- sum(vapply(sc1, function(s) {
  x <- ibs[ibs$scenario == s & ibs$model %in% coxtype, ]
  nrow(x) > 1 && diff(range(x$mean)) <= 2 * max(x$mcse)
}, TRUE))
H$H1 <- c("Large n, low p: Cox-type models predict equally well",
          "IBS range across Cox-type models <= 2 x MCSE", verdict(k, length(sc1)))

# H2: small n or few events -> normal-prior and ridge better calibrated than Cox
cs <- get_m("calib_slope")
sc2 <- design$id[design$block == "A" & design$dgm == "sparse" & (design$n == 200 | design$event_frac == 0.1)]
k2 <- 0; n2 <- 0
for (s in sc2) {
  x <- cs[cs$scenario == s, ]; cox <- x[x$model == "S1_cox", ]
  if (!nrow(cox)) next                          # Cox not estimable: counts as support below
  n2 <- n2 + 1
  better <- all(abs(x$mean[x$model %in% c("S2_ridge", "S3_bayes_normal")] - 1) < abs(cox$mean - 1))
  k2 <- k2 + better
}
H$H2 <- c("Small n / few events: shrinkage better calibrated than Cox",
          "|calibration slope - 1| smaller for ridge and Bayes-normal than for Cox (scenarios where Cox is estimable)",
          verdict(k2, n2))

# H3: sparse high-dim -> horseshoe matches/beats lasso on MCC, coverage >= 0.925
mcc <- get_m("mcc")
sc3 <- design$id[design$block %in% c("B", "C") & design$dgm == "sparse"]
k3 <- 0; n3 <- 0
for (s in sc3) {
  a <- mcc[mcc$scenario == s & mcc$model == "S5_bayes_horseshoe", ]
  b <- mcc[mcc$scenario == s & mcc$model == "S4_lasso_enet", ]
  if (!nrow(a) || !nrow(b) || is.null(coef_summary)) next
  cv <- coef_summary[coef_summary$scenario == s & coef_summary$model == "S5_bayes_horseshoe" &
                       coef_summary$type == "true", "coverage"]
  n3 <- n3 + 1
  k3 <- k3 + (a$mean >= b$mean - 2 * sqrt(a$mcse^2 + b$mcse^2) && length(cv) && cv >= 0.925)
}
H$H3 <- c("Sparse high-dim: horseshoe >= lasso in selection, near-nominal coverage",
          "MCC(horseshoe) >= MCC(lasso) - 2 MCSE and coverage of true log-HRs >= 0.925", verdict(k3, n3))

# H4: nonlinear / non-PH -> flexible models beat linear ones on IBS
sc4 <- design$id[design$dgm %in% c("nonlinear", "nonph")]
k4 <- sum(vapply(sc4, function(s) {
  x <- ibs[ibs$scenario == s, ]
  fl <- x[x$model %in% c("S6_rsf", "S7_xgb", "S9_bart"), ]; li <- x[x$model %in% coxtype, ]
  nrow(fl) && nrow(li) && min(fl$mean) < min(li$mean) - 2 * max(c(fl$mcse, li$mcse))
}, TRUE))
H$H4 <- c("Non-linear / non-PH: RSF, XGBoost or BART beat linear models",
          "best flexible IBS < best linear IBS - 2 MCSE", verdict(k4, length(sc4)))

# H5: prior sensitivity by EPV (Block D prior-grid scenarios)
sc5 <- design$id[design$variants]
h5 <- do.call(base::rbind, lapply(sc5, function(s) {
  x <- ibs[ibs$scenario == s & ibs$model %in% c("S3_bayes_normal", "S5_bayes_horseshoe", PRIOR_VARIANTS), ]
  if (!nrow(x)) return(NULL)
  dsg <- design[design$id == s, ]
  data.frame(scenario = s, epv = dsg$n * dsg$event_frac / dsg$p,
             rel_range = diff(range(x$mean)) / min(x$mean))
}))
if (!is.null(h5) && nrow(h5)) {
  ok <- (h5$epv < 10 & h5$rel_range > 0.05) | (h5$epv > 20 & h5$rel_range <= 0.05)
  H$H5 <- c("Prior-sensitive when EPV < 10, robust when EPV > 20",
            "relative IBS range across the prior grid > 5% (EPV < 10) or <= 5% (EPV > 20)",
            verdict(sum(ok), nrow(h5)))
} else H$H5 <- c("Prior-sensitive when EPV < 10, robust when EPV > 20", "-", verdict(0, 0))

# H6: select-then-fit: exaggerated HRs, under-coverage, optimistic apparent C
sc6 <- design$id[design$block %in% c("B", "C")]
k6 <- 0; n6 <- 0
if (!is.null(coef_summary)) for (s in sc6) {
  st <- coef_summary[coef_summary$scenario == s & coef_summary$model == "S8_select_fit", ]
  if (!nrow(st)) next
  app <- get_m("apparent_c"); tru <- get_m("harrell_c")
  opt <- app$mean[app$scenario == s & app$model == "S8_select_fit"] -
         tru$mean[tru$scenario == s & tru$model == "S8_select_fit"]
  exag <- st$exaggeration[st$type == "true"]
  cov_all <- weighted.mean(st$coverage, st$n)
  n6 <- n6 + 1
  k6 <- k6 + (length(exag) && exag > 1 && cov_all < 0.90 && length(opt) && opt > 0)
}
H$H6 <- c("Select-then-fit: exaggerated HRs, under-coverage, optimism",
          "exaggeration ratio > 1, coverage of selected log-HRs < 0.90, apparent C > true C",
          verdict(k6, n6))

hyp <- data.frame(hypothesis = names(H), statement = sapply(H, `[`, 1),
                  criterion = sapply(H, `[`, 2), result = sapply(H, `[`, 3), row.names = NULL)
write.csv(hyp, file.path(paper_dir, "hypotheses.csv"), row.names = FALSE)
print(hyp[, c("hypothesis", "result")], row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Figures
# -----------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10))
pal_models <- function(x) model_lab(x)
model_scale <- list(scale_colour_manual(values = MODEL_COLOURS),
                    guides(colour = guide_legend(nrow = 3, byrow = TRUE)))

# 7a. Prediction: MSE of S(t*50) relative to the best model in each scenario
mse <- get_m("mse_S", 2)
mse <- do.call(base::rbind, lapply(split(mse, mse$scenario), function(d) transform(d, rel = mean / min(mean))))
mse$label <- factor(lab[as.character(mse$scenario)], levels = lab)
p1 <- ggplot(mse[mse$model != "S8leak_select_fit", ],
             aes(rel, label, colour = pal_models(model))) +
  geom_vline(xintercept = 1, colour = "grey60") + geom_point(size = 2) +
  scale_x_log10() + model_scale +
  labs(x = "MSE of S(t*50) relative to the best model (log scale; 1 = best)", y = NULL, colour = NULL,
       title = "Prediction accuracy across scenarios") +
  theme(legend.position = "bottom")
ggsave(file.path(paper_dir, "fig_prediction_relative_mse.png"), p1,
       width = 10, height = 2.5 + 0.3 * length(unique(mse$label)), dpi = 150)

# 7b. Calibration slope (1 = ideal)
cs$label <- factor(lab[as.character(cs$scenario)], levels = lab)
p2 <- ggplot(cs[cs$model != "S8leak_select_fit", ],
             aes(mean, label, colour = pal_models(model))) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbar(aes(xmin = mean - 2 * mcse, xmax = mean + 2 * mcse), width = 0, orientation = "y") +
  geom_point(size = 2) + model_scale +
  labs(x = "Calibration slope at t*50 (1 = ideal; < 1 overfitting; > 1 over-shrinkage)", y = NULL,
       colour = NULL, title = "Calibration across scenarios") + theme(legend.position = "bottom")
ggsave(file.path(paper_dir, "fig_calibration_slope.png"), p2,
       width = 10, height = 2.5 + 0.3 * length(unique(cs$label)), dpi = 150)

# 7c. Coverage of 95% intervals for true log-HRs (inference)
if (!is.null(coef_summary)) {
  cv <- coef_summary[coef_summary$type == "true" & !is.na(coef_summary$coverage), ]
  if (nrow(cv)) {
    cv$label <- factor(lab[as.character(cv$scenario)], levels = lab)
    p3 <- ggplot(cv, aes(coverage, label, colour = pal_models(model))) +
      annotate("rect", xmin = 0.925, xmax = 0.975, ymin = -Inf, ymax = Inf, alpha = 0.12) +
      geom_vline(xintercept = 0.95, colour = "grey60") +
      geom_errorbar(aes(xmin = coverage - 2 * coverage_mcse, xmax = coverage + 2 * coverage_mcse), width = 0, orientation = "y") +
      geom_point(size = 2) + model_scale +
      labs(x = "Coverage of 95% intervals for true log hazard ratios (band: 92.5-97.5%)", y = NULL,
           colour = NULL, title = "Honest uncertainty? Interval coverage") + theme(legend.position = "bottom")
    ggsave(file.path(paper_dir, "fig_coverage.png"), p3,
           width = 10, height = 2.5 + 0.3 * length(unique(cv$label)), dpi = 150)
  }
}

# 7d. H6: apparent vs split-sample vs true C for select-then-fit
opt <- S[S$model %in% c("S8_select_fit", "S8leak_select_fit") &
           S$metric %in% c("apparent_c", "local_split_c", "harrell_c"), ]
if (nrow(opt)) {
  opt$label <- factor(lab[as.character(opt$scenario)], levels = lab)
  opt$metric <- factor(opt$metric, c("harrell_c", "local_split_c", "apparent_c"),
                       c("True (large test set)", "80/20 split estimate", "Apparent (training)"))
  p4 <- ggplot(opt, aes(mean, label, colour = metric, shape = model_lab(model))) +
    geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
    labs(x = "Harrell's C", y = NULL, colour = NULL, shape = NULL,
         title = "Select-then-fit: reported vs true discrimination") + theme(legend.position = "bottom")
  ggsave(file.path(paper_dir, "fig_h6_optimism.png"), p4,
         width = 9, height = 2 + 0.35 * length(unique(opt$label)), dpi = 150)
}

# -----------------------------------------------------------------------------
# 8. Practitioner guide
# -----------------------------------------------------------------------------
main <- design[design$block %in% c("A", "B", "C"), ]
guide <- build_guide(M[M$scenario %in% main$id, ], C[C$scenario %in% main$id, , drop = FALSE], main)
write.csv(guide, file.path(paper_dir, "practitioner_guide.csv"), row.names = FALSE)
md <- guide_markdown(guide, title = if (MODE == "demo") "Practitioner guide (DEMO RUN: not for publication)"
                                    else "Practitioner guide")
writeLines(md, file.path(paper_dir, "practitioner_guide.md"))
guide_heatmap(guide, file.path(paper_dir, "fig_guide_heatmap.png"))
cat(md, sep = "\n")

log_msg("done. Outputs in ", paper_dir)
