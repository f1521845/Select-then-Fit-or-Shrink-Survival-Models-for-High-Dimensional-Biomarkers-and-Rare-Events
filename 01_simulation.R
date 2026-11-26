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
#    hypotheses.csv          H1-H6 (H3a, H3b): registered decision rules, verdict per hypothesis
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
# 6. Hypotheses: registered decision rules (OSF registration, Section 3)
# -----------------------------------------------------------------------------
# Paired comparisons use per-replication differences between models fitted to
# the same simulated datasets; MCSE of a difference = SD(differences)/sqrt(R).
# A scenario "fails" for a comparison if any model involved has fewer than 80%
# of the replications available for that scenario.
# Support: criterion met in >= 70% of designated scenarios (H6: each of parts
# a-c); < 30% = contradicted; otherwise mixed.
coxtype  <- c("S1_cox", "S2_ridge", "S3_bayes_normal", "S4_lasso_enet", "S5_bayes_horseshoe")
flexible <- c("S6_rsf", "S7_xgb", "S9_bart")
design$epv <- design$n * design$event_frac / design$p
nrep_sc <- tapply(M$rep, M$scenario, function(r) length(unique(r)))

# per-replication values of one metric for one model and scenario
pr <- function(metric, model, s, t = NA) {
  x <- M[M$scenario == s & M$model == model & M$metric == metric, ]
  x <- if (is.na(t)) x[is.na(x$t_index), ] else x[!is.na(x$t_index) & x$t_index == t, ]
  setNames(x$value, x$rep)
}
enough <- function(v, s) length(v) >= 0.8 * nrep_sc[[as.character(s)]]
# paired difference a - b: mean and MCSE over common replications
pdiff <- function(a, b) {
  r <- intersect(names(a), names(b)); d <- a[r] - b[r]
  c(mean = mean(d), mcse = sd(d) / sqrt(length(d)), n = length(d))
}
# "clearly greater than zero": difference > 2 x MCSE
gt0 <- function(d) is.finite(d[["mean"]]) && d[["mean"]] > 2 * d[["mcse"]]

# per-replication coefficient summaries
coef_rep <- function(model, s, fun) {
  if (is.null(C)) return(numeric(0))
  x <- C[C$scenario == s & C$model == model, ]
  if (!nrow(x)) return(numeric(0))
  v <- vapply(split(x, x$rep), fun, numeric(1))
  v[is.finite(v)]
}
is_sel <- function(d) if (all(is.na(d$lo))) !is.na(d$est) & d$est != 0 else d$lo > 0 | d$hi < 0
mse_beta    <- function(d) mean((d$est - d$true)^2)
bias_sel    <- function(d) { k <- is_sel(d) & d$true != 0                 # signed: > 0 = away from zero
                             if (!any(k)) NA else mean(sign(d$true[k]) * (d$est[k] - d$true[k])) }
cover_sel   <- function(d) { k <- is_sel(d) & !is.na(d$lo)
                             if (!any(k)) NA else mean(d$lo[k] <= d$true[k] & d$true[k] <= d$hi[k]) }
cover_true  <- function(d) { k <- d$true != 0 & !is.na(d$lo)
                             if (!any(k)) NA else mean(d$lo[k] <= d$true[k] & d$true[k] <= d$hi[k]) }

verdict <- function(met, label = "") {
  met <- met[!is.na(met)]; n <- length(met)
  if (n == 0) return("not testable in this run")
  f <- mean(met)
  sprintf("%s%s (%d of %d scenarios, %.0f%%)", label,
          if (f >= 0.7) "supported" else if (f < 0.3) "contradicted" else "mixed", sum(met), n, 100 * f)
}
H <- list(); detail <- list()

# H1 (secondary): Block A sparse, n = 5,000; all pairwise IBS differences
# among S1-S5 within +/-0.005 including their Monte Carlo 95% interval
sc <- design$id[design$block == "A" & design$dgm == "sparse" & design$n == 5000]
met <- vapply(sc, function(s) {
  v <- lapply(coxtype, function(m) pr("ibs", m, s))
  if (!all(vapply(v, enough, TRUE, s = s))) return(FALSE)
  pr2 <- combn(length(v), 2)
  all(apply(pr2, 2, function(ij) { d <- pdiff(v[[ij[1]]], v[[ij[2]]])
    abs(d[["mean"]]) + 2 * d[["mcse"]] < 0.005 }))
}, TRUE)
H$H1 <- c("Large n, low dimension: Cox-type models predict about equally well",
          "every pairwise IBS difference among S1-S5 within +/-0.005 incl. 2 x MCSE", verdict(met))

# H2 (core): Block A sparse with 5 <= EPV <= 20 (6 scenarios); S2 and S3 both
# have calibration slope closer to 1 and lower RMSE(beta) than S1
sc <- design$id[design$block == "A" & design$dgm == "sparse" & design$epv >= 5 & design$epv <= 20]
met <- vapply(sc, function(s) {
  dev <- function(m) abs(pr("calib_slope", m, s) - 1)
  mb  <- function(m) coef_rep(m, s, mse_beta)
  if (!all(vapply(c("S1_cox", "S2_ridge", "S3_bayes_normal"), function(m) enough(pr("ibs", m, s), s), TRUE)))
    return(FALSE)
  all(vapply(c("S2_ridge", "S3_bayes_normal"), function(m)
    gt0(pdiff(dev("S1_cox"), dev(m))) && gt0(pdiff(mb("S1_cox"), mb(m))), TRUE))
}, TRUE)
H$H2 <- c("Few events per variable: ridge and normal-prior Cox beat classical Cox",
          "S2 and S3: |calibration slope - 1| and MSE(beta) lower than S1 by > 2 x MCSE (EPV 5-20)",
          verdict(met))

# H3a / H3b (core): sparse scenarios in Blocks B and C (12)
sc <- design$id[design$block %in% c("B", "C") & design$dgm == "sparse"]
met <- vapply(sc, function(s) {
  a <- pr("mcc", "S5_bayes_horseshoe", s); b <- pr("mcc", "S4_lasso_enet", s)
  if (!enough(a, s) || !enough(b, s)) return(FALSE)
  d <- pdiff(a, b); d[["mean"]] - 2 * d[["mcse"]] > -0.05
}, TRUE)
H$H3a <- c("Sparse high-dim: horseshoe selects at least as well as lasso/elastic net",
           "MCC(S5) - MCC(S4) > -0.05, including 2 x MCSE", verdict(met))
met <- vapply(sc, function(s) {
  cv <- coef_rep("S5_bayes_horseshoe", s, cover_true)
  if (!enough(cv, s)) return(FALSE)
  m <- mean(cv); m >= 0.925 && m <= 0.975
}, TRUE)
H$H3b <- c("Sparse high-dim: horseshoe intervals close to nominal",
           "coverage of true non-zero log HRs between 92.5% and 97.5%", verdict(met))

# H4 (secondary): Block A nonlinear and non-PH (36); best flexible model has
# lower IBS than best linear model, and BART coverage of S(t*50) in [92.5, 97.5]%
sc <- design$id[design$block == "A" & design$dgm %in% c("nonlinear", "nonph")]
met <- vapply(sc, function(s) {
  mean_ibs <- function(ms) { v <- vapply(ms, function(m) { x <- pr("ibs", m, s)
                                if (enough(x, s)) mean(x) else NA }, 1); v[!is.na(v)] }
  fl <- mean_ibs(flexible); li <- mean_ibs(coxtype)
  if (!length(fl) || !length(li)) return(FALSE)
  d <- pdiff(pr("ibs", names(which.min(li)), s), pr("ibs", names(which.min(fl)), s))
  cv <- pr("cover_S", "S9_bart", s, t = 2)
  gt0(d) && enough(cv, s) && mean(cv) >= 0.925 && mean(cv) <= 0.975
}, TRUE)
H$H4 <- c("Nonlinear / non-PH: flexible models predict better; BART near-nominal",
          "best flexible IBS < best linear IBS by > 2 x MCSE and BART coverage of S(t*) in [92.5, 97.5]%",
          verdict(met))

# H5 (core): Block D prior-grid references; largest relative change in RMSE(beta)
# or IBS across the prior grid >= 10% when EPV < 10 and < 5% when EPV > 20
sc <- design$id[design$variants]
fam <- list(S3_bayes_normal = c("S3v_fixed1", "S3v_scaled", "S3v_hyper2.5"),
            S5_bayes_horseshoe = c("S5v_p0_1", "S5v_p0_5", "S5v_p0_50"))
met <- vapply(sc, function(s) {
  ch <- unlist(lapply(names(fam), function(def) {
    base_ibs <- mean(pr("ibs", def, s)); base_rm <- sqrt(mean(coef_rep(def, s, mse_beta)))
    unlist(lapply(fam[[def]], function(v) c(
      abs(mean(pr("ibs", v, s)) - base_ibs) / base_ibs,
      abs(sqrt(mean(coef_rep(v, s, mse_beta)))  - base_rm) / base_rm)))
  }))
  ch <- ch[is.finite(ch)]
  if (!length(ch)) return(NA)
  epv <- design$epv[design$id == s]
  if (epv < 10) max(ch) >= 0.10 else if (epv > 20) max(ch) < 0.05 else NA
}, TRUE)
H$H5 <- c("Prior matters when EPV < 10, not when EPV > 20",
          "max relative change in RMSE(beta) or IBS across prior grid >= 10% (EPV < 10) / < 5% (EPV > 20)",
          verdict(met))

# H6 (core): all Blocks B and C (20). (a) bias of selected true-signal log HRs
# away from zero; (b) S8 coverage of selected coefficients < 90%; (c) apparent
# and internal-split C exceed true C; and S4/S5 reduce bias and optimism, S5
# has higher coverage. Each part must be met in >= 70% of scenarios.
sc <- design$id[design$block %in% c("B", "C")]
opt <- function(m, s, which = "apparent_c") { a <- pr(which, m, s); t <- pr("harrell_c", m, s)
  r <- intersect(names(a), names(t)); setNames(a[r] - t[r], r) }
absd <- function(a, b) {   # paired difference of absolute values |a| - |b|, per replication
  r <- intersect(names(a), names(b)); d <- abs(a[r]) - abs(b[r])
  c(mean = mean(d), mcse = sd(d) / sqrt(length(d)), n = length(d)) }
parts <- t(vapply(sc, function(s) {
  ok <- all(vapply(c("S8_select_fit", "S4_lasso_enet", "S5_bayes_horseshoe"),
                   function(m) enough(pr("harrell_c", m, s), s), TRUE))
  if (!ok) return(c(a = FALSE, b = FALSE, c = FALSE))
  b8 <- coef_rep("S8_select_fit", s, bias_sel)
  b4 <- coef_rep("S4_lasso_enet", s, bias_sel); b5 <- coef_rep("S5_bayes_horseshoe", s, bias_sel)
  bias8 <- c(mean = mean(b8), mcse = sd(b8) / sqrt(length(b8)))
  pa <- gt0(bias8) && gt0(absd(b8, b4)) && gt0(absd(b8, b5))
  c8 <- coef_rep("S8_select_fit", s, cover_sel); c5 <- coef_rep("S5_bayes_horseshoe", s, cover_sel)
  pb <- length(c8) > 1 && mean(c8) + 2 * sd(c8) / sqrt(length(c8)) < 0.90 && gt0(pdiff(c5, c8))
  o8 <- opt("S8_select_fit", s); o8s <- opt("S8_select_fit", s, "local_split_c")
  o4 <- opt("S4_lasso_enet", s); o5 <- opt("S5_bayes_horseshoe", s)
  pc <- gt0(c(mean = mean(o8), mcse = sd(o8) / sqrt(length(o8)))) &&
        gt0(c(mean = mean(o8s), mcse = sd(o8s) / sqrt(length(o8s)))) &&
        gt0(absd(o8, o4)) && gt0(absd(o8, o5))
  c(a = pa, b = pb, c = pc)
}, logical(3)))
if (length(sc)) {
  fr <- colMeans(parts)
  H$H6 <- c("Select-then-fit: biased HRs, under-coverage, optimism; shrinkage reduces them",
            "(a) bias away from 0, (b) S8 coverage < 90%, (c) apparent and internal-split C > true C; S4/S5 better",
            if (all(fr >= 0.7)) sprintf("supported (a %.0f%%, b %.0f%%, c %.0f%% of %d scenarios)", 100*fr[1], 100*fr[2], 100*fr[3], length(sc))
            else if (all(fr < 0.3)) sprintf("contradicted (a %.0f%%, b %.0f%%, c %.0f%% of %d scenarios)", 100*fr[1], 100*fr[2], 100*fr[3], length(sc))
            else sprintf("mixed (a %.0f%%, b %.0f%%, c %.0f%% of %d scenarios)", 100*fr[1], 100*fr[2], 100*fr[3], length(sc)))
} else H$H6 <- c("Select-then-fit: biased HRs, under-coverage, optimism; shrinkage reduces them", "-",
                 "not testable in this run")

H <- H[c("H6", "H2", "H3a", "H3b", "H5", "H4", "H1")]
hyp <- data.frame(hypothesis = names(H),
                  role = c("core", "core", "core", "core", "core", "secondary", "secondary"),
                  statement = sapply(H, `[`, 1), criterion = sapply(H, `[`, 2),
                  result = sapply(H, `[`, 3), row.names = NULL)
write.csv(hyp, file.path(paper_dir, "hypotheses.csv"), row.names = FALSE)
print(hyp[, c("hypothesis", "role", "result")], row.names = FALSE)

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
cs <- get_m("calib_slope")
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
