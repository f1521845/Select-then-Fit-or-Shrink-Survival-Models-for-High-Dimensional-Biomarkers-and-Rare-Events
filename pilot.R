# Pilot (protocol Section 8): checks code and runtime ONLY.
# Pilot results must not be used to choose methods or priors on performance.
#
# Usage: Rscript pilot.R correctness   # 5 replications of 3 representative scenarios
#        Rscript pilot.R timing        # 1 replication of the slowest scenarios

source("R/run.R")
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "correctness"
out <- file.path("results", "pilot", phase)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Protocol MCMC settings: 4 chains, 1,000 warm-up + 1,000 sampling (2 cores here)
stan_full <- list(chains = 4, iter = 2000, warmup = 1000, cores = 2)

if (phase == "correctness") {
  scenarios <- list(
    make_scenario(101, n = 200,  p = 50,  event_frac = 0.3, dgm = "sparse",    block_name = "A"),
    make_scenario(102, n = 1000, p = 10,  event_frac = 0.3, dgm = "nonlinear", block_name = "A"),
    make_scenario(103, n = 200,  p = 500, event_frac = 0.4, dgm = "sparse",    block_name = "B"))
  reps <- 1:5
  for (sc in scenarios) {
    sc <- calibrate(sc)
    saveRDS(sc, file.path(out, sprintf("scenario_%d.rds", sc$id)))
    for (r in reps) {
      f <- file.path(out, sprintf("s%d_r%02d.rds", sc$id, r))
      if (file.exists(f)) next
      message(sprintf("[%s] scenario %d rep %d", format(Sys.time(), "%H:%M:%S"), sc$id, r))
      res <- run_rep(sc, r, stan_args = stan_full, xgb_configs = 10)
      saveRDS(res, f)
    }
  }
}

if (phase == "timing") {
  # Slowest scenario of Block A (all models, full settings)
  scA <- calibrate(make_scenario(201, n = 5000, p = 50, event_frac = 0.7, dgm = "nonph", block_name = "A"))
  message("Block A slowest: n=5000, p=50")
  res <- run_rep(scA, 1, stan_args = stan_full, xgb_configs = 10)
  saveRDS(res, file.path(out, "blockA_slowest.rds"))

  # Slowest scenarios of Blocks B and C: non-Bayesian models in full; Bayesian
  # models timed with a short single-chain run to extrapolate per-iteration cost.
  short <- list(chains = 1, iter = 100, warmup = 50, cores = 1)
  for (cfg in list(list(id = 202, n = 1000, p = 2000, ev = 0.4, nc = 0, blk = "B", nt = 10000),
                   list(id = 203, n = 20000, p = 2000, ev = 0.005, nc = 3, blk = "C", nt = 20000))) {
    sc <- calibrate(make_scenario(cfg$id, n = cfg$n, p = cfg$p, event_frac = cfg$ev,
                                  dgm = "sparse", n_clin = cfg$nc, block_name = cfg$blk,
                                  n_test = cfg$nt), n_cal = 5e4)
    message(sprintf("Block %s slowest: n=%d, p=%d", cfg$blk, cfg$n, cfg$p))
    fast <- c("S2_ridge", "S4_lasso_enet", "S6_rsf", "S7_xgb", "S8_select_fit")
    res <- run_rep(sc, 1, models = fast, xgb_configs = 10)
    saveRDS(res, file.path(out, sprintf("block%s_fast.rds", cfg$blk)))
    res_b <- run_rep(sc, 1, models = c("S3_bayes_normal", "S5_bayes_horseshoe"), stan_args = short)
    saveRDS(res_b, file.path(out, sprintf("block%s_bayes_short.rds", cfg$blk)))
  }
}
message("done: ", phase)
