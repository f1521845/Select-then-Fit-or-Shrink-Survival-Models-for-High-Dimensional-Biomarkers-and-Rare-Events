# Block C timing (biobank-like): n = 20,000, p = 2,000 + 3 clinical, 0.5% events.
# Full RSF tuning is too slow for this 2-core pilot machine, so RSF is timed
# with 50 trees at one setting and extrapolated linearly to the full budget
# (4 x 200 tuning trees + 500 final trees = 1,300 trees).
source("R/run.R")
out <- "results/pilot/timing"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
sc <- calibrate(make_scenario(203, n = 20000, p = 2000, event_frac = 0.005, dgm = "sparse",
                              n_clin = 3, block_name = "C", n_test = 20000), n_cal = 5e4)
message("calibrated: lambda=", signif(sc$lambda, 3), " t*=", paste(round(sc$tstar, 2), collapse = ","))
res <- run_rep(sc, 1, models = c("S2_ridge", "S4_lasso_enet", "S7_xgb", "S8_select_fit", "S8leak_select_fit"),
               xgb_configs = 10)
saveRDS(res, file.path(out, "blockC_fast.rds"))

RNGkind("L'Ecuyer-CMRG"); set.seed(1e6 * sc$id + 1); tr <- sim_data(sc, sc$n)
suppressPackageStartupMessages(library(randomForestSRC)); options(rf.cores = 2)
df <- data.frame(time = tr$time, status = tr$status, tr$X, check.names = FALSE)
t_rsf <- system.time(rfsrc(Surv(time, status) ~ ., df, ntree = 50, mtry = round(sqrt(ncol(tr$X))), nodesize = 15))[["elapsed"]]
message(sprintf("RSF 50 trees (mtry = sqrt(p)): %.1fs -> extrapolated full budget (1,300 trees): %.0fs", t_rsf, t_rsf * 1300 / 50))
saveRDS(list(rsf_50_trees = t_rsf), file.path(out, "blockC_rsf_time.rds"))

short <- list(chains = 1, iter = 100, warmup = 50, cores = 1, rerun = FALSE)
res_b <- run_rep(sc, 1, models = c("S5_bayes_horseshoe"), stan_args = short)
saveRDS(res_b, file.path(out, "blockC_bayes_short.rds"))
message("done: blockC")
