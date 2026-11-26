# Summarise replication results: mean and Monte Carlo SE per scenario x model.
# Usage: Rscript summarise.R results/pilot/correctness

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results/pilot/correctness"
files <- list.files(dir, pattern = "^s\\d+_r\\d+\\.rds$", full.names = TRUE)
R <- lapply(files, readRDS)
M <- do.call(base::rbind, lapply(R, `[[`, "metrics"))
C <- do.call(base::rbind, lapply(R, `[[`, "coefs"))
D <- do.call(base::rbind, lapply(R, `[[`, "diags"))
I <- do.call(base::rbind, lapply(R, `[[`, "info"))

mcse <- function(x) sd(x) / sqrt(length(x))
fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)

cat("Replications per scenario:\n"); print(table(I$scenario))
cat("\nTraining events per replication (mean):\n"); print(tapply(I$n_events, I$scenario, mean))

key <- c("mse_S", "ibs", "harrell_c", "uno_c", "calib_slope", "ici", "apparent_c",
         "local_split_c", "cover_S", "tpr", "fdr", "mcc", "fit_seconds")
M2 <- M[M$metric %in% key & (is.na(M$t_index) | M$t_index == 2), ]
for (sid in sort(unique(M2$scenario))) {
  cat(sprintf("\n=== Scenario %d: mean (MCSE) over replications; S(t*) metrics at t*50 ===\n", sid))
  s <- M2[M2$scenario == sid, ]
  tab <- tapply(seq_len(nrow(s)), list(s$model, s$metric), function(i) fmt(mean(s$value[i]), mcse(s$value[i])))
  print(noquote(tab[, intersect(key, colnames(tab)), drop = FALSE]))
}

fails <- M[M$metric %in% c("failed", "not_estimable"), ]
if (nrow(fails)) { cat("\nFailures / not estimable:\n"); print(aggregate(value ~ scenario + model + metric, fails, sum)) }

if (!is.null(C) && nrow(C)) {
  C$cover <- C$lo <= C$true & C$true <= C$hi
  C$type <- ifelse(C$true != 0, "true", "null")
  cat("\n=== Coefficients (linear-PH scenarios): bias, RMSE, 95% coverage by predictor type ===\n")
  agg <- aggregate(cbind(bias = est - true, sq = (est - true)^2, cover = cover) ~ scenario + model + type,
                   C, mean, na.action = na.pass)
  agg$rmse <- sqrt(agg$sq); agg$sq <- NULL
  print(agg[order(agg$scenario, agg$type, agg$model), ], digits = 3, row.names = FALSE)
  # post-selection bias for true predictors: |estimate| vs |truth| when selected
  sel <- C[C$model %in% c("S8_select_fit", "S8b_select_bayes", "S8leak_select_fit", "S1_cox", "S3_bayes_normal", "S5_bayes_horseshoe") & C$type == "true", ]
  if (nrow(sel)) {
    cat("\nExaggeration ratio |estimate| / |true| for true predictors (1 = unbiased):\n")
    print(aggregate(ratio ~ scenario + model, transform(sel, ratio = abs(est) / abs(true)), mean), digits = 3, row.names = FALSE)
  }
}

if (!is.null(D) && nrow(D)) {
  dd <- D[D$key %in% c("max_rhat", "min_ess", "divergences"), ]
  dd$value <- as.numeric(dd$value)
  cat("\n=== MCMC diagnostics (worst over replications) ===\n")
  print(aggregate(value ~ scenario + model + key, dd, function(v) if (all(v >= 1) && max(v) < 2) max(v) else c(min = min(v), max = max(v))), digits = 4)
}
