# One simulation replication (protocol Section 7, "Procedure per replication")
# ---------------------------------------------------------------------------

source("R/dgm.R"); source("R/models.R"); source("R/metrics.R")

ALL_MODELS <- c("S1_cox", "S2_ridge", "S3_bayes_normal", "S4_lasso_enet",
                "S5_bayes_horseshoe", "S6_rsf", "S7_xgb", "S8_select_fit",
                "S8b_select_bayes", "S8leak_select_fit", "S9_bart")
PRIOR_VARIANTS <- c("S3v_fixed1", "S3v_scaled", "S3v_hyper2.5", "S5v_p0_1", "S5v_p0_5", "S5v_p0_50")
# "Bayesian arms": run for the first `reps_bayes` replications only (protocol)
BAYES_MODELS <- c("S3_bayes_normal", "S5_bayes_horseshoe", "S8b_select_bayes", "S9_bart", PRIOR_VARIANTS)

#' Run one replication of one scenario.
#' @param sc calibrated scenario (see calibrate())
#' @param rep replication number; with sc$id it fixes the random-number stream
#' @param models which strategies to run
#' @param stan_args list of extra arguments for fit_bayes_cox (iter, chains, ...)
run_rep <- function(sc, rep, models = ALL_MODELS, stan_args = list(),
                    xgb_configs = 10, bart_args = list(), verbose = TRUE) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1e6 * sc$id + rep)
  tr <- sim_data(sc, sc$n)
  te <- sim_data(sc, sc$n_test)
  loc <- sim_data(sc, max(50, round(sc$n / 4)))   # "local test" split for H6
  grid <- seq(sc$tstar[3] / 50, sc$tstar[3], length.out = 50)
  times <- sort(unique(c(sc$tstar, grid)))
  S_true <- true_surv(sc, te$lp, times)
  it <- match(sc$tstar, times); ig <- match(grid, times)
  all_vars <- colnames(tr$X)
  linear_ph <- sc$dgm %in% c("sparse", "dense")

  metrics <- list(); coefs <- list(); diags <- list()
  fits <- list()
  add <- function(model, name, value, t = NA_integer_)
    metrics[[length(metrics) + 1]] <<- data.frame(model = model, metric = name,
                                                   t_index = t, value = unname(value))

  fitters <- list(
    S1_cox = function() fit_cox(tr$X, tr$time, tr$status),
    S2_ridge = function() fit_ridge(tr$X, tr$time, tr$status),
    S3_bayes_normal = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "normal"), stan_args)),
    S4_lasso_enet = function() fit_lasso(tr$X, tr$time, tr$status),
    S5_bayes_horseshoe = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "horseshoe"), stan_args)),
    S6_rsf = function() fit_rsf(tr$X, tr$time, tr$status),
    S7_xgb = function() fit_xgb(tr$X, tr$time, tr$status, n_configs = xgb_configs),
    S8_select_fit = function() fit_select_then_fit(tr$X, tr$time, tr$status, fits$S7_xgb),
    S8b_select_bayes = function() do.call(fit_select_then_fit, c(list(tr$X, tr$time, tr$status, fits$S7_xgb, second = "bayes"), stan_args)),
    S8leak_select_fit = function() fit_select_leak(tr$X, tr$time, tr$status, fits$S7_xgb, loc$X, loc$time, loc$status),
    S9_bart = function() do.call(fit_bart, c(list(tr$X, tr$time, tr$status), bart_args)),
    # prior-sensitivity variants (Block D, H5)
    S3v_fixed1 = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "normal", sigma = 1, hier = FALSE), stan_args)),
    S3v_scaled = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "normal", sigma = "scaled", hier = FALSE), stan_args)),
    S3v_hyper2.5 = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "normal", sigma = 2.5, hier = TRUE), stan_args)),
    S5v_p0_1 = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "horseshoe", p0 = 1), stan_args)),
    S5v_p0_5 = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "horseshoe", p0 = 5), stan_args)),
    S5v_p0_50 = function() do.call(fit_bayes_cox, c(list(tr$X, tr$time, tr$status, prior = "horseshoe", p0 = 50), stan_args))
  )
  # S8 variants need S7 first
  if (any(grepl("^S8", models)) && !"S7_xgb" %in% models) models <- c(models, "S7_xgb")
  models <- intersect(c(ALL_MODELS, PRIOR_VARIANTS), models)

  for (m in models) {
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(withCallingHandlers(fitters[[m]](), error = function(e) {
      calls <- sys.calls()
      assign(".last_err_calls", vapply(tail(calls, 6), function(cl) paste(deparse(cl, nlines = 1L), collapse = ""), ""), envir = globalenv())
    }), error = function(e) list(error = paste0(conditionMessage(e), " | in: ",
                                                 paste(get0(".last_err_calls", envir = globalenv(), ifnotfound = ""), collapse = " <- "))))
    fit_sec <- proc.time()[["elapsed"]] - t0
    fits[[m]] <- f
    if (!is.null(f$not_estimable)) { add(m, "not_estimable", 1); next }
    if (!is.null(f$error)) {
      add(m, "failed", 1)
      diags[[m]] <- data.frame(model = m, key = "error", value = f$error)
      if (verbose) message(sprintf("  %-20s FAILED: %s", m, f$error))
      next
    }
    t1 <- proc.time()[["elapsed"]]
    S <- tryCatch(f$surv(te$X, times), error = function(e) NULL)
    if (is.null(S)) { add(m, "failed", 1); next }
    risk <- 1 - S[, it[2]]
    # prediction
    for (j in 1:3) add(m, "mse_S", mean((S[, it[j]] - S_true[, it[j]])^2), j)
    add(m, "ibs", ibs(te$time, te$status, S[, ig], grid))
    add(m, "harrell_c", harrell_c(te$time, te$status, risk))
    add(m, "uno_c", uno_c(te$time, te$status, risk, sc$tstar[3]))
    cal <- calib_truth(S[, it[2]], S_true[, it[2]])
    add(m, "calib_slope", cal["calib_slope"]); add(m, "ici", cal["ici"])
    # optimism: apparent (training) C, and for S8 variants the local split C
    if (!is.null(f$lp_train)) add(m, "apparent_c", harrell_c(tr$time, tr$status, f$lp_train))
    if (grepl("^S8", m)) {
      Sl <- f$surv(loc$X, sc$tstar[2])
      add(m, "local_split_c", harrell_c(loc$time, loc$status, 1 - Sl[, 1]))
    }
    # interval coverage for S(t*) where intervals exist (Bayesian models, BART)
    if (!is.null(f$surv_int)) {
      ci <- f$surv_int(te$X[seq_len(min(2000, sc$n_test)), , drop = FALSE], sc$tstar)
      st <- S_true[seq_len(min(2000, sc$n_test)), it, drop = FALSE]
      for (j in 1:3) add(m, "cover_S", mean(ci$lo[, j] <= st[, j] & st[, j] <= ci$hi[, j]), j)
    }
    # coefficients (linear PH truth only)
    if (!is.null(f$coef) && linear_ph) {
      tb <- sc$truth$beta[f$coef$var]
      coefs[[m]] <- data.frame(model = m, f$coef, true = tb, row.names = NULL)
    }
    if (!is.null(f$selected) || m == "S4_lasso_enet") {
      sm <- selection_metrics(f$selected, sc$truth$active, all_vars)
      for (k in names(sm)) add(m, k, sm[k])
    }
    add(m, "fit_seconds", fit_sec)
    add(m, "predict_seconds", proc.time()[["elapsed"]] - t1)
    if (!is.null(f$diag)) {
      d <- unlist(f$diag)
      diags[[m]] <- data.frame(model = m, key = names(d), value = as.character(d))
    }
    if (verbose) message(sprintf("  %-20s fit %6.1fs  C=%.3f  IBS=%.4f", m, fit_sec,
                                 harrell_c(te$time, te$status, risk),
                                 ibs(te$time, te$status, S[, ig], grid)))
  }
  info <- data.frame(scenario = sc$id, rep = rep, n_events = sum(tr$status),
                     event_frac = mean(tr$status))
  list(info = info,
       metrics = cbind(scenario = sc$id, rep = rep, do.call(base::rbind, metrics)),
       coefs = if (length(coefs)) cbind(scenario = sc$id, rep = rep, do.call(base::rbind, coefs)) else NULL,
       diags = if (length(diags)) cbind(scenario = sc$id, rep = rep, do.call(base::rbind, diags)) else NULL)
}
