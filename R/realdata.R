# Real-data engine (protocol Section 5)
# ---------------------------------------------------------------------------
# Repeated K-fold cross-validation with all data-dependent steps inside each
# training fold (imputation, screening, standardisation, tuning), plus external
# validation (train on one dataset, test on another). Models are the same
# fit_* wrappers used in the simulation (R/models.R).
#
# Performance is measured against OBSERVED (censored) outcomes:
#   Harrell's C, Uno's C (truncated at t*75), IPCW integrated Brier score up to
#   t*75, cumulative/dynamic time-dependent AUC at t*50 (Uno et al. 2007),
#   calibration slope at t*50 (Cox regression of the outcome on
#   cloglog(1 - S_pred)) and ICI at t*50 (Austin, Harrell & van Klaveren 2020,
#   using a Cox model with a penalised spline of cloglog(1 - S_pred)).

source("R/models.R"); source("R/metrics.R")
check_packages(c("mice", "ggplot2"))
suppressPackageStartupMessages(library(survival))
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- metrics against observed outcomes ----------------------------------------

td_auc <- function(time, status, risk, t) {
  km_c <- survfit(Surv(time, 1 - status) ~ 1)
  G <- stepfun(km_c$time, c(1, km_c$surv))
  cases <- which(time <= t & status == 1); ctrls <- which(time > t)
  if (!length(cases) || !length(ctrls)) return(NA_real_)
  w <- 1 / pmax(G(time[cases] - 1e-10), 1e-8)
  rc <- sort(risk[ctrls])
  below <- findInterval(risk[cases], rc, left.open = TRUE)    # controls with lower risk
  ties <- findInterval(risk[cases], rc) - below
  sum(w * (below + 0.5 * ties)) / (sum(w) * length(ctrls))
}

calib_observed <- function(time, status, S_t, t) {
  F <- pmin(pmax(1 - S_t, 1e-6), 1 - 1e-6)
  cll <- log(-log(1 - F))
  slope <- tryCatch(unname(coef(coxph(Surv(time, status) ~ cll))), error = function(e) NA_real_)
  ici <- tryCatch({
    fit <- coxph(Surv(time, status) ~ pspline(cll, df = 3))
    bh <- basehaz(fit, centered = FALSE)
    H0 <- approx(bh$time, bh$hazard, xout = t, method = "constant", rule = 2, f = 0)$y
    obs <- 1 - exp(-H0 * exp(predict(fit, type = "lp", reference = "zero")))
    mean(abs(F - obs))
  }, error = function(e) NA_real_)
  c(calib_slope = slope, ici = ici)
}

observed_metrics <- function(time, status, S, times, tstar, grid) {
  it <- match(tstar, times); ig <- match(grid, times)
  risk <- 1 - S[, it[2]]
  cal <- calib_observed(time, status, S[, it[2]], tstar[2])
  c(harrell_c = harrell_c(time, status, risk),
    uno_c = uno_c(time, status, risk, tstar[3]),
    td_auc = td_auc(time, status, risk, tstar[2]),
    ibs = ibs(time, status, S[, ig, drop = FALSE], grid),
    cal)
}

# ---- preprocessing inside a fold ------------------------------------------------

# Fast univariate Cox score statistics for every column (screening).
cox_score_screen <- function(X, time, status, K) {
  if (ncol(X) <= K) return(colnames(X))
  o <- order(time, decreasing = TRUE)
  Xs <- scale(X[o, , drop = FALSE]); d <- status[o]
  Xs[is.na(Xs)] <- 0
  n_risk <- seq_along(o)                               # risk set = 1..i (descending time)
  S1 <- apply(Xs, 2, cumsum) / n_risk                  # risk-set means
  S2 <- apply(Xs^2, 2, cumsum) / n_risk
  ev <- d == 1
  U <- colSums(Xs[ev, , drop = FALSE] - S1[ev, , drop = FALSE])
  V <- colSums(S2[ev, , drop = FALSE] - S1[ev, , drop = FALSE]^2)
  score <- U^2 / pmax(V, 1e-12)
  colnames(X)[order(score, decreasing = TRUE)[seq_len(K)]]
}

# Impute predictors with mice, fitted on the training rows only (test rows are
# imputed from the training-row models via `ignore`). The outcome is NOT used,
# so test outcomes can never leak into test covariates.
impute_fold <- function(df_tr, df_te, m = 5, seed = 1) {
  if (!anyNA(df_tr) && !anyNA(df_te)) return(rep(list(list(tr = df_tr, te = df_te)), 1))
  if (ncol(df_tr) > 200) {
    # high-dimensional data: simple training-fold median / mode imputation of the
    # (few) incomplete columns; mice with thousands of columns is impractical
    fill <- lapply(df_tr, function(v) if (is.numeric(v)) median(v, na.rm = TRUE)
                   else names(which.max(table(v))))
    for (v in names(df_tr)) {
      if (anyNA(df_tr[[v]])) df_tr[[v]][is.na(df_tr[[v]])] <- fill[[v]]
      if (anyNA(df_te[[v]])) df_te[[v]][is.na(df_te[[v]])] <- fill[[v]]
    }
    return(list(list(tr = df_tr, te = df_te)))
  }
  all <- rbind(df_tr, df_te)
  ign <- c(rep(FALSE, nrow(df_tr)), rep(TRUE, nrow(df_te)))
  imp <- suppressWarnings(mice::mice(all, m = m, maxit = 5, ignore = ign,   # pmm / logreg / polyreg defaults
                               printFlag = FALSE, seed = seed))
  lapply(seq_len(m), function(k) {
    cmp <- mice::complete(imp, k)   # mice is not attached: attaching it masks base::rbind
    list(tr = cmp[!ign, , drop = FALSE], te = cmp[ign, , drop = FALSE])
  })
}

# Factors -> treatment-coded dummies, with levels fixed on the full dataset.
encode <- function(df, levels_ref) {
  for (v in names(levels_ref)) df[[v]] <- factor(df[[v]], levels = levels_ref[[v]])
  X <- model.matrix(~ ., df)[, -1, drop = FALSE]
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  X
}

# ---- one train/test evaluation ----------------------------------------------------

#' Fit `models` on a training set and evaluate them on a test set.
#' @param df_tr,df_te data.frames of predictors (may contain NA and factors)
#' @param spec list(m_imp, screen_K, tstar, grid, stan_args, xgb_configs, bart_args, levels)
eval_split <- function(df_tr, y_tr, df_te, y_te, models, spec, seed = 1) {
  if (any(grepl("^S8", models))) models <- union("S7_xgb", models)   # S8 reuses the S7 fit
  models <- intersect(names(MODEL_LABELS_RD), models)                  # fixed order
  sets <- impute_fold(df_tr, df_te, m = spec$m_imp, seed = seed)
  times <- sort(unique(c(spec$tstar, spec$grid)))
  preds <- list(); app <- list(); sel <- list(); secs <- list(); fails <- character(0)
  for (k in seq_along(sets)) {
    Xtr <- encode(sets[[k]]$tr, spec$levels); Xte <- encode(sets[[k]]$te, spec$levels)
    keep <- apply(Xtr, 2, sd) > 0
    Xtr <- Xtr[, keep, drop = FALSE]; Xte <- Xte[, keep, drop = FALSE]
    if (!is.null(spec$screen_K) && ncol(Xtr) > spec$screen_K) {
      # always keep clinical (non-gene) columns; screen the rest
      clin <- grep(spec$protect_regex %||% "^$", colnames(Xtr), value = TRUE)
      cand <- setdiff(colnames(Xtr), clin)
      ks <- cox_score_screen(Xtr[, cand, drop = FALSE], y_tr$time, y_tr$status, spec$screen_K)
      Xtr <- Xtr[, c(clin, ks), drop = FALSE]; Xte <- Xte[, c(clin, ks), drop = FALSE]
    }
    fits <- list()
    for (m in models) {
      t0 <- proc.time()[["elapsed"]]
      f <- tryCatch(switch(m,
        S1_cox = fit_cox(Xtr, y_tr$time, y_tr$status),
        S2_ridge = fit_ridge(Xtr, y_tr$time, y_tr$status),
        S3_bayes_normal = do.call(fit_bayes_cox, c(list(Xtr, y_tr$time, y_tr$status, prior = "normal"), spec$stan_args)),
        S4_lasso_enet = fit_lasso(Xtr, y_tr$time, y_tr$status),
        S5_bayes_horseshoe = do.call(fit_bayes_cox, c(list(Xtr, y_tr$time, y_tr$status, prior = "horseshoe"), spec$stan_args)),
        S6_rsf = fit_rsf(Xtr, y_tr$time, y_tr$status),
        S7_xgb = fit_xgb(Xtr, y_tr$time, y_tr$status, n_configs = spec$xgb_configs),
        S8_select_fit = fit_select_then_fit(Xtr, y_tr$time, y_tr$status, fits$S7_xgb),
        S8b_select_bayes = do.call(fit_select_then_fit, c(list(Xtr, y_tr$time, y_tr$status, fits$S7_xgb, second = "bayes"), spec$stan_args)),
        S8leak_select_fit = fit_select_leak(Xtr, y_tr$time, y_tr$status, fits$S7_xgb, Xte, y_te$time, y_te$status),
        S9_bart = do.call(fit_bart, c(list(Xtr, y_tr$time, y_tr$status), spec$bart_args)),
        stop("unknown model ", m)),
        error = function(e) list(error = conditionMessage(e)))
      fits[[m]] <- f
      secs[[m]] <- (secs[[m]] %||% 0) + proc.time()[["elapsed"]] - t0
      if (!is.null(f$error) || !is.null(f$not_estimable)) {
        fails <- c(fails, paste0(m, ": ", f$error %||% "not estimable")); next
      }
      ev <- tryCatch({
        S <- f$surv(Xte, times)
        if (anyNA(S)) stop("NA predictions")
        preds[[m]] <- if (is.null(preds[[m]])) S / length(sets) else preds[[m]] + S / length(sets)
        if (!is.null(f$lp_train) && k == 1) app[[m]] <- harrell_c(y_tr$time, y_tr$status, f$lp_train)
        if (!is.null(f$selected) && k == 1) sel[[m]] <- f$selected
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(ev)) fails <- c(fails, paste0(m, ": prediction failed: ", ev))
    }
  }
  rows <- list()
  for (m in names(preds)) {
    v <- observed_metrics(y_te$time, y_te$status, preds[[m]], times, spec$tstar, spec$grid)
    if (!is.null(app[[m]])) v <- c(v, apparent_c = app[[m]])
    v <- c(v, fit_seconds = secs[[m]])
    rows[[m]] <- data.frame(model = m, metric = names(v), value = unname(v))
  }
  list(metrics = do.call(base::rbind, rows), selected = sel, failures = unique(fails), preds = preds)
}

# ---- repeated K-fold cross-validation -------------------------------------------------

make_folds <- function(status, K = 5, repeats = 10, seed = 2026) {
  set.seed(seed)
  lapply(seq_len(repeats), function(r) {
    f <- integer(length(status))
    for (s in 0:1) {                                  # stratified by event status
      idx <- which(status == s)
      f[idx] <- sample(rep_len(seq_len(K), length(idx)))
    }
    f
  })
}

#' @param df predictors (data.frame), @param y data.frame(time, status)
#' @param bayes_folds Bayesian arms run only in the first `bayes_folds` outer folds
run_cv <- function(df, y, models, spec, K = 5, repeats = 10, bayes_folds = Inf,
                   cache_dir = NULL, verbose = TRUE) {
  folds <- make_folds(y$status, K, repeats)
  out <- list(); sels <- list(); fails <- list(); j <- 0
  for (r in seq_len(repeats)) for (k in seq_len(K)) {
    j <- j + 1
    f_cache <- if (!is.null(cache_dir)) file.path(cache_dir, sprintf("fold_r%02d_k%d.rds", r, k))
    if (!is.null(f_cache) && file.exists(f_cache)) {
      res <- readRDS(f_cache)
    } else {
      te <- folds[[r]] == k
      ms <- if (j > bayes_folds) setdiff(models, BAYES_MODELS_RD) else models
      if (verbose) message(sprintf("[%s] repeat %d fold %d (%d models)", format(Sys.time(), "%H:%M:%S"), r, k, length(ms)))
      set.seed(1000 * r + k)
      res <- eval_split(df[!te, , drop = FALSE], y[!te, ], df[te, , drop = FALSE], y[te, ],
                        ms, spec, seed = 1000 * r + k)
      res$preds <- NULL
      if (!is.null(f_cache)) { dir.create(cache_dir, FALSE, TRUE); saveRDS(res, f_cache) }
    }
    out[[j]] <- cbind(repeat_id = r, fold = k, res$metrics)
    sels[[j]] <- res$selected; fails[[j]] <- res$failures
  }
  list(metrics = do.call(base::rbind, out), selected = sels, failures = unlist(fails))
}

BAYES_MODELS_RD <- c("S3_bayes_normal", "S5_bayes_horseshoe", "S8b_select_bayes", "S9_bart")

# ---- summaries -------------------------------------------------------------------------

MODEL_LABELS_RD <- c(
  S1_cox = "Cox", S2_ridge = "Ridge Cox", S3_bayes_normal = "Bayesian Cox (normal)",
  S4_lasso_enet = "Lasso/EN Cox", S5_bayes_horseshoe = "Bayesian Cox (horseshoe)",
  S6_rsf = "RSF", S7_xgb = "XGBoost-Cox", S8_select_fit = "Select-then-fit",
  S8b_select_bayes = "Select-then-fit (Bayes)", S8leak_select_fit = "Select-then-fit (leak)",
  S9_bart = "BART")

summarise_cv <- function(M, metrics = c("ibs", "harrell_c", "uno_c", "td_auc", "calib_slope",
                                         "ici", "apparent_c", "fit_seconds")) {
  M <- M[M$metric %in% metrics, ]
  a <- aggregate(value ~ model + metric, M, function(x) c(mean = mean(x, na.rm = TRUE),
                                                          sd = sd(x, na.rm = TRUE), folds = sum(!is.na(x))))
  a <- cbind(a[, 1:2], as.data.frame(a$value))
  w <- reshape(a[, c("model", "metric", "mean", "sd")], idvar = "model", timevar = "metric",
               direction = "wide")
  w$model <- factor(w$model, names(MODEL_LABELS_RD))
  w <- w[order(w$model), ]
  w$label <- MODEL_LABELS_RD[as.character(w$model)]
  w
}

# Pretty table: "mean (SD)" for each metric
format_cv_table <- function(w, metrics = c("ibs", "harrell_c", "uno_c", "td_auc", "calib_slope", "ici")) {
  out <- data.frame(Model = w$label)
  for (m in metrics) {
    mu <- w[[paste0("mean.", m)]]; s <- w[[paste0("sd.", m)]]
    if (!is.null(mu)) out[[m]] <- ifelse(is.na(mu), "-", sprintf("%.3f (%.3f)", mu, s))
  }
  names(out) <- c("Model", c(ibs = "IBS", harrell_c = "Harrell C", uno_c = "Uno C", td_auc = "AUC(t*50)",
                             calib_slope = "Calib. slope", ici = "ICI")[metrics[metrics %in% names(out)]])
  out
}

# Selection stability across folds: selection frequency and mean pairwise Jaccard
selection_stability <- function(sels, model) {
  s <- lapply(sels, `[[`, model); s <- s[!vapply(s, is.null, TRUE)]
  if (length(s) < 2) return(NULL)
  freq <- sort(table(unlist(s)) / length(s), decreasing = TRUE)
  pairs <- combn(length(s), 2)
  jac <- apply(pairs, 2, function(p) {
    a <- s[[p[1]]]; b <- s[[p[2]]]
    if (!length(union(a, b))) NA_real_ else length(intersect(a, b)) / length(union(a, b))
  })
  # Jaccard is undefined for two empty sets (e.g. a horseshoe that selects nothing)
  list(freq = freq, jaccard = if (all(is.na(jac))) NA_real_ else mean(jac, na.rm = TRUE),
       mean_size = mean(lengths(s)), folds = length(s))
}

# Guide lookup: what does the simulation-based practitioner guide recommend for
# a dataset with this n, event fraction and p? Picks the nearest guide cell
# (log scale for n and p, linear for event fraction).
guide_lookup <- function(guide, n, event_frac, p, linear = TRUE) {
  g <- guide[guide$linear == linear, ]
  if (!nrow(g)) return(NULL)
  d <- (log(g$n) - log(n))^2 / log(10)^2 + (log(g$p) - log(p))^2 / log(10)^2 +
       ((g$event_frac - event_frac) / 0.3)^2
  g[which.min(d), ]
}

# ---- full-data fits (associations) and plots ---------------------------------------------

#' Fit Cox-type models on the full dataset (single imputation) and return
#' hazard ratios per SD with 95% intervals where the method provides them.
fit_full_coefs <- function(df, y, spec, models = c("S1_cox", "S2_ridge", "S3_bayes_normal",
                                                    "S4_lasso_enet", "S5_bayes_horseshoe")) {
  set.seed(7)
  imp <- impute_fold(df, df[0, , drop = FALSE], m = 1)[[1]]$tr
  X <- encode(imp, spec$levels); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  sdx <- apply(X, 2, sd)
  out <- list()
  for (m in models) {
    f <- tryCatch(switch(m,
      S1_cox = fit_cox(X, y$time, y$status, min_epv = 0),
      S2_ridge = fit_ridge(X, y$time, y$status),
      S3_bayes_normal = do.call(fit_bayes_cox, c(list(X, y$time, y$status, prior = "normal"), spec$stan_args)),
      S4_lasso_enet = fit_lasso(X, y$time, y$status),
      S5_bayes_horseshoe = do.call(fit_bayes_cox, c(list(X, y$time, y$status, prior = "horseshoe"), spec$stan_args))),
      error = function(e) NULL)
    if (is.null(f) || is.null(f$coef)) next
    cf <- f$coef
    # per-SD hazard ratios (coefficients are on the original scale)
    out[[m]] <- data.frame(model = m, var = cf$var,
                           hr = exp(cf$est * sdx[cf$var]),
                           lo = exp(cf$lo * sdx[cf$var]), hi = exp(cf$hi * sdx[cf$var]))
  }
  do.call(base::rbind, out)
}

plot_cv_results <- function(M, title = NULL) {
  suppressPackageStartupMessages(library(ggplot2))
  d <- M[M$metric %in% c("ibs", "harrell_c", "calib_slope"), ]
  a <- aggregate(value ~ model + metric, d, function(x) c(m = mean(x, na.rm = TRUE), s = sd(x, na.rm = TRUE)))
  a <- cbind(a[, 1:2], as.data.frame(a$value))
  a$label <- factor(MODEL_LABELS_RD[a$model], rev(MODEL_LABELS_RD))
  a$metric <- factor(a$metric, c("ibs", "harrell_c", "calib_slope"),
                     c("Integrated Brier score (lower = better)", "Harrell's C (higher = better)",
                       "Calibration slope (1 = ideal)"))
  ref <- data.frame(metric = factor(levels(a$metric)[3], levels(a$metric)), x = 1)
  ggplot(a, aes(m, label)) +
    geom_vline(data = ref, aes(xintercept = x), colour = "grey60") +
    geom_errorbar(aes(xmin = m - s, xmax = m + s), width = 0, orientation = "y", colour = "grey45") +
    geom_point(size = 2.3, colour = "#2F5D8A") +
    facet_wrap(~metric, scales = "free_x") +
    scale_x_continuous(breaks = function(l) pretty(l, n = 4)) +
    labs(x = "Mean across outer folds (bars: ± 1 SD)", y = NULL, title = title) +
    theme_minimal(base_size = 11) + theme(panel.spacing = grid::unit(2, "lines"))
}

plot_coefs <- function(cf, title = NULL) {
  suppressPackageStartupMessages(library(ggplot2))
  cf$label <- MODEL_LABELS_RD[cf$model]
  ggplot(cf, aes(hr, var, colour = label)) +
    geom_vline(xintercept = 1, colour = "grey60") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", position = position_dodge(width = 0.7), na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.7), size = 1.8) +
    scale_x_log10() +
    labs(x = "Hazard ratio per SD (log scale; bars: 95% interval where available)", y = NULL,
         colour = NULL, title = title) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
}

# Read the practitioner guide written by 01_simulation.R (collect > demo).
read_guide <- function() {
  for (m in c("collect", "local", "demo")) {
    f <- file.path("results", m, "paper", "practitioner_guide.csv")
    if (file.exists(f)) {
      g <- read.csv(f, stringsAsFactors = FALSE)
      for (v in c("recommended", "also_fine", "avoid", "reliable_intervals", "best_mean"))
        if (v %in% names(g)) g[[v]] <- ifelse(is.na(g[[v]]), "", as.character(g[[v]]))
      return(structure(g, source = m))
    }
  }
  NULL
}

guide_check <- function(guide, w, n, event_frac, p) {
  if (is.null(guide)) return(cat("No practitioner guide found yet: run 01_simulation.R first.\n"))
  cell <- guide_lookup(guide, n, event_frac, p)
  ok <- as.character(w$model) != "S8leak_select_fit"   # the leak variant is illustrative, never a winner
  best <- w$label[ok][which.min(w$mean.ibs[ok])]
  cat(sprintf("Guide source: %s run%s\n", attr(guide, "source"),
              if (attr(guide, "source") == "demo") " (DEMO: not for interpretation)" else ""))
  cat(sprintf("Nearest guide cell: n = %s, events = %g%%, p = %s (EPV %s)\n",
              format(cell$n, big.mark = ","), 100 * cell$event_frac, format(cell$p, big.mark = ","), cell$epv))
  cat("  Guide recommends: ", cell$recommended, "\n", sep = "")
  if (nzchar(cell$also_fine)) cat("  Also fine:        ", cell$also_fine, "\n", sep = "")
  cat("  Best here (IBS):  ", best, "\n", sep = "")
  rec <- trimws(strsplit(paste(cell$recommended, cell$also_fine, sep = ", "), ",")[[1]])
  cat(if (best %in% rec) "  -> Consistent with the guide.\n" else "  -> Not what the guide predicts: discuss.\n")
  invisible(list(cell = cell, best = best))
}
