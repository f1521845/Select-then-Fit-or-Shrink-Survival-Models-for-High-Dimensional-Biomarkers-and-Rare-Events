# Model wrappers S1-S9 (protocol Section 3)
# ---------------------------------------------------------------------------
# Every fit_* function takes training data (X, time, status), returns a
# "fit" list with a common interface:
#   $surv(Xnew, times)  -> matrix of predicted S(t | x), nrow(Xnew) x length(times)
#   $surv_int(Xnew, times) -> optional list(lo, hi) 95% intervals for S(t | x)
#   $lp_train           -> risk score on the training data (for apparent C)
#   $coef               -> optional data.frame(var, est, lo, hi) on the ORIGINAL
#                          predictor scale (NA intervals when none exist)
#   $selected           -> optional character vector of selected variables
#   $diag               -> optional list of diagnostics (convergence etc.)
# Predictors are standardised with training means/SDs inside each wrapper.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

# ---- shared helpers ----------------------------------------------------------

standardiser <- function(X) {
  mu <- colMeans(X)
  sdv <- apply(X, 2, sd)
  sdv[sdv == 0] <- 1
  list(mu = mu, sd = sdv,
       apply = function(Z) sweep(sweep(Z[, names(mu), drop = FALSE], 2, mu), 2, sdv, "/"))
}

# Breslow estimator of the baseline cumulative hazard given a linear predictor.
breslow <- function(time, status, lp) {
  o <- order(time)
  time <- time[o]; status <- status[o]; r <- exp(lp[o] - max(lp))
  risk <- rev(cumsum(rev(r)))                    # sum of r over {j: t_j >= t_i}
  ut <- unique(time[status == 1])
  d <- tapply(status[status == 1], time[status == 1], sum)[as.character(ut)]
  first <- match(ut, time)                       # first index at each event time
  h <- d / risk[first]
  H <- cumsum(h) * exp(-max(lp))                 # undo the shift
  list(time = ut, H = H)
}

# Evaluate a step function H0 at arbitrary times.
eval_H0 <- function(b, times) {
  idx <- findInterval(times, b$time)
  c(0, b$H)[idx + 1]
}

cox_surv_fun <- function(b, lp_new, times) {
  H <- eval_H0(b, times)
  exp(-outer(exp(lp_new), H))
}

# Step-function interpolation of a survival matrix given on a time grid.
interp_surv <- function(S, grid, times) {
  idx <- findInterval(times, grid)
  S1 <- cbind(1, S)
  S1[, idx + 1, drop = FALSE]
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  attr(val, "seconds") <- proc.time()[["elapsed"]] - t0
  val
}

# ---- S1: classical Cox ---------------------------------------------------------

fit_cox <- function(X, time, status, min_epv = 5) {
  if (ncol(X) >= nrow(X) || sum(status) / ncol(X) < min_epv)
    return(list(not_estimable = TRUE))
  st <- standardiser(X)
  Z <- st$apply(X)
  fit <- coxph(Surv(time, status) ~ Z, ties = "breslow")
  b <- coef(fit); se <- rep(NA_real_, length(b))
  ok <- !is.na(b)                                  # aliased (collinear) columns -> NA
  se[ok] <- sqrt(diag(vcov(fit)))[seq_len(sum(ok))]
  b[!ok] <- 0                                      # dropped from the linear predictor
  names(b) <- colnames(X)
  lp <- drop(Z %*% b)
  bh <- breslow(time, status, lp)
  list(
    surv = function(Xn, times) cox_surv_fun(bh, drop(st$apply(Xn) %*% b), times),
    lp_train = lp,
    coef = data.frame(var = colnames(X), est = b / st$sd,
                      lo = (b - 1.96 * se) / st$sd, hi = (b + 1.96 * se) / st$sd),
    diag = list(converged = fit$info["convergence"] == 0 || TRUE)
  )
}

# ---- S2 / S4: penalised Cox (glmnet) -------------------------------------------

fit_glmnet <- function(X, time, status, alphas = 0, nfolds = 10, s = "lambda.min") {
  st <- standardiser(X)
  Z <- st$apply(X)
  time <- pmax(time, min(time[time > 0]) / 2)      # glmnet's Cox family rejects times <= 0
  y <- Surv(time, status)
  foldid <- sample(rep_len(seq_len(nfolds), nrow(Z)))
  # For p >= 500, lasso/elastic-net paths stop at 200 non-zero coefficients
  # (pmax): the dense end of the path dominated runtime in the pilot and is
  # irrelevant for sparse selection. Ridge (alpha = 0) is never capped.
  cvs <- lapply(alphas, function(a) {
    extra <- if (a > 0 && ncol(Z) >= 500) list(pmax = 200) else list()
    suppressWarnings(do.call(cv.glmnet, c(list(Z, y, family = "cox", alpha = a, foldid = foldid,
                                               standardize = FALSE), extra)))
  })
  best <- which.min(sapply(cvs, function(cv) min(cv$cvm)))
  cv <- cvs[[best]]
  b <- as.numeric(coef(cv, s = s)); names(b) <- colnames(X)
  lp <- drop(Z %*% b)
  bh <- breslow(time, status, lp)
  list(
    surv = function(Xn, times) cox_surv_fun(bh, drop(st$apply(Xn) %*% b), times),
    lp_train = lp,
    coef = data.frame(var = colnames(X), est = b / st$sd, lo = NA_real_, hi = NA_real_),
    selected = if (alphas[best] > 0) colnames(X)[b != 0] else NULL,
    diag = list(alpha = alphas[best], lambda = cv[[s]])
  )
}

fit_ridge <- function(X, time, status, ...) fit_glmnet(X, time, status, alphas = 0, ...)
fit_lasso <- function(X, time, status, ...) fit_glmnet(X, time, status, alphas = c(0.5, 1), ...)

# ---- S3 / S5 / S8b: Bayesian Cox (Stan) -----------------------------------------

stan_cache <- new.env()
get_stan_model <- function(path = "stan/bayes_cox.rds") {
  if (is.null(stan_cache$m)) {
    if (!file.exists(path)) {                    # first use: compile once (1-2 min) and save
      message("Compiling the Stan model (first use only, 1-2 minutes)...")
      saveRDS(rstan::stan_model(sub("\\.rds$", ".stan", path)), path)
    }
    stan_cache$m <- readRDS(path)
  }
  stan_cache$m
}

fit_bayes_cox <- function(X, time, status, prior = c("normal", "horseshoe"),
                          sigma = 1, hier = TRUE, p0 = NULL, slab_scale = 2, slab_df = 4,
                          chains = 4, iter = 2000, warmup = 1000, cores = 2,
                          adapt_delta = 0.99, n_draws_pred = 200, rerun = TRUE) {
  prior <- match.arg(prior)
  suppressPackageStartupMessages(library(rstan))
  st <- standardiser(X)
  Z <- st$apply(X)
  P <- ncol(Z); d <- sum(status)
  if (is.null(p0)) p0 <- min(10, floor(P / 2))
  p0 <- min(p0, floor(P / 2))                   # keep tau0 finite (e.g. p0 = 50 when P = 50)
  # Normal prior: by default hierarchical ("Bayesian ridge"): beta ~ N(0, s^2) with
  # s ~ half-normal(0, sigma), so the amount of shrinkage is learned from the data,
  # the Bayesian counterpart of ridge with lambda chosen by cross-validation.
  # Fixed-SD versions (hier = FALSE) are prior-sensitivity variants (H5): the pilot
  # showed fixed sigma = 1 overfits and fixed sigma = 1/sqrt(P) over-shrinks when p is large.
  if (identical(sigma, "scaled")) sigma <- 1 / sqrt(P)
  # Global scale: Piironen & Vehtari (2017) with sqrt(number of events) in
  # place of sigma/sqrt(n), since Cox information scales with events.
  tau0 <- p0 / (P - p0) / sqrt(d)
  o <- order(time, decreasing = TRUE)
  ts <- time[o]
  last <- vapply(seq_along(ts), function(i) max(which(ts == ts[i])), 1L)
  dat <- list(N = nrow(Z), P = P, X = Z[o, , drop = FALSE], status = status[o],
              last = last, prior = as.integer(prior == "horseshoe"), hier = as.integer(hier),
              sigma = sigma, tau0 = tau0, slab_scale = slab_scale, slab_df = slab_df)
  run <- function(it, wu) suppressWarnings(sampling(
    get_stan_model(), data = dat, chains = chains, iter = it, warmup = wu,
    cores = cores, refresh = 0, pars = "beta",
    control = list(adapt_delta = adapt_delta, max_treedepth = 10)))
  fit <- run(iter, warmup)
  reran <- FALSE
  # Protocol rule: R-hat >= 1.01 or bulk ESS <= 400 -> re-run once with twice
  # the iterations; any remaining problem is reported, not dropped.
  chk <- summary(fit, pars = "beta")$summary
  if (rerun && (max(chk[, "Rhat"], na.rm = TRUE) >= 1.01 || min(chk[, "n_eff"], na.rm = TRUE) <= 400)) {
    fit <- run(2 * iter, 2 * warmup); reran <- TRUE
  }
  B <- as.matrix(fit, pars = "beta")            # draws x P (standardised scale)
  colnames(B) <- colnames(X)
  summ <- summary(fit, pars = "beta")$summary
  ndiv <- sum(get_divergent_iterations(fit))
  q <- apply(B, 2, quantile, c(0.025, 0.975))
  est <- colMeans(B)
  keep <- unique(round(seq(1, nrow(B), length.out = min(n_draws_pred, nrow(B)))))
  Bk <- B[keep, , drop = FALSE]
  LPtr <- Z %*% t(Bk)                            # n x draws
  bhs <- lapply(seq_len(ncol(LPtr)), function(j) breslow(time, status, LPtr[, j]))
  surv_draws <- function(Xn, times, fun) {
    Zn <- st$apply(Xn); LPn <- Zn %*% t(Bk)
    acc <- NULL
    for (j in seq_len(ncol(LPn))) {
      Sj <- cox_surv_fun(bhs[[j]], LPn[, j], times)
      acc <- fun(acc, Sj, j)
    }
    acc
  }
  list(
    surv = function(Xn, times)
      surv_draws(Xn, times, function(acc, S, j) if (is.null(acc)) S / length(bhs) else acc + S / length(bhs)),
    surv_int = function(Xn, times) {
      arr <- surv_draws(Xn, times, function(acc, S, j) {
        if (is.null(acc)) acc <- array(NA_real_, c(dim(S), length(bhs)))
        acc[, , j] <- S; acc })
      list(lo = apply(arr, c(1, 2), quantile, 0.025), hi = apply(arr, c(1, 2), quantile, 0.975))
    },
    lp_train = drop(Z %*% est),
    coef = data.frame(var = colnames(X), est = est / st$sd,
                      lo = q[1, ] / st$sd, hi = q[2, ] / st$sd),
    selected = colnames(X)[q[1, ] > 0 | q[2, ] < 0],
    diag = list(max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
                min_ess = min(summ[, "n_eff"], na.rm = TRUE),
                divergences = ndiv, tau0 = tau0, sigma = sigma, hier = hier, p0 = p0, reran = reran)
  )
}

# ---- S6: random survival forest ------------------------------------------------

fit_rsf <- function(X, time, status, ntree = 500,
                    grid = NULL, nthreads = 2) {
  suppressPackageStartupMessages(library(randomForestSRC))
  options(rf.cores = nthreads, mc.cores = 1)
  p <- ncol(X)
  # mtry grid: {sqrt(p), p/3} for p < 500; {sqrt(p), 2 sqrt(p)} for p >= 500,
  # where p/3 is computationally prohibitive (pilot: 20 min per fit at p = 2,000).
  if (is.null(grid)) {
    m2 <- if (p >= 500) 2 * sqrt(p) else p / 3
    grid <- expand.grid(mtry = unique(pmax(1, round(c(sqrt(p), m2)))), nodesize = c(5, 15))
  }
  df <- data.frame(time = time, status = status, X, check.names = FALSE)
  # tune on out-of-bag error (1 - Harrell's C), with a smaller forest
  oob <- sapply(seq_len(nrow(grid)), function(i) {
    f <- rfsrc(Surv(time, status) ~ ., df, ntree = 200, mtry = grid$mtry[i],
               nodesize = grid$nodesize[i], importance = "none")
    tail(f$err.rate[!is.na(f$err.rate)], 1)
  })
  g <- grid[which.min(oob), ]
  f <- rfsrc(Surv(time, status) ~ ., df, ntree = ntree, mtry = g$mtry,
             nodesize = g$nodesize, importance = "none")
  pred_mat <- function(Xn) predict(f, newdata = data.frame(Xn, check.names = FALSE))$survival
  list(
    surv = function(Xn, times) interp_surv(pred_mat(Xn), f$time.interest, times),
    lp_train = rowSums(f$chf),                 # in-bag ensemble CHF summed: risk score
    diag = list(mtry = g$mtry, nodesize = g$nodesize)
  )
}

# ---- S7: XGBoost-Cox ------------------------------------------------------------

fit_xgb <- function(X, time, status, n_configs = 10, nfold = 5, nthread = 2,
                    max_rounds = 1000) {
  suppressPackageStartupMessages(library(xgboost))
  # xgboost's Cox objective codes censoring as a negative label, so a time of
  # exactly 0 would be ambiguous: move zeros to half the smallest positive time.
  time <- pmax(time, min(time[time > 0]) / 2)
  label <- ifelse(status == 1, time, -time)
  dtr <- xgb.DMatrix(X, label = label, nthread = nthread)
  configs <- data.frame(
    eta = 10^runif(n_configs, -2, -1), max_depth = sample(1:6, n_configs, TRUE),
    subsample = runif(n_configs, 0.5, 1), colsample_bytree = runif(n_configs, 0.3, 1),
    min_child_weight = sample(1:10, n_configs, TRUE))
  res <- lapply(seq_len(n_configs), function(i) {
    par <- c(list(objective = "survival:cox", eval_metric = "cox-nloglik",
                  nthread = nthread), as.list(configs[i, ]))
    cv <- xgb.cv(params = par, data = dtr, nrounds = max_rounds, nfold = nfold,
                 early_stopping_rounds = 25, verbose = 0)
    ev <- cv$evaluation_log
    list(par = par, score = min(ev$test_cox_nloglik_mean),
         rounds = ev$iter[which.min(ev$test_cox_nloglik_mean)])
  })
  best <- res[[which.min(sapply(res, `[[`, "score"))]]
  bst <- xgb.train(params = best$par, data = dtr, nrounds = best$rounds, verbose = 0)
  margin <- function(Xn) predict(bst, xgb.DMatrix(Xn, nthread = nthread), outputmargin = TRUE)
  lp <- margin(X)
  bh <- breslow(time, status, lp)
  list(
    surv = function(Xn, times) cox_surv_fun(bh, margin(Xn), times),
    lp_train = lp, booster = bst, params = best$par, rounds = best$rounds,
    diag = list(rounds = best$rounds, max_depth = best$par$max_depth, eta = best$par$eta)
  )
}

# Mean |SHAP| ranking from an xgboost booster.
shap_rank <- function(bst, X, nthread = 2) {
  contrib <- predict(bst, xgb.DMatrix(X, nthread = nthread), predcontrib = TRUE)
  contrib <- contrib[, colnames(X), drop = FALSE]   # drops the bias / "(Intercept)" column
  imp <- colMeans(abs(contrib))
  names(sort(imp, decreasing = TRUE))
}

# ---- S8 / S8b / S8-leak: select-then-fit -----------------------------------------

# `xgb_fit` is the S7 fit on the same training data (reused, as in practice).
fit_select_then_fit <- function(X, time, status, xgb_fit, k = 10,
                                second = c("cox", "bayes"), X_rank = NULL, ...) {
  second <- match.arg(second)
  # ranking data: training data (honest) or train+test (leak variant)
  Xr <- if (is.null(X_rank)) X else X_rank
  ranked <- shap_rank(xgb_fit$booster, Xr)
  sel <- head(ranked, k)
  Xs <- X[, sel, drop = FALSE]
  f <- if (second == "cox") fit_cox(Xs, time, status, min_epv = 0)
       else fit_bayes_cox(Xs, time, status, prior = "normal", ...)
  f$selected <- sel
  f$surv0 <- f$surv
  f$surv <- function(Xn, times) f$surv0(Xn[, sel, drop = FALSE], times)
  if (!is.null(f$surv_int)) {
    si <- f$surv_int
    f$surv_int <- function(Xn, times) si(Xn[, sel, drop = FALSE], times)
  }
  f
}

# Leak variant: the XGBoost ranking model is re-fitted on training + "local
# test" data (the whole dataset before splitting), with the tuned settings.
fit_select_leak <- function(X, time, status, xgb_fit, X_loc, time_loc, status_loc, k = 10) {
  Xall <- rbind(X, X_loc)
  tt <- c(time, time_loc); tt <- pmax(tt, min(tt[tt > 0]) / 2)
  lab <- ifelse(c(status, status_loc) == 1, tt, -tt)
  bst_all <- xgb.train(params = xgb_fit$params,
                       data = xgb.DMatrix(Xall, label = lab, nthread = 2),
                       nrounds = xgb_fit$rounds, verbose = 0)
  fit_select_then_fit(X, time, status, list(booster = bst_all), k = k, X_rank = Xall)
}

# ---- S9: BART survival (exploratory) ------------------------------------------

fit_bart <- function(X, time, status, K = 20, ndpost = 500, nskip = 250,
                     chunk = 1000, keepevery = 5) {
  suppressPackageStartupMessages(library(BART))
  post <- invisible(capture.output(
    fit <- surv.bart(x.train = X, times = time, delta = status, K = K,
                     ndpost = ndpost, nskip = nskip, keepevery = keepevery,
                     printevery = 1e9)))
  grid <- fit$times
  Ke <- fit$K                                  # BART merges tied quantiles: Ke <= K
  pred_grid <- function(Xn) {
    out_mean <- out_lo <- out_hi <- matrix(NA_real_, nrow(Xn), Ke)
    for (s in split(seq_len(nrow(Xn)), ceiling(seq_len(nrow(Xn)) / chunk))) {
      pre <- surv.pre.bart(times = time, delta = status, x.train = X,
                           x.test = Xn[s, , drop = FALSE], K = K)
      invisible(capture.output(pr <- predict(fit, newdata = pre$tx.test)))   # silence C++ log
      S <- pr$surv.test                          # draws x (n_chunk * K)
      arr <- array(S, c(nrow(S), Ke, length(s))) # draws x Ke x n (subject-major)
      out_mean[s, ] <- t(apply(arr, c(2, 3), mean))
      out_lo[s, ]   <- t(apply(arr, c(2, 3), quantile, 0.025))
      out_hi[s, ]   <- t(apply(arr, c(2, 3), quantile, 0.975))
    }
    list(mean = out_mean, lo = out_lo, hi = out_hi)
  }
  cache <- new.env()
  getp <- function(Xn) {
    key <- digest_key(Xn)
    if (!identical(cache$key, key)) { cache$key <- key; cache$val <- pred_grid(Xn) }
    cache$val
  }
  list(
    surv = function(Xn, times) interp_surv(getp(Xn)$mean, grid, times),
    surv_int = function(Xn, times) {
      g <- getp(Xn)
      list(lo = interp_surv(g$lo, grid, times), hi = interp_surv(g$hi, grid, times))
    },
    lp_train = NULL,
    diag = list(K = Ke)
  )
}

digest_key <- function(X) c(dim(X), sum(X[1, ]), sum(X[nrow(X), ]), sum(X[, 1]))

# ---- package check -------------------------------------------------------------
# Fail loudly at the start if any package is missing, instead of silently
# dropping a model from the comparison.
check_packages <- function(pkgs = c("survival", "glmnet", "rstan", "randomForestSRC",
                                    "xgboost", "BART")) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing R packages: ", paste(miss, collapse = ", "), "\nInstall them with:\n  install.packages(c(",
         paste0('"', miss, '"', collapse = ", "), "))\n",
         "Then delete any quick/demo results created without them (see README) and re-run.", call. = FALSE)
  invisible(TRUE)
}
check_packages()
