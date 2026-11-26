# Data-generating mechanisms (protocol Section 4, "D")
# ---------------------------------------------------------------------------
# Weibull proportional hazards (Bender, Augustin & Blettner 2005) with
# block-correlated normal predictors, four effect structures, and censoring
# calibrated so that a target fraction of subjects has an observed event.
#
# All "truth" needed for evaluation is returned: true coefficients, true
# active set, and a function giving the true S(t | x) for any new data.

# ---- predictors ------------------------------------------------------------

# Block-exchangeable correlation: blocks of `block` variables with
# within-block correlation rho, independent between blocks.
gen_x <- function(n, p, rho = 0.3, block = 20) {
  X <- matrix(0, n, p)
  starts <- seq(1, p, by = block)
  for (s in starts) {
    idx <- s:min(s + block - 1, p)
    k <- length(idx)
    # exchangeable: x = sqrt(rho) * z0 + sqrt(1 - rho) * z
    z0 <- rnorm(n)
    X[, idx] <- sqrt(rho) * z0 + sqrt(1 - rho) * matrix(rnorm(n * k), n, k)
  }
  colnames(X) <- sprintf("x%04d", seq_len(p))
  X
}

# Clinical covariates for the biobank-like block: age (standardised),
# sex (0/1) and one lab value (standardised), placed in front of X.
gen_clinical <- function(n) {
  cbind(c_age = rnorm(n), c_sex = rbinom(n, 1, 0.5), c_lab = rnorm(n))
}

# ---- effect structures -----------------------------------------------------

# Positions of the true predictors: one per block (first variable of each
# block) so that true signals are not all mutually correlated.
true_positions <- function(p, s, block = 20) {
  pos <- seq(1, p, by = block)
  if (length(pos) < s) pos <- seq_len(p)  # small p: take the first s
  pos[seq_len(s)]
}

# Fixed (per-scenario) coefficients. `seed` keeps the dense effects identical
# across replications, so every replication estimates the same estimand.
make_beta <- function(p, dgm, s = 10, seed = 1L, n_clin = 0) {
  beta <- numeric(p)
  s <- min(s, floor(p / 2))  # p = 10 -> 5 true predictors (protocol note)
  pos <- true_positions(p, s)
  if (dgm %in% c("sparse", "nonlinear", "nonph")) {
    hr <- seq(1.2, 2.0, length.out = s)
    sgn <- rep(c(1, -1), length.out = s)
    beta[pos] <- sgn * log(hr)
  } else if (dgm == "dense") {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    beta <- rnorm(p, 0, 0.1)
    pos <- seq_len(p)
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  } else stop("unknown dgm: ", dgm)
  names(beta) <- sprintf("x%04d", seq_len(p))
  clin <- if (n_clin > 0) c(c_age = log(1.8), c_sex = log(1.3), c_lab = -log(1.5)) else NULL
  list(beta = c(clin, beta), active = c(names(clin), names(beta)[pos]), pos = pos)
}

# Extra (non-linear / interaction) part of the log hazard, using true predictors.
nonlinear_part <- function(X, pos) {
  a <- pos[1]; b <- pos[2]; c1 <- pos[3]; d1 <- pos[4]; e1 <- pos[5]
  f1 <- pos[min(6, length(pos))]
  0.4 * (X[, a]^2 - 1) +            # quadratic
    0.7 * (X[, b] > 0.5) +           # threshold
    0.4 * X[, c1] * X[, d1] +        # interaction 1
    0.4 * X[, e1] * X[, f1]          # interaction 2
}

# Time-varying part for non-PH: g(x) multiplies log(t), i.e. the log hazard
# ratio of two predictors changes linearly in log time.
nonph_part <- function(X, pos) {
  g <- 0.3 * X[, pos[1]] - 0.3 * X[, pos[2]]
  pmin(pmax(g, -1), 1)  # keep shape + g > 0 (shape = 1.5)
}

# ---- scenario definition ----------------------------------------------------

#' A scenario: list(n, p, event_frac, dgm, rho, censoring, n_clin, ...)
make_scenario <- function(id, n, p, event_frac, dgm = "sparse", rho = 0.3,
                          censoring = "random", n_clin = 0, block_name = "A",
                          shape = 1.5, tau = 10, cens_prop = 0.2,
                          n_test = 10000) {
  sc <- list(id = id, block = block_name, n = n, p = p, event_frac = event_frac,
             dgm = dgm, rho = rho, censoring = censoring, n_clin = n_clin,
             shape = shape, tau = tau, cens_prop = cens_prop, n_test = n_test)
  sc$truth <- make_beta(p, dgm, seed = 1000L + id, n_clin = n_clin)
  sc
}

# Design matrix (clinical + X) for a scenario.
gen_design <- function(sc, n) {
  X <- gen_x(n, sc$p, sc$rho)
  if (sc$n_clin > 0) X <- cbind(gen_clinical(n), X)
  X
}

# Linear predictor and time-varying multiplier for given design.
lin_pred <- function(sc, X) {
  b <- sc$truth$beta
  eta <- drop(X[, names(b), drop = FALSE] %*% b)
  Xp <- X[, grep("^x", colnames(X)), drop = FALSE]
  if (sc$dgm == "nonlinear") eta <- eta + nonlinear_part(Xp, sc$truth$pos)
  g <- if (sc$dgm == "nonph") nonph_part(Xp, sc$truth$pos) else rep(0, nrow(X))
  list(eta = eta, g = g)
}

# Cumulative hazard H(t | x) = lambda * k * exp(eta) * t^(k+g) / (k+g)
cumhaz <- function(t, lp, lambda, k) {
  kk <- k + lp$g
  lambda * k * exp(lp$eta) * t^kk / kk
}

# Draw event times by inversion: H(T) = -log(U)
draw_times <- function(lp, lambda, k) {
  kk <- k + lp$g
  u <- runif(length(lp$eta))
  (-log(u) * kk / (lambda * k * exp(lp$eta)))^(1 / kk)
}

# Censoring: exponential random censoring (rate chosen so that `cens_prop`
# would be censored before tau if no events occurred), plus administrative
# censoring at tau. "covariate" = rate depends on one signal and one noise
# predictor (sensitivity analysis).
draw_censoring <- function(sc, X) {
  rate <- -log(1 - sc$cens_prop) / sc$tau
  if (sc$censoring == "covariate") {
    Xp <- X[, grep("^x", colnames(X)), drop = FALSE]
    noise <- setdiff(seq_len(ncol(Xp)), sc$truth$pos)[1]
    rate <- rate * exp(0.5 * Xp[, sc$truth$pos[1]] + 0.5 * Xp[, noise])
  }
  pmin(rexp(nrow(X), rate), sc$tau)
}

# ---- calibration ------------------------------------------------------------

# Find lambda so that the expected observed event fraction equals the target,
# using a large Monte Carlo sample with common random numbers. Also returns
# the evaluation times t* = 25/50/75th percentiles of observed event times.
calibrate <- function(sc, n_cal = 1e5, seed = 99L) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed + sc$id)
  X <- gen_design(sc, n_cal)
  lp <- lin_pred(sc, X)
  kk <- sc$shape + lp$g
  e <- -log(runif(n_cal))            # common random numbers for T
  C <- draw_censoring(sc, X)
  ev_frac <- function(loglam) {
    T <- (e * kk / (exp(loglam) * sc$shape * exp(lp$eta)))^(1 / kk)
    mean(T <= C)
  }
  f <- function(loglam) ev_frac(loglam) - sc$event_frac
  root <- uniroot(f, c(-30, 10), tol = 1e-6)$root
  lambda <- exp(root)
  T <- (e * kk / (lambda * sc$shape * exp(lp$eta)))^(1 / kk)
  obs <- T[T <= C]
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  sc$lambda <- lambda
  sc$achieved_event_frac <- mean(T <= C)
  sc$tstar <- unname(quantile(obs, c(0.25, 0.5, 0.75)))
  sc
}

# ---- simulate one dataset ---------------------------------------------------

sim_data <- function(sc, n) {
  X <- gen_design(sc, n)
  lp <- lin_pred(sc, X)
  T <- draw_times(lp, sc$lambda, sc$shape)
  C <- draw_censoring(sc, X)
  list(X = X, time = pmin(T, C), status = as.integer(T <= C), lp = lp)
}

# True survival matrix for a simulated dataset at given times (n x length(times))
true_surv <- function(sc, lp, times) {
  sapply(times, function(t) exp(-cumhaz(t, lp, sc$lambda, sc$shape)))
}
