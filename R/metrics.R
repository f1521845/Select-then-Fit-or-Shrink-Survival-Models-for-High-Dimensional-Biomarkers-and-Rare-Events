# Performance measures (protocol Section 6)
# ---------------------------------------------------------------------------

harrell_c <- function(time, status, risk) {
  survival::concordance(survival::Surv(time, status) ~ risk, reverse = TRUE)$concordance
}

uno_c <- function(time, status, risk, tau) {
  survival::concordance(survival::Surv(time, status) ~ risk, reverse = TRUE,
                        timewt = "n/G2", ymax = tau)$concordance
}

# IPCW Brier score on a grid of times and its integral (Graf et al. 1999).
ibs <- function(time, status, S, grid) {
  km_c <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  G <- stepfun(km_c$time, c(1, km_c$surv))
  G_tminus <- G(time - 1e-10)
  bs <- vapply(seq_along(grid), function(j) {
    t <- grid[j]
    w_event <- (time <= t & status == 1) / pmax(G_tminus, 1e-8)
    w_alive <- (time > t) / max(G(t), 1e-8)
    mean(w_event * S[, j]^2 + w_alive * (1 - S[, j])^2)
  }, 0)
  # trapezoidal integral from grid[1] to grid[end], scaled by the range
  sum(diff(grid) * (head(bs, -1) + tail(bs, -1)) / 2) / (max(grid) - min(grid))
}

# Calibration against the known truth at one time point:
#  - slope: lm(cloglog true risk ~ cloglog predicted risk); 1 = ideal,
#    < 1 = predictions too extreme (overfitting)
#  - ICI: mean |predicted - smoothed true| (Austin et al. 2020 idea, using truth)
calib_truth <- function(S_pred, S_true, n_sub = 3000) {
  cl <- function(s) log(-log(pmin(pmax(s, 1e-8), 1 - 1e-8)))
  i <- if (length(S_pred) > n_sub) sample(length(S_pred), n_sub) else seq_along(S_pred)
  slope <- unname(coef(lm(cl(S_true[i]) ~ cl(S_pred[i])))[2])
  sm <- tryCatch(loess(S_true[i] ~ S_pred[i], span = 0.75, degree = 2),
                 error = function(e) NULL)
  ici <- if (is.null(sm)) NA_real_ else mean(abs(S_pred[i] - fitted(sm)))
  c(calib_slope = slope, ici = ici)
}

selection_metrics <- function(selected, active, all_vars) {
  sel <- all_vars %in% selected; act <- all_vars %in% active
  tp <- sum(sel & act); fp <- sum(sel & !act); fn <- sum(!sel & act); tn <- sum(!sel & !act)
  mcc_den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  c(tpr = tp / max(1, tp + fn),
    fdr = if (tp + fp == 0) 0 else fp / (tp + fp),
    mcc = if (mcc_den == 0) 0 else (tp * tn - fp * fn) / mcc_den,
    n_selected = tp + fp)
}
