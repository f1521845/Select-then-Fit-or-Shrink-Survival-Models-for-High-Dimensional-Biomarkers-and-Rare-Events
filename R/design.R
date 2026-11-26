# Full simulation design (protocol Section 4), after the pilot decisions:
# 83 scenarios = A 54, B 16, C 4, D 9
# ---------------------------------------------------------------------------
# Returns one row per scenario, with the replication counts:
#   reps_main  = replications for the non-Bayesian strategies
#   reps_bayes = replications for the Bayesian arms (S3, S5, S8b, S9 and the
#                prior variants): 200 when n = 5,000 in Block A and in Blocks
#                B and C (coverage MCSE about 1.5 points); otherwise 500.

source("R/dgm.R")

design_table <- function() {
  A <- expand.grid(n = c(200, 1000, 5000), p = c(10, 50), event_frac = c(0.1, 0.3, 0.7),
                   dgm = c("sparse", "nonlinear", "nonph"), stringsAsFactors = FALSE)
  A$block <- "A"; A$rho <- 0.3; A$censoring <- "random"; A$n_clin <- 0; A$n_test <- 10000
  A$variants <- FALSE

  B <- expand.grid(n = c(200, 1000), p = c(500, 2000), event_frac = c(0.1, 0.4),
                   dgm = c("sparse", "dense"), stringsAsFactors = FALSE)
  B$block <- "B"; B$rho <- 0.3; B$censoring <- "random"; B$n_clin <- 0; B$n_test <- 10000
  B$variants <- FALSE

  # Block C, redesigned after the pilot: n = 10,000 (was 20,000), event fractions
  # 1% and 4% so the event counts (~100 and ~400) match the original design.
  C <- expand.grid(n = 10000, p = 2000, event_frac = c(0.01, 0.04), rho = c(0.3, 0.7),
                   stringsAsFactors = FALSE)
  C$dgm <- "sparse"; C$block <- "C"; C$censoring <- "random"; C$n_clin <- 3
  C$n_test <- 50000; C$variants <- FALSE

  # Block D: one factor at a time from a reference scenario in A and in B.
  refA <- data.frame(n = 1000, p = 50, event_frac = 0.3, dgm = "sparse")
  refB <- data.frame(n = 200, p = 500, event_frac = 0.4, dgm = "sparse")
  refH <- data.frame(n = 5000, p = 10, event_frac = 0.3, dgm = "sparse")   # high EPV (150), for H5
  mk <- function(ref, rho, cens, variants)
    cbind(ref, block = "D", rho = rho, censoring = cens, n_clin = 0, n_test = 10000,
          variants = variants)
  D <- rbind(
    mk(refA, 0.3, "random", TRUE),       # prior grid on reference A
    mk(refB, 0.3, "random", TRUE),       # prior grid on reference B
    mk(refH, 0.3, "random", TRUE),       # prior grid at high EPV (H5: robust when EPV > 20)
    mk(refA, 0.0, "random", FALSE), mk(refA, 0.7, "random", FALSE),
    mk(refB, 0.0, "random", FALSE), mk(refB, 0.7, "random", FALSE),
    mk(refA, 0.3, "covariate", FALSE), mk(refB, 0.3, "covariate", FALSE))

  cols <- c("block", "n", "p", "event_frac", "dgm", "rho", "censoring", "n_clin", "n_test", "variants")
  d <- rbind(A[cols], B[cols], C[cols], D[cols])
  d$id <- seq_len(nrow(d))
  d$reps_main <- 500
  d$reps_bayes <- ifelse(d$block %in% c("B", "C") | d$n == 5000 |
                           (d$block == "D" & d$p == 500), 200, 500)
  d[c("id", cols, "reps_main", "reps_bayes")]
}

scenario_from_row <- function(r) {
  make_scenario(r$id, n = r$n, p = r$p, event_frac = r$event_frac, dgm = r$dgm,
                rho = r$rho, censoring = r$censoring, n_clin = r$n_clin,
                block_name = r$block, n_test = r$n_test)
}
