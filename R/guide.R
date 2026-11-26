# Practitioner guide: which model to use for which n, event rate and p
# ---------------------------------------------------------------------------
# Built automatically from the simulation results, so the guide in the paper is
# exactly what the data say. Rules (pre-specified, protocol Section 6):
#   * Prediction: primary measure = MSE of S(t*50 | x) against the truth.
#     In each scenario, the best model is the one with the lowest mean. Any model
#     whose PAIRED difference from the best is within 2 x MCSE is "as good as best".
#   * Inference: a model is "reliable for effect estimation" in a scenario when
#     its 95% interval coverage for true log-HRs is >= 0.925 (500 reps: about
#     2.5 x MCSE below nominal).
#   * Cells of the guide pool scenarios with the same (n, event fraction, p)
#     across the linear-PH data-generating mechanisms; nonlinear and non-PH
#     scenarios form a separate "if you suspect non-linear or time-varying
#     effects" part of the guide.

MODEL_LABELS <- c(
  S1_cox = "Cox", S2_ridge = "Ridge Cox", S3_bayes_normal = "Bayesian Cox (normal)",
  S4_lasso_enet = "Lasso/EN Cox", S5_bayes_horseshoe = "Bayesian Cox (horseshoe)",
  S6_rsf = "RSF", S7_xgb = "XGBoost-Cox", S8_select_fit = "Select-then-fit",
  S8b_select_bayes = "Select-then-fit (Bayes)", S8leak_select_fit = "Select-then-fit (leak)",
  S9_bart = "BART")

MODEL_COLOURS <- c(
  "Cox" = "#7F7F7F", "Ridge Cox" = "#9ECAE1", "Bayesian Cox (normal)" = "#3182BD",
  "Lasso/EN Cox" = "#FDAE6B", "Bayesian Cox (horseshoe)" = "#E6550D", "RSF" = "#74C476",
  "XGBoost-Cox" = "#31A354", "Select-then-fit" = "#BCBDDC", "Select-then-fit (Bayes)" = "#756BB1",
  "Select-then-fit (leak)" = "#D9D9D9", "BART" = "#E7BA52")

# Models eligible for recommendation (the leak variant is illustrative only)
ELIGIBLE <- setdiff(names(MODEL_LABELS), "S8leak_select_fit")

# Per-scenario: best model and the set within 2 x MCSE (paired), for one metric.
best_sets <- function(M, metric = "mse_S", t_index = 2, lower_better = TRUE) {
  m <- M[M$metric == metric & (is.na(M$t_index) | M$t_index == t_index) &
           M$model %in% ELIGIBLE, c("scenario", "rep", "model", "value")]
  out <- lapply(split(m, m$scenario), function(s) {
    w <- reshape(s, idvar = "rep", timevar = "model", direction = "wide")
    vals <- as.matrix(w[, grep("^value\\.", names(w)), drop = FALSE])
    colnames(vals) <- sub("^value\\.", "", colnames(vals))
    sgn <- if (lower_better) 1 else -1
    means <- colMeans(vals, na.rm = TRUE)
    best <- names(which.min(sgn * means))
    diffs <- sgn * (vals - vals[, best])          # >= 0 means worse than best
    md <- colMeans(diffs, na.rm = TRUE)
    se <- apply(diffs, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    se[best] <- 0
    asgood <- names(md)[md <= 2 * pmax(se, 1e-12)]
    data.frame(scenario = s$scenario[1], best = best,
               as_good = paste(sort(asgood), collapse = ";"),
               n_models = ncol(vals))
  })
  do.call(base::rbind, out)
}

# Coverage-based reliability for effect estimation (true predictors only).
# Only models that estimate EVERY coefficient with an interval are eligible:
# select-then-fit reports intervals only for the predictors it selected, so its
# coverage is conditional on selection (reported separately, H6), not a basis
# for recommending it for effect estimation.
INFERENCE_MODELS <- c("S1_cox", "S3_bayes_normal", "S5_bayes_horseshoe")

reliable_inference <- function(C, threshold = 0.925) {
  if (is.null(C) || !nrow(C)) return(NULL)
  C <- C[C$true != 0 & !is.na(C$lo) & C$model %in% INFERENCE_MODELS, ]
  C$cover <- C$lo <= C$true & C$true <= C$hi
  a <- aggregate(cover ~ scenario + model, C, mean)
  a$reliable <- a$cover >= threshold
  a
}

#' Build the guide.
#' @param M metrics data.frame (all replications), @param C coefficients,
#' @param design design table (id, block, n, p, event_frac, dgm, ...)
build_guide <- function(M, C, design) {
  bs <- best_sets(M)
  bs <- merge(bs, design[, c("id", "block", "n", "p", "event_frac", "dgm")],
              by.x = "scenario", by.y = "id")
  inf <- reliable_inference(C)
  bs$linear <- bs$dgm %in% c("sparse", "dense")
  # pool scenarios into guide cells
  cell_key <- function(d) paste(d$n, d$event_frac, d$p, d$linear)
  cells <- lapply(split(bs, cell_key(bs)), function(d) {
    tab <- table(unlist(strsplit(d$as_good, ";")))
    share <- tab / nrow(d)                         # share of scenarios where "as good as best"
    rec <- names(share)[share == max(share)]
    also <- names(share)[share >= 0.5 & !(names(share) %in% rec)]
    present <- intersect(ELIGIBLE, unique(M$model[M$scenario %in% d$scenario & M$metric == "mse_S"]))
    avoid <- setdiff(present, names(share))       # never as good as the best in this cell
    inf_ok <- if (!is.null(inf)) {
      ii <- inf[inf$scenario %in% d$scenario, ]
      ok <- tapply(ii$reliable, ii$model, all)
      names(ok)[ok %in% TRUE]
    } else character(0)
    n_events <- round(d$n[1] * d$event_frac[1])
    best1 <- names(sort(table(d$best), decreasing = TRUE))[1]   # lowest mean most often
    n_tied <- length(rec) - (best1 %in% rec)
    data.frame(
      n = d$n[1], event_frac = d$event_frac[1], p = d$p[1], linear = d$linear[1],
      events = n_events, epv = round(n_events / d$p[1], 1), scenarios = nrow(d),
      best_mean = MODEL_LABELS[best1], n_tied = n_tied,
      recommended = paste(MODEL_LABELS[rec], collapse = ", "),
      also_fine = paste(MODEL_LABELS[also], collapse = ", "),
      avoid = paste(MODEL_LABELS[avoid], collapse = ", "),
      reliable_intervals = paste(MODEL_LABELS[inf_ok], collapse = ", "),
      stringsAsFactors = FALSE)
  })
  g <- do.call(base::rbind, cells)
  g[order(!g$linear, g$p, g$n, g$event_frac), ]
}

# Markdown version of the guide for the paper / supplement.
guide_markdown <- function(g, title = "Practitioner guide") {
  big <- function(x) formatC(x, format = "d", big.mark = ",")
  fmt_row <- function(r) {
    eq <- setdiff(trimws(strsplit(r$recommended, ",")[[1]]), r$best_mean)
    sprintf("| %s | %s | %s | %s | %s | %s | %s |",
      big(r$n), paste0(round(100 * r$event_frac, 1), "% (", r$events, ")"), big(r$p), r$epv,
      r$best_mean, if (length(eq)) paste(eq, collapse = ", ") else "-",
      ifelse(r$reliable_intervals == "", "none", r$reliable_intervals))
  }
  hdr <- c("| n | Events | p | EPV | Lowest prediction error | Statistically equivalent (within 2 MCSE) | Reliable 95% intervals for effects |",
           "| --- | --- | --- | --- | --- | --- | --- |")
  lin <- g[g$linear, ]; nl <- g[!g$linear, ]
  out <- c(paste("#", title), "",
           "Built automatically from the simulation results with pre-specified rules (R/guide.R).",
           "EPV = events per candidate predictor. Prediction is judged on the MSE of S(t*50 | x);",
           "\"reliable intervals\" means 95% interval coverage of true log hazard ratios >= 92.5% in every scenario of the cell,",
           "among methods that estimate all effects (Cox, Bayesian Cox normal / horseshoe).",
           "", "## Linear proportional-hazards settings", "", hdr,
           vapply(seq_len(nrow(lin)), function(i) fmt_row(lin[i, ]), ""))
  if (nrow(nl)) out <- c(out, "", "## If you suspect non-linear or time-varying effects", "", hdr,
                         vapply(seq_len(nrow(nl)), function(i) fmt_row(nl[i, ]), ""))
  out
}

# Heatmap: which model family wins in each cell.
guide_heatmap <- function(g, file) {
  suppressPackageStartupMessages(library(ggplot2))
  big <- function(x) formatC(x, format = "d", big.mark = ",")
  g$row <- sprintf("p = %s, n = %s", big(g$p), big(g$n))
  g$row <- factor(g$row, levels = unique(g$row[order(g$p, g$n)]))
  g$col <- factor(paste0(100 * g$event_frac, "% events"),
                  levels = paste0(sort(unique(100 * g$event_frac)), "% events"))
  g$panel <- factor(ifelse(g$linear, "Linear PH truth", "Non-linear / non-PH truth"),
                    c("Linear PH truth", "Non-linear / non-PH truth"))
  g$txt <- sprintf("%s\nEPV %s%s", sub(" Cox$", "", g$best_mean), g$epv,
                   ifelse(g$n_tied > 0, sprintf(" | %d tied", g$n_tied), ""))
  p <- ggplot(g, aes(col, row, fill = best_mean)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = txt), size = 3, colour = "grey10", lineheight = 0.9) +
    facet_grid(. ~ panel, scales = "free", space = "free") +
    scale_fill_manual(values = MODEL_COLOURS) +
    labs(x = NULL, y = NULL, fill = "Lowest prediction error",
         title = "Which model predicts best, by sample size, event rate and dimension",
         subtitle = "Tile = model with the lowest MSE of S(t*50); 'tied' = other models within 2 Monte Carlo SE") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))
  ggsave(file, p, width = 10, height = 2.2 + 0.55 * length(unique(g$row)), dpi = 150)
  invisible(p)
}
