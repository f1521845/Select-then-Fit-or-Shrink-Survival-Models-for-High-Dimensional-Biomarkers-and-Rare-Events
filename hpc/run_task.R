# HPC runner for the full simulation.
#
#   Rscript hpc/run_task.R calibrate        # once: calibrate every scenario, write the task list
#   Rscript hpc/run_task.R task <k>         # array job k: one scenario x one batch of replications
#   Rscript hpc/run_task.R status           # how many replications are done per scenario
#
# Run from the project root. Each replication is saved as its own file and is
# skipped if it already exists, so failed or timed-out tasks can be resubmitted.

source("R/design.R"); source("R/run.R")
args <- commandArgs(trailingOnly = TRUE)
cmd <- if (length(args)) args[1] else "status"
BATCH <- 10                                    # replications per array task
sc_dir <- "results/scenarios"; out_dir <- "results/main"
dir.create(sc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stan_full <- list(chains = 4, iter = 2000, warmup = 1000,
                  cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

if (cmd == "calibrate") {
  d <- design_table()
  write.csv(d, "results/design.csv", row.names = FALSE)
  for (i in seq_len(nrow(d))) {
    f <- file.path(sc_dir, sprintf("sc_%03d.rds", d$id[i]))
    if (!file.exists(f)) saveRDS(calibrate(scenario_from_row(d[i, ])), f)
  }
  tasks <- do.call(base::rbind, lapply(seq_len(nrow(d)), function(i) {
    starts <- seq(1, d$reps_main[i], by = BATCH)
    data.frame(scenario = d$id[i], rep_from = starts, rep_to = pmin(starts + BATCH - 1, d$reps_main[i]))
  }))
  tasks$task <- seq_len(nrow(tasks))
  write.csv(tasks, "results/tasks.csv", row.names = FALSE)
  message(nrow(d), " scenarios, ", nrow(tasks), " array tasks (use --array=1-", nrow(tasks), ")")
}

if (cmd == "task") {
  k <- as.integer(if (length(args) > 1) args[2] else Sys.getenv("SLURM_ARRAY_TASK_ID"))
  tasks <- read.csv("results/tasks.csv"); d <- read.csv("results/design.csv")
  tk <- tasks[tasks$task == k, ]; row <- d[d$id == tk$scenario, ]
  sc <- readRDS(file.path(sc_dir, sprintf("sc_%03d.rds", tk$scenario)))
  nthreads <- stan_full$cores
  for (r in tk$rep_from:tk$rep_to) {
    f <- file.path(out_dir, sprintf("s%03d_r%03d.rds", sc$id, r))
    if (file.exists(f)) next
    models <- ALL_MODELS
    if (isTRUE(row$variants)) models <- c(models, PRIOR_VARIANTS)
    if (r > row$reps_bayes) models <- setdiff(models, BAYES_MODELS)
    message(sprintf("[%s] scenario %d rep %d (%d models)", format(Sys.time(), "%H:%M:%S"), sc$id, r, length(models)))
    res <- run_rep(sc, r, models = models, stan_args = stan_full, xgb_configs = 50)
    saveRDS(res, f)
  }
}

if (cmd == "status") {
  d <- read.csv("results/design.csv")
  done <- table(as.integer(sub("^s(\\d+)_r.*", "\\1", list.files(out_dir, "^s\\d+_r\\d+\\.rds$"))))
  d$done <- as.integer(done[as.character(d$id)]); d$done[is.na(d$done)] <- 0
  print(d[, c("id", "block", "n", "p", "event_frac", "dgm", "reps_main", "done")], row.names = FALSE)
  message(sum(d$done), " of ", sum(d$reps_main), " replications done")
}
