# Select-then-Fit or Shrink? Survival Models for High-Dimensional Biomarkers and Rare Events

Analysis code for a preregistered neutral comparison study: a simulation study (83 scenarios), five real-data applications and a practitioner guide.

- **Author:** Raghad Alhuthil
- **Preregistration:** Open Science Framework (link to be added)
- **Status:** code frozen at registration, before the main simulation. The release tagged `v1.0-registration` is the version registered on OSF. Any later change will be listed in the deviations table of the paper.

## Paper workflow

| Step | File | What it produces |
| --- | --- | --- |
| 1. Simulation + guide | `01_simulation.R` | Summary tables with Monte Carlo SEs, the hypothesis table (registered decision rules for H1, H2, H3a, H3b, H4, H5 and H6), figures, and the **practitioner guide** (CSV, Markdown, heatmap) in `results/<mode>/paper/` |
| 2. SUPPORT2 | `02_SUPPORT2.R` | Large clinical cohort with EHR missingness (H1, H4) |
| 3. PBC | `03_PBC.R` | Small trial, about 7 events per predictor (H2, H5) |
| 4. DLBCL | `04_DLBCL.R` | p = 7,399 genes, n = 240 (H3, H6) |
| 5. CHOP → R-CHOP | `05_CHOP_RCHOP.R` | External validation across treatment eras (H6 optimism gap, H3) |
| 6. GSE1992 | `06_GSE1992.R` | Extreme low-event stress test (35 events, descriptive) |

**How to run (RStudio):** set the working directory to the `survsim` folder (Session > Set Working Directory > To Source File Location), open a script and click **Source**. Run `01_simulation.R` first, because each real-data script ends by checking its result against the practitioner guide. Tables print to the Console, figures appear in the Plots pane, and everything is saved in `results/realdata/<dataset>/`.

### Modes

- `01_simulation.R`: `SURVSIM_MODE=demo` (small local run that exercises every step; not for the paper), `local` (a chosen subset of scenarios), or `collect` (summarise the HPC run).
- Real-data scripts: `QUICK <- TRUE` (default in the R scripts) for a fast check of about 5–10 minutes; set `QUICK <- FALSE` for the full protocol analysis (5-fold CV × 10 repeats, m = 5 imputations, full MCMC; hours to days). Every outer fold is cached, so an interrupted run resumes. Quick-mode results are saved separately (`results/realdata/<dataset>/quick/`) and are not for interpretation.

## Layout

| File | What it does |
| --- | --- |
| `R/dgm.R` | Data-generating mechanisms: block-correlated predictors, Weibull PH times (sparse, dense, nonlinear, non-PH), censoring calibrated to a target event fraction, true S(t \| x) |
| `R/models.R` | Strategies S1–S9: Cox, ridge, lasso/elastic net, Bayesian Cox (hierarchical normal and regularised horseshoe, Stan), RSF, XGBoost-Cox, select-then-fit (+ Bayesian and leak variants), BART survival |
| `R/metrics.R` | Harrell's and Uno's C, IPCW integrated Brier score, truth-based calibration slope and ICI, selection metrics |
| `R/run.R` | `run_rep()`: one simulation replication, with per-replication random-number streams |
| `R/design.R` | Full design: 83 scenarios (A 54, B 16, C 4, D 9) with replication counts |
| `R/guide.R` | Practitioner guide built from the simulation results (pre-specified rules) |
| `R/realdata.R` | Real-data engine: repeated CV with imputation, screening and tuning inside folds; external validation; observed-outcome metrics (C, Uno's C, time-dependent AUC, IBS, calibration slope, ICI); selection stability; plots |
| `stan/bayes_cox.stan` | Bayesian Cox partial likelihood with a normal (fixed or hierarchical) or regularised horseshoe prior |
| `data/prepare_data.py` | Exports SUPPORT2, DLBCL and GSE1992 from SurvSet to CSV (already done: the CSV files are included in `data/`) |
| `hpc/run_task.R`, `hpc/submit.slurm` | HPC runner: calibrate once, then one SLURM array task per scenario × 10 replications; resumable |
| `pilot.R`, `summarise.R` | Pilot study (code and runtime checks) |

## Requirements

- R ≥ 4.3 with `survival`, `glmnet`, `rstan` (≥ 2.32), `randomForestSRC`, `xgboost` (≥ 3.0), `BART`, `mice`, `ggplot2`.
- Python ≥ 3.9 with `SurvSet`, only if you want to re-export the data (not needed: the CSV files are included).

The Stan model compiles automatically the first time a Bayesian model is fitted (1–2 minutes), and is then reused.

Check the packages (the code stops with an install command if any is missing):

```r
install.packages(c("survival", "glmnet", "rstan", "randomForestSRC", "xgboost", "BART", "mice", "ggplot2"))
```

Results are cached (simulation replications in `results/<mode>/reps/`, real-data folds in `results/realdata/<dataset>/<quick|full>/`). If a run happened while a package was missing, delete that folder before re-running, or the cached results without that model are reused.

## Data

All data are public and de-identified, and are included in `data/`:

| File | Source |
| --- | --- |
| `support2.csv`, `DLBCL.csv`, `gse1992.csv` | SurvSet repository, exported with `data/prepare_data.py` |
| `chop.rda`, `rchop.rda` | R package `bujar` |
| PBC | loaded directly from the R package `survival` |

## Running the full simulation on a cluster

```sh
Rscript hpc/run_task.R calibrate      # writes results/design.csv and results/tasks.csv
mkdir -p logs && sbatch hpc/submit.slurm   # edit --array to the task count first
Rscript hpc/run_task.R status         # progress
SURVSIM_MODE=collect Rscript 01_simulation.R   # tables, figures, practitioner guide
```

## Running the pilot

```sh
Rscript pilot.R correctness
Rscript summarise.R results/pilot/correctness
Rscript pilot.R timing
```

Every replication is seeded from `1e6 * scenario_id + rep` (L'Ecuyer-CMRG), so any single replication can be re-run in isolation.

## Implementation decisions taken during the pilot (to report as protocol clarifications)

1. **Bayesian Cox likelihood:** partial likelihood (Breslow ties) in Stan, with a Breslow baseline computed per posterior draw for predictions. This replaces the M-spline or piecewise-exponential options listed in the protocol: it needs no baseline hazard parameters and was fast enough in the pilot. Consequence: intervals for S(t | x) reflect uncertainty in β but not in the baseline hazard.
2. **Horseshoe global scale:** τ₀ = p₀ / (p − p₀) / √(number of events). Piironen & Vehtari use σ/√n; for Cox models the information scales with events, not n.
3. **p = 10 scenarios:** the number of true predictors is min(10, p/2), so 5 when p = 10 (otherwise every predictor would be true and null coefficients could not be studied). The horseshoe's p₀ follows the same rule.
4. **Evaluation times t\*:** the 25th, 50th and 75th percentiles of *observed event times* (from a 100,000-subject calibration sample). This replaces "percentiles of the true event-time distribution", which lie beyond follow-up when events are rare.
5. **Calibration with known truth:** the slope is from regressing cloglog(true risk) on cloglog(predicted risk) at t\*50, and the ICI is the mean |predicted − loess-smoothed true| survival. Both use the known truth instead of censored outcomes.
6. **H6 "local split":** each replication also draws a held-out set of size n/4. This mimics an 80/20 split. S8 and S8-leak report their C-index on it, alongside the true C-index from the large test set.
7. **MCMC:** 4 chains, 1,000 + 1,000 iterations, adapt_delta = 0.99. If R-hat ≥ 1.01 or ESS ≤ 400, the fit is re-run once with twice the iterations. `rstan`'s `n_eff` is used for ESS.
8. **XGBoost SHAP:** the bias column is dropped before ranking. xgboost 3.x names it "(Intercept)", not "BIAS".
9. **Normal prior (S3, S8b): hierarchical.** β ~ N(0, s²) with s ~ half-normal(0, 1), so the amount of shrinkage is learned from the data. This is the Bayesian counterpart of ridge tuned by CV. *Informed by the pilot:* a fixed σ = 1 overfitted when p ≥ events (p = 500: apparent C = 1.000, calibration slope 0.36–0.44). A fixed σ = 1/√p over-shrank (slope 1.7). The hierarchical prior gave slope 1.03 in a one-replication check. The two fixed priors are kept as H5 variants (`S3v_fixed1`, `S3v_scaled`), along with a wider hyperprior (`S3v_hyper2.5`).
10. **Replications:** 500 for all strategies, except that the Bayesian arms (S3, S5, S8b, S9 and the prior variants) run only the first 200 replications when n = 5,000 (Block A) and in Blocks B, C and the p = 500 Block D scenarios.
11. **RSF tuning:** the mtry grid is {√p, p/3} for p < 500, and {√p, 2√p} for p ≥ 500.
12. **Block C:** n = 10,000 (not 20,000), with event fractions of 1% and 4%, keeping roughly the original 100 and 400 events.
13. **Lasso/elastic-net path cap:** for p ≥ 500, glmnet paths stop at 200 non-zero coefficients (`pmax = 200`). Uncapped paths at n = 20,000 did not finish in 20 minutes. With the cap, one path at n = 10,000 takes under 1 second. Ridge is never capped. A 20-λ ridge path at n = 10,000 takes about 80 seconds, so ridge is the costliest frequentist model in Block C.
14. **Prior-sensitivity p₀:** capped at p/2 (so p₀ = 50 becomes 25 when p = 50).
15. **Block D high-EPV prior-grid reference** (n = 5,000, p = 10, 30% events; EPV 150), added so that the "robust when EPV > 20" half of H5 can be tested.
16. **Real data, imputation:** `mice` within training folds; test rows are imputed from the training-row models (`ignore`), and the outcome is not used, so test outcomes cannot leak. When there are more than 200 predictors, the few incomplete columns get training-fold median/mode imputation.
17. **Real data, screening:** the top K genes by univariate Cox score test are selected inside each training fold and shared by all models. Clinical covariates are never screened out (GSE1992).
18. **Zero survival times** (CHOP) are moved to half the smallest positive time for XGBoost only, whose Cox objective encodes censoring as a negative label.

