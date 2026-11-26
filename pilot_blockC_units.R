# Block C unit timings (one fit of each component) for extrapolation.
source("R/run.R")
sc <- calibrate(make_scenario(203, n = 20000, p = 2000, event_frac = 0.005, dgm = "sparse",
                              n_clin = 3, block_name = "C", n_test = 20000), n_cal = 5e4)
RNGkind("L'Ecuyer-CMRG"); set.seed(1e6 * sc$id + 1); tr <- sim_data(sc, sc$n)
cat("events:", sum(tr$status), "\n")
st <- standardiser(tr$X); Z <- st$apply(tr$X); y <- Surv(tr$time, tr$status)
tm <- function(label, expr) { t <- system.time(force(expr))[["elapsed"]]; cat(sprintf("%-45s %8.1fs\n", label, t)); t }
out <- list()
out$lasso_path <- tm("lasso path, 1 fit (100 lambdas)", glmnet(Z, y, family = "cox", alpha = 1, standardize = FALSE))
out$ridge_path20 <- tm("ridge path, 1 fit (20 lambdas)", glmnet(Z, y, family = "cox", alpha = 0, nlambda = 20, standardize = FALSE))
suppressPackageStartupMessages(library(xgboost))
lab <- ifelse(tr$status == 1, tr$time, -tr$time); d <- xgb.DMatrix(tr$X, label = lab, nthread = 2)
out$xgb_200 <- tm("xgboost 200 rounds, depth 3, eta 0.05", xgb.train(list(objective = "survival:cox", max_depth = 3, eta = 0.05, nthread = 2), d, nrounds = 200))
suppressPackageStartupMessages(library(randomForestSRC)); options(rf.cores = 2)
df <- data.frame(time = tr$time, status = tr$status, tr$X, check.names = FALSE)
out$rsf_50 <- tm("RSF 50 trees, mtry sqrt(p), nodesize 15", rfsrc(Surv(time, status) ~ ., df, ntree = 50, mtry = 45, nodesize = 15))
suppressPackageStartupMessages(library(rstan))
o <- order(tr$time, decreasing = TRUE); ts <- tr$time[o]
dat <- list(N = nrow(Z), P = ncol(Z), X = Z[o, ], status = tr$status[o], last = seq_along(ts),
            prior = 1L, sigma = 1, tau0 = 10 / (ncol(Z) - 10) / sqrt(sum(tr$status)), slab_scale = 2, slab_df = 4)
out$hs_40iter <- tm("horseshoe Stan, 1 chain, 20 warmup + 20 iter", sampling(get_stan_model(), data = dat, chains = 1, iter = 40, warmup = 20, refresh = 0, pars = "beta"))
saveRDS(out, "results/pilot/timing/blockC_units.rds")
