// Bayesian Cox model on the partial likelihood (Breslow handling of ties)
// with either a normal prior or a regularised horseshoe prior
// (Piironen & Vehtari 2017), non-centred parameterisation.
//
// Data must be sorted by DEcreasing time, so that the risk set of subject i
// is {1, ..., last[i]}, where last[i] is the last index of i's tie group.
data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;               // standardised predictors
  array[N] int<lower=0, upper=1> status;
  array[N] int<lower=1, upper=N> last;
  int<lower=0, upper=1> prior;  // 0 = normal, 1 = regularised horseshoe
  int<lower=0, upper=1> hier;   // normal prior: 1 = SD learned (half-normal hyperprior), 0 = fixed
  real<lower=0> sigma;          // normal prior: fixed SD (hier = 0) or hyperprior scale (hier = 1)
  real<lower=0> tau0;           // horseshoe global scale
  real<lower=0> slab_scale;     // horseshoe slab scale
  real<lower=0> slab_df;        // horseshoe slab df
}
transformed data {
  int H = prior == 1 ? 1 : 0;   // horseshoe parameters only when needed
  int S = (prior == 0 && hier == 1) ? 1 : 0;  // learned SD only when needed
}
parameters {
  vector[P] z;
  array[H] real<lower=0> tau_raw;
  vector<lower=0>[H * P] lambda;
  array[H] real<lower=0> caux;
  array[S] real<lower=0> sig;
}
transformed parameters {
  vector[P] beta;
  if (prior == 0) {
    if (S == 1) beta = sig[1] * z;
    else beta = sigma * z;
  } else {
    real tau = tau0 * tau_raw[1];
    real c = slab_scale * sqrt(caux[1]);
    vector[P] lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda)));
    beta = z .* lambda_tilde * tau;
  }
}
model {
  vector[N] eta = X * beta;
  real m = max(eta);
  vector[N] log_risk = m + log(cumulative_sum(exp(eta - m)));
  z ~ std_normal();
  if (S == 1) sig[1] ~ normal(0, sigma);      // half-normal hyperprior on the SD
  if (prior == 1) {
    tau_raw[1] ~ cauchy(0, 1);
    lambda ~ cauchy(0, 1);
    caux[1] ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  }
  for (i in 1:N) {
    if (status[i] == 1) target += eta[i] - log_risk[last[i]];
  }
}
