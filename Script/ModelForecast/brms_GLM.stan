// generated with brms 2.14.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  vector[N] offsets;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
}
transformed parameters {
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N) + offsets;
    target += poisson_log_glm_lpmf(Y | X, mu, b);
  }
  // priors including all constants
}
generated quantities {
}