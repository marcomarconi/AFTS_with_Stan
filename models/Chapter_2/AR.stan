data {
  int<lower=0> K;
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real intercept;
  real ar[K];
  real<lower=0> sigma;
}
model {
  intercept ~ normal(0, 1);
  ar ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  for (n in (K+1):N) {
    real mu = intercept;
    for (k in 1:K)
      mu += ar[k] * y[n-k];
    y[n] ~ normal(mu, sigma);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_hat;
  {
    for (n in 1:K){
      log_lik[n] = normal_lpdf(y[n] | intercept / (1 - sum(ar)), sigma);
      y_hat[n] = normal_rng(intercept / (1 - sum(ar)), sigma);
    }
    
    for (n in (K+1):N)  {
      real mu = intercept;
      for (k in 1:K)
        mu += ar[k] * y[n-k];
      log_lik[n] = normal_lpdf(y[n] | mu, sigma);
      y_hat[n] = normal_rng(mu, sigma);
    }
  }
}
