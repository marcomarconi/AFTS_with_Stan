// this is a frankenstein between STAN manual AR and GARCH

data {
  int<lower=0> T;
  real r[T];
  int u[T];
  int ar_order;
}
parameters {
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
  real<lower=0> sigma1;
  real gamma;
  real ar0;
  vector[ar_order] ar;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(
                    alpha1 * pow(r[t-1] - ar0, 2)
                    + beta1 * pow(sigma[t-1], 2)
                    + gamma * (1 - u[t])

                    );
}
model {
  alpha1 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  sigma1 ~ cauchy(0, 1);
  ar0 ~ normal(0, 1);
  ar ~ normal(0, 1);
  for (n in (ar_order+1):T) {
    real mu = ar0;
    for (k in 1:ar_order)
      mu += ar[k] * r[n-k];
    r[n] ~ normal(mu, sigma[n]);    
  }
  
}
generated quantities {
  vector[T-ar_order] r_hat;
  vector[T-ar_order] log_lik;
  // for (n in 1:ar_order) {
  //   r_hat[n] = normal_rng(ar0, alpha0);
  //   log_lik[n] = normal_lpdf(r[n] | ar0, alpha0);
  // }
  for (n in (ar_order+1):T) {
    real mu = ar0;
    for (k in 1:ar_order)
      mu += ar[k] * r[n-k];
    r_hat[n-ar_order] = normal_rng(mu, sigma[n]);
    log_lik[n-ar_order] = normal_lpdf(r[n] | mu, sigma[n]);
  }
}
