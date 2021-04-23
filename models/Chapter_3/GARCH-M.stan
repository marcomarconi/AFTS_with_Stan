data {
  int<lower=0> T;
  real r[T];
}
parameters {
  real mu;
  real premium;
  real<lower=0> sigma1;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(r[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));
}
model {
  premium ~ normal(0, 1);
  sigma1 ~ cauchy(0, 1);
  for(t in 2:T)
    r[t] ~ normal(mu + premium*pow(sigma[t], 2) , sigma[t]);
}
generated quantities {
  vector[T-1] r_hat;
  vector[T-1] log_lik;
  for (t in 1:(T-1)) {
      r_hat[t] = normal_rng(mu + premium*pow(sigma[t], 2), sigma[t]);
      log_lik[t] = normal_lpdf(r[t] | mu + premium*pow(sigma[t], 2), sigma[t]);
  }
}
