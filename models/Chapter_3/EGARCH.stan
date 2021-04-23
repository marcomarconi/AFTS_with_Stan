data {
  int<lower=0> T;
  real r[T];
  real<lower=0> sigma1;
}
parameters {
  real mu;
  //real<lower=0> sigma1;
  real alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
  real gamma;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T) {
    real shock = r[t-1] - mu;
    sigma[t] = sqrt(exp(alpha0
                    + alpha1 * ((fabs(shock) + gamma*shock) / sigma[t-1] )
                    + beta1 * log(pow(sigma[t-1], 2))));
  }  
}
model {
  mu ~ normal(0, 10);
  alpha0 ~ normal(0, 2);
  gamma ~ normal(0, 2);
  r ~ normal(mu, sigma);
}
generated quantities {
  vector[T-1] r_hat;
  vector[T-1] log_lik;
  for (t in 1:(T-1)) {
      r_hat[t] = normal_rng(mu, sigma[t]);
      log_lik[t] = normal_lpdf(r[t] | mu, sigma[t]);
  }
}
