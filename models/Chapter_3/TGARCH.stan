data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
parameters {
  real mu;
  real<lower=0> sigma1;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
  real gamma;
}
transformed parameters {
  real<lower=0> sigma[N];
  sigma[1] = sigma1;
  for (t in 2:N)
    sigma[t] = sqrt(alpha0
                    + (alpha1 + gamma*(y[t-1] - mu < 0)) * pow(y[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));
}
model {
  sigma1 ~ cauchy(0, 1);
  gamma ~ normal(0, 1);
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N+M] r_hat;
  vector[N+M] log_lik;
  {
    vector[M+N] z = append_row(y, m);
    vector[M+N] sigma_hat;
    r_hat[1] = normal_rng(mu, sigma1);
    log_lik[1] = normal_lpdf(z[1] | mu, sigma1);
    sigma_hat[1] = sigma[1];
    for (t in 2:(N+M)) {
      sigma_hat[t] = sqrt(alpha0
                    + (alpha1 + gamma*(z[t-1] - mu < 0)) * pow(z[t-1] - mu, 2)
                    + beta1 * pow(sigma_hat[t-1], 2));
      r_hat[t] = normal_rng(mu, sigma_hat[t]);
      log_lik[t] = normal_lpdf(z[t] | mu, sigma_hat[t]);    
    }              
  }
}
