// this is a frankenstein between SNAN manual AR and GARCH

data {
  int<lower=0> N;
  vector[N] y;
  int Kar;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
parameters {
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  real ar0;
  vector[Kar] ar;
  
}
transformed parameters {
  real<lower=0> sigma[N];
  sigma[1] = sigma1;
  for (t in 2:N)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-1] - ar0, 2)
                    + beta1 * pow(sigma[t-1], 2));
}
model {
  alpha0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  sigma1 ~ cauchy(0, 1);
  ar0 ~ normal(0, 1);
  ar ~ normal(0, 1);
  for (n in (Kar+1):N) {
    real mu = ar0;
    for (k in 1:Kar)
      mu += ar[k] * y[n-k];
    y[n] ~ normal(mu, sigma[n]);    
  }
  
}
generated quantities {
  vector[N+M] log_lik;
  {
    vector[M+N] z = append_row(y, m);
    vector[M+N] sigma_hat;
    real log_lik_mu = 0;

    for (t in 1:Kar) 
      sigma_hat[t] = sigma[t];
    for (t in (Kar+1):(N+M)) {
      real mu_hat = ar0;
      sigma_hat[t] = sqrt(alpha0
                    + alpha1 * pow(z[t-1] - ar0, 2)
                    + beta1 * pow(sigma_hat[t-1], 2));
      ;
      for (k in 1:Kar)
        mu_hat += ar[k] * z[t-k];
      log_lik[t] = normal_lpdf(z[t] | mu_hat, sigma_hat[t]);  
      log_lik_mu += log_lik[t];
    }       
    log_lik_mu /= M+N-(Kar+1);
    for (t in 1:Kar) 
      log_lik[t] = log_lik_mu;    
    
  }
}
