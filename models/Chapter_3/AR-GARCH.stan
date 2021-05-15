data {
  int<lower=0> N;
  vector[N] y;
  int K;

}
parameters {
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  real ar0;
  vector[K] ar;
  
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
  for (n in (K+1):N) {
    real mu = ar0;
    for (k in 1:K)
      mu += ar[k] * y[n-k];
    y[n] ~ normal(mu, sigma[n]);    
  }
  
}
generated quantities {
  vector[N] log_lik;  
  {
    real log_lik_mu = 0;
    for (t in 1:K)
      log_lik[t] = normal_lpdf(y[t] | ar0 / (1 - sum(ar)), sigma[t]);  
    for (t in (K+1):N) {
      real mu = ar0;
      for (k in 1:K)
        mu += ar[k] * y[t-k];
      log_lik[t] = normal_lpdf(y[t] | mu, sigma[t]);  
      log_lik_mu += log_lik[t];
    }       
    log_lik_mu /= N-(K+1);
    for (t in 1:K) 
      log_lik[t] = log_lik_mu;
  }
}
