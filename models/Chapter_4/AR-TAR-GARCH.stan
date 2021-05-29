data {
  int<lower=0> N;
  vector[N] y;
  int K;
  int threshold;
  real sigma1;
  real mu1;
}
parameters {
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=-1,upper=1> phi_shock;                 
  real<lower=-1,upper=1> phi_sigma;                  
  real ar0;
  vector<lower=0,upper=1>[K] ar;
  
}
transformed parameters {
  real<lower=0> sigma[N];
  real mu[N];
  for (n in 1:K) {
    sigma[n] = sigma1;
    mu[n] = mu1;
  }
  for (t in (K+1):N) {
    mu[t] = ar0;
    for (k in 1:K) 
      mu[t] += ar[k] * y[t-k];
    sigma[t] = sqrt(alpha0
              + alpha1 * pow(y[t-1] - mu[t-1], 2)
              + beta1 * pow(sigma[t-1], 2)
              + (phi_shock *  pow(y[t-1] - mu[t-1], 2) + phi_sigma * pow(sigma[t-1], 2)) * ((y[t-1] - mu[t-1]) < threshold)
              );
  }
}
model {
  alpha0 ~ normal(0, 1);
  alpha1 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  sigma1 ~ cauchy(0, 1);
  phi_shock ~ normal(0, 2);
  phi_sigma ~ normal(0, 2);
  ar0 ~ normal(0, 1);
  ar ~ normal(0, 1);
  y ~ normal(mu, sigma);    
  
}
generated quantities {
  vector[N] log_lik;
  {
    real log_lik_mu = 0;
    for (t in (K+1):N) {
      log_lik[t] = normal_lpdf(y[t] | mu[t], sigma[t]);  
      log_lik_mu += log_lik[t];
    }       
    log_lik_mu /= N-(K+1);
    for (t in 1:K) 
      log_lik[t] = log_lik_mu;    
  }
}
