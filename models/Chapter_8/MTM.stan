data {
  int<lower=1> N;            // num observations
  int K;
  vector[K] y[N];
  real lambda1;
  real lambda2;
}
transformed data {
  matrix[K,K] Sigma = diag_matrix(to_vector({1,1}));
}
parameters {
  vector[K] c1;                   // mean coeff
  vector[K] c2;                   // mean coeff
  vector[K] c3;                   // mean coeff
  vector[K] beta1;                   // mean coeff
  vector[K] beta2;                   // mean coeff
  vector[K] beta3;                   // mean coeff
  matrix[K,K] phi1;                  // autoregression coeff
  matrix[K,K] phi2;                  // autoregression coeff
  matrix[K,K] phi3;                  // autoregression coeff
  vector[N] z;
  real<lower=0> sigma_z;
  real<lower=-1,upper=1> phi_z;

  //cov_matrix[K] Sigma;
}
model {
  vector[K] err[N];
  c1 ~ normal(0, 1);
  c2 ~ normal(0, 1);
  c3 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  for(i in 1:K)
    for(j in 1:K) {
        phi1[i,j] ~ normal(0, 1);
        phi2[i,j] ~ normal(0, 1);
        phi3[i,j] ~ normal(0, 1);
    }
    
  for (t in 2:N) {
    if(z[t] <= lambda1)
      err[t,] = c1 + phi1 * to_vector(y[t-1,]) + beta1 * z[t-1];
    else if(z[t] > lambda1 && z[t] <= lambda2)
      err[t,] = c2 + phi2 * to_vector(y[t-1,]) + beta2 * z[t-1];
    else if(z[t] > lambda2)
      err[t,] = c3 + phi3 * to_vector(y[t-1,]) + beta3 * z[t-1];  
  }
  z[2:N] ~ normal(z[1:(N-1)]*0.85, 1);
  y[2:N] ~ multi_normal(err[2:N], Sigma);
  
}
/*
generated quantities {
  vector[N] log_lik;
  {
      vector[K] err[N];
      err[1] = to_vector(y[1,]) - mu + phi * mu;
      log_lik[1] = multi_normal_lpdf(err[1] | rep_vector(0, 2), Sigma);
      for (t in 2:N) {
        err[t] = to_vector(y[t]) - (mu + phi * to_vector(y[t-1]) + theta * err[t-1]);
        log_lik[t] = multi_normal_lpdf(err[t] | rep_vector(0, 2), Sigma);
      }
  }

}
*/
