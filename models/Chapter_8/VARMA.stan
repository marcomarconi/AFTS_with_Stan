data {
  int<lower=1> N;            // num observations
  int K;                    // Matrices size
  matrix[N,K] y;                // observed outputs

}
parameters {
  vector[K] mu;                   // mean coeff
  matrix[K,K] phi;                  // autoregression coeff
  matrix[K,K] theta;                // moving avg coeff
  cov_matrix[K] Sigma;
}
model {
  vector[K] err[N];
  mu ~ normal(0, 10);
  for(i in 1:K)
    for(j in 1:K) {
        phi[i,j] ~ normal(0, 1);
        theta[i,j] ~ normal(0, 1);
    }
  Sigma ~ inv_wishart(3, diag_matrix(rep_vector(2, K)));
  err[1] = to_vector(y[1,]) - mu + phi * mu;
  err[1] ~ multi_normal(rep_vector(0, 2), Sigma);
  for (t in 2:N) 
    err[t] = to_vector(y[t]) - (mu + phi * to_vector(y[t-1]) + theta * err[t-1]);
  err ~ multi_normal(rep_vector(0, 2), Sigma);
  
}

generated quantities {
  vector[N] log_lik;
  vector[K] err[N];
  {
      err[1] = to_vector(y[1,]) - mu + phi * mu;
      log_lik[1] = multi_normal_lpdf(err[1] | rep_vector(0, 2), Sigma);
      for (t in 2:N) {
        err[t] = to_vector(y[t]) - (mu + phi * to_vector(y[t-1]) + theta * err[t-1]);
        log_lik[t] = multi_normal_lpdf(err[t] | rep_vector(0, 2), Sigma);
      }
  }

}

