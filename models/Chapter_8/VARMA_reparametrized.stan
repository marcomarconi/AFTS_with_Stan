data {
  int<lower=1> N;            // num observations
  int K;                  // Matrices size 
  matrix[N,K] y;                // observed outputs

}

transformed data {
  matrix[K,N] y_vector; 
  y_vector = y';
}

parameters {
  vector[K] mu;                   // mean coeff
  matrix[K,K] phi;                  // autoregression coeff
  matrix[K,K] theta;                // moving avg coeff
  corr_matrix[K] Omega;
  vector<lower=0>[K] tau;
}

transformed parameters {
  vector[K] epsilon[N];
  matrix[K,K] Sigma;
  Sigma = quad_form_diag(Omega, tau);
  epsilon[1] = y_vector[,1] - mu + phi * mu;
  for (t in 2:N) 
    epsilon[t] = y_vector[,t] - (mu + phi * y_vector[,t-1] + theta * epsilon[t-1]);
}
model {
  mu ~ normal(0, 10);
  for(i in 1:K)
    for(j in 1:K) {
        phi[i,j] ~ normal(0, 1);
        theta[i,j] ~ normal(0, 1);
    }
  tau ~ exponential(1);
  Omega ~ lkj_corr(2);
  epsilon ~ multi_normal(rep_vector(0, 2), Sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[K] err[N];
  {
      err[1] = y_vector[,1] - mu + phi * mu;
      log_lik[1] = multi_normal_lpdf(err[1] | rep_vector(0, 2), Sigma);
      for (t in 2:N) {
        err[t] = y_vector[,t] - (mu + phi * y_vector[,t-1] + theta * err[t-1]);
        log_lik[t] = multi_normal_lpdf(err[t] | rep_vector(0, 2), Sigma);
      }
  }

}

