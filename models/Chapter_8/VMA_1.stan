
data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] y;
}
transformed data {
  vector[K] y_vector[N];
  for(i in 1:N)
    y_vector[i] = to_vector(y[i,]);
}

parameters {
  vector[K] mu;
  cov_matrix[K] Sigma;
  matrix[K,K] A;     

}
transformed parameters {
  matrix[K,N] epsilon;       // error terms
  epsilon[,1] = y_vector[1] - mu;
  for (t in 2:N) {
    epsilon[,t] = y_vector[t] - mu  - A * epsilon[,t - 1];
  }
}
model {
  vector[K] eta[N];
  mu ~ normal(0,5);
  Sigma ~ inv_wishart(3, diag_matrix(rep_vector(2, K)));
  for(i in 1:K)
    for(j in 1:K)
        A[i,j] ~ normal(0, 1);
  for(t in 2:N) 
    eta[t] = (mu +  A *  epsilon[,t - 1]);
   y_vector[2:N] ~ multi_normal(eta[2:N], Sigma);
  
}
generated quantities {
  vector[N] log_lik;
  {
    log_lik[1] = multi_normal_lpdf(y[1,] | mu , Sigma);
    for (t in 2:N) {
      vector[K] mu_tmp = mu;
      mu_tmp += A * to_vector(y[t-1,]);
      log_lik[t] = multi_normal_lpdf(y[t,] | mu_tmp, Sigma);
    }
  }
}
