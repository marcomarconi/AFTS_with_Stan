data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] y;
  int p;
}
parameters {
  vector[K] mu;
  matrix[K,K] phi[p];
  cov_matrix[K] Sigma;
}
model {
  mu ~ normal(0,1);
  Sigma ~ inv_wishart(2, diag_matrix(rep_vector(1, K)));
  for(t in (p+1):N) {
    vector[K] mu_tmp = to_vector(mu);
    for(i in 1:p)
      mu_tmp += phi[p] * to_vector(y[t-p,]);  
    y[t,] ~ multi_normal(mu_tmp, Sigma);
  }
}
generated quantities {
  vector[N] log_lik;
  {
    for(i in 1:p)
      log_lik[i] = multi_normal_lpdf(y[i,] | mu , Sigma);
    for (t in (p+1):(N)) {
      vector[K] mu_tmp = mu;
      for(i in 1:p)
        mu_tmp += phi[i] * to_vector(y[t-i,]);
      log_lik[t] = multi_normal_lpdf(y[t,] | mu_tmp, Sigma);
    }
  }
}
