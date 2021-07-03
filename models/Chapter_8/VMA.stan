
data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] y;
  int Q;
}

transformed data {
  vector[K] y_vector[N];
  for(i in 1:N)
    y_vector[i] = to_vector(y[i,]);
}

parameters {
  vector[K] mu;
  cov_matrix[K] Sigma;
  matrix[K,K] Theta[Q];     
}

transformed parameters {
  matrix[K,N] epsilon;       // error terms
  for (t in 1:N) {
    epsilon[,t] = y_vector[t] - mu;
    for(q in 1:min(t-1, Q))
      epsilon[,t] = epsilon[,t] -  Theta[q] * epsilon[,t-q];
  }
}
model {
  vector[K] eta[N];
  mu ~ normal(0,5);
  Sigma ~ inv_wishart(3, diag_matrix(rep_vector(2, K)));
  for(k in 1:Q)
    for(i in 1:K)
      for(j in 1:K)
        Theta[k][i,j] ~ normal(0, 1);
  for(t in 1:N) {
    eta[t] = mu;
    for(q in 1:min(t-1, Q))
      eta[t] += Theta[q] * epsilon[,t - q];  
  }
  y_vector ~ multi_normal(eta, Sigma);
}

generated quantities {
  vector[N] log_lik;
  {
    for(i in 1:Q)
      log_lik[i] = multi_normal_lpdf(y[i,] | mu , Sigma);
    for (t in (Q+1):N) {
      vector[K] mu_tmp = mu;
      for(i in 1:Q)
        mu_tmp += Theta[i] * to_vector(y[t-i,]);
      log_lik[t] = multi_normal_lpdf(y[t,] | mu_tmp, Sigma);
    }
  }
}
