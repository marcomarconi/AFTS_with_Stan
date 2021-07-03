
data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] y;
  int p;
}
transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(append_row(rep_row_vector(0, 2), y[1:(N-1),])) * sqrt(N - 1);
  R_ast = qr_thin_R(append_row(rep_row_vector(0, 2), y[1:(N-1),])) / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
}

parameters {
  vector[K] mu;
  cov_matrix[K] Sigma;
  matrix[K,K] theta[p];     

}

model {
  mu ~ normal(0,1);
  Sigma ~ inv_wishart(2, diag_matrix(rep_vector(1, K)));
  for(t in (p+1):N) {
    row_vector[K] mu_tmp = to_row_vector(mu);
    for(i in 1:p)
      mu_tmp += Q_ast[t+1-i,] * theta[i]; // QR
    y[t,] ~ multi_normal(mu_tmp, Sigma);
  }
}

generated quantities {
  matrix[K,K] phi[p];
  vector[N] log_lik;

  for(i in 1:p)
    phi[i] = R_ast_inverse * theta[i]; 
    
  {
    for(i in 1:p)
      log_lik[i] = multi_normal_lpdf(y[i,] | mu , Sigma);
    for (t in (p+1):(N)) {
      row_vector[K] mu_tmp = to_row_vector(mu);
      for(i in 1:p)
        mu_tmp += Q_ast[t+1-i,] * theta[i];
      log_lik[t] = multi_normal_lpdf(y[t,] | mu_tmp, Sigma);
    }
  }
}
