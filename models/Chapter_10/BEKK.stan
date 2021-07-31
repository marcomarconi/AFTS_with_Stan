data {
  int<lower=0> N;
  int K;
  vector[K] y[N];
  matrix[K,K] Sigma1;


}
parameters {
  vector[K] mu;
  cholesky_factor_cov[K] A;
  matrix[K,K] ARCH;
  matrix[K,K] GARCH;
}
transformed parameters {
  matrix[K,K] Sigma[N];
  Sigma[1] = Sigma1;
  for(t in 2:N) {
      vector[K] a = y[t-1] - mu;
      Sigma[t] =  multiply_lower_tri_self_transpose(A) + quad_form(a * to_row_vector(a), ARCH') + quad_form(Sigma[t-1], GARCH');
  }
}
model {
  mu ~ normal(0,1);
  for(k in 1:K)
    for(j in 1:K)
      [A[k,j], ARCH[k,j], GARCH[k,j]] ~ normal(0,1);
  for(t in 2:N) 
     y[t] ~ multi_normal(mu, Sigma[t]);

}
generated quantities {
}
