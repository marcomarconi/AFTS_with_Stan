data {
  int<lower=0> N;
  vector[2] y[N];
  matrix[2,2] Sigma1;
}
parameters {
  vector[2] mu;
  real<lower=-1,upper=1> rho;
  vector<lower=0>[2] alpha0;
  cov_matrix[2] Alpha1;
  cov_matrix[2] Beta1;
}
transformed parameters {
  vector[2] Xi[N];
  Xi[1][1] = Sigma1[1,1];
  Xi[1][2] = Sigma1[2,2];
  for (t in 2:N) {
    vector[2] a_1 = y[t-1] - mu;
    Xi[t] = alpha0
                    + Alpha1 * [a_1[1]*a_1[1], a_1[2]*a_1[2]]'
                    + Beta1 * Xi[t-1];
  }
}
model {
  matrix[2,2] Rho ;
  mu ~ normal(0, 1);
  alpha0 ~ normal(0, 1);
  Rho = [[1, rho], [rho, 1]];
  for (t in 2:N) 
    y[t] ~ multi_normal(mu, quad_form_diag(Rho, Xi[t]));
}
generated quantities {
  vector[N] log_lik;
  matrix[2,2] Rho = [[1, rho], [rho, 1]];
  for (t in 1:N)
    log_lik[t] = multi_normal_lpdf(y[t] | mu, quad_form_diag(Rho, Xi[t]));
}
