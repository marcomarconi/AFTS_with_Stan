data {
  int<lower=0> N;
  vector[2] y[N];
  matrix[2,2] Sigma1;
}
parameters {
  vector[2] mu;
  vector<lower=0>[2] alpha0;
  cov_matrix[2] Alpha1;
  cov_matrix[2] Beta1;
  real omega0;
  real omega1;
  real omega2;
}
transformed parameters {
  vector[2] Xi[N];
  vector[N] q;
  vector[N] rho;
  Xi[1][1] = Sigma1[1,1];
  Xi[1][2] = Sigma1[2,2];
  rho[1] = 0.5;
  q[1] = omega0 / (1 - omega1 - omega2);
  for (t in 2:N) {
    vector[2] a_1 = y[t-1] - mu;
    Xi[t] = alpha0
                    + Alpha1 * [a_1[1]*a_1[1], a_1[2]*a_1[2]]'
                    + Beta1 * Xi[t-1];
    q[t] = omega0 + omega1 * rho[t-1] + omega2 * (a_1[1]*a_1[2]) / sqrt(Xi[t-1,1] * Xi[t-1,2]);
    rho[t] = (exp(q[t]) - 1.0) / (exp(q[t]) + 1.0);
  }
  
}
model {
  mu ~ normal(0, 2);
  alpha0 ~ normal(0, 2);
  [omega0, omega1, omega2] ~ normal(0, 5);
  Alpha1 ~ inv_wishart(2, diag_matrix(rep_vector(1,2)));
  Beta1 ~ inv_wishart(2, diag_matrix(rep_vector(1,2)));

  for (t in 2:N) {
    matrix[2,2] Corr = diag_matrix(rep_vector(1,2));
    Corr[1,2] = rho[t];
    Corr[2,1] = rho[t];
    y[t] ~ multi_normal(mu, quad_form_diag(Corr, Xi[t,]));
  }
}
generated quantities {
  vector[N] log_lik;
  matrix[2,2] Corr = diag_matrix(rep_vector(1,2));
  for (t in 1:N) {
    Corr[1,2] = rho[t];
    Corr[2,1] = rho[t];
    log_lik[t] = multi_normal_lpdf(y[t] | mu, quad_form_diag(Corr, Xi[t,]));
  }
}
