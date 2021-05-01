
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] x;
}

parameters {
  real<lower=0> sigma_y;
  cov_matrix[2] Sigma;
  vector[2] beta_init;
  vector[2] err_beta[N];
}
transformed parameters {
  matrix[2,N] beta;
  vector[N] mu;
  
  beta[,1] = beta_init + err_beta[1];
  mu[1] = beta[,1]' * [1, x[1]]';
  for(t in 2:N) {
    beta[,t] = beta[,t-1] + err_beta[t];
    mu[t] =  beta[,t]' * [1, x[t]]';
  }
}

model {
  beta_init ~ cauchy(0, 5);
  sigma_y ~ student_t(4, 0, 1);
  Sigma ~ inv_wishart(2, diag_matrix(rep_vector(1,2)));

  err_beta ~ multi_normal([0,0]', Sigma);
  y ~ normal(mu, sigma_y);
}
