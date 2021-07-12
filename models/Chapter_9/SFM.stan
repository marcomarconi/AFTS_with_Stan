data {
  int<lower=0> T;
  int<lower=0> K;
  matrix[T,K] r;
  vector[T] m;
}
parameters {
  vector[K]  alpha;
  vector[K]  beta;
  corr_matrix[K] Omega;
  vector<lower=0>[K] tau;
}
transformed parameters{
  matrix[K,K] Sigma;
  Sigma = quad_form_diag(Omega, tau);
}
model {
  tau ~ exponential(1);
  Omega ~ lkj_corr(2);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 1);
  for(t in 1:T)
    r[t,] ~ multi_normal(m[t]*beta + alpha, Sigma);
}

