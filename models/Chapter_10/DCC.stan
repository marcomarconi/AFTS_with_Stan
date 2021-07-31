data {
  int<lower=0> N;
  int K;
  vector[K] y[N];
  matrix[K,K] Rho_hat;
  matrix[K,K] Psi[N];
  int prior;
}
parameters {
  real<lower=0> nu;
  vector<lower=1e-12>[K] alpha0;
  vector<lower=0,upper=1>[K] alpha1;
  vector<lower=0,upper=1>[K] alpha2;
  real<lower=0,upper=1> theta1;
  real<lower=0,upper=(1-theta1)> theta2;
}

transformed parameters {
  matrix[K,K] D[N];
  matrix[K,K] Rho_m[N];
  matrix[N,K] a;

  matrix[K,K]  Alpha0 = diag_matrix(alpha0);
  matrix[K,K]  Alpha1 = diag_matrix(alpha1);
  matrix[K,K]  Alpha2 = diag_matrix(alpha2);
  
  D[1] = diag_matrix([1, 1, 1, 1]');
  Rho_m[1] = diag_matrix([1, 1, 1, 1]');
  a[1,]  = (y[1])';
  for (t in 2:N) {
    a[t,] = (y[t])';
    D[t] = diag_matrix([1, 1, 1, 1]');
    for(k in 1:K)
      D[t,k,k] = sqrt(alpha0[k] + alpha1[k] * D[t-1,k,k]^2 + alpha2[k] * (a[t-1,k] * a[t-1,k]));
    Rho_m[t] = (1-theta2-theta1) * Rho_hat + theta1 * Psi[t-1] + theta2 * Rho_m[t-1];
  }
  
}
model {
  matrix [K,K] Sigma;
  nu ~ normal(0, 10);
  alpha0 ~ normal(0,.1);
  alpha1 ~ normal(0,1);
  alpha2 ~ normal(0,.1);
  [theta1, theta2] ~ normal(0, 1);

  for (t in 1:N) {
    Sigma = D[t] * Rho_m[t] * D[t];
    if(!prior) 
      y[t] ~ multi_student_t(nu, rep_vector(0, K), Sigma);
    
  }
}
generated quantities {
  vector[N] log_lik;
  for (t in 1:N)  {
    matrix [K,K] Sigma = D[t] * Rho_m[t] * D[t];
    log_lik[t] = multi_student_t_lpdf(y[t] | nu, rep_vector(0, K), Sigma);
  }
}
