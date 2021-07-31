data {
  int<lower=0> N;
  int K;
  vector[K] y[N];
  matrix[K,K] Sigma1;
}
parameters {
  vector[K] mu;
  real<lower=0> A11;
  real<lower=0> A21;
  real<lower=0> A22;
  real<lower=0> arch11;
  real<lower=0> arch21;
  real<lower=0> arch22;
  real<lower=0> garch11;
  real<lower=0> garch21;
  real<lower=0> garch22;
}
transformed parameters {
  matrix[K,K] Sigma[N];
  Sigma[1] = Sigma1;
  for(t in 2:N) {
      vector[K] a = y[t-1] - mu;
      Sigma[t][1,1] = A11 + arch11 * a[1]*a[1] + garch11 * Sigma[t-1][1,1];
      Sigma[t][2,1] = A21 + arch21 * a[1]*a[2] + garch21 * Sigma[t-1][2,1];
      Sigma[t][1,2] = Sigma[t][2,1];
      Sigma[t][2,2] = A22 + arch22 * a[2]*a[2] + garch22 * Sigma[t-1][2,2];
  }
}
model {
  mu ~ normal(0,1);
  [A11, A21, A22] ~ normal(0,0.001);
  [arch11, arch21, arch22, garch11, garch21, garch22] ~ normal(0,1);
  for(t in 2:N) 
     y[t] ~ multi_normal(mu, Sigma[t]);

}

