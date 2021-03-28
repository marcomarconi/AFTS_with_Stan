data {
  int<lower=1> N;            // num observations
    vector[N] y;                 // observed outputs
    int M;
    vector[M] m;
}
parameters {
  real mu;                   // mean coeff
  real<lower = -1, upper = 1> phi;                  // autoregression coeff
  real<lower = -1, upper = 1> theta;                // moving avg coeff
  real<lower=0> sigma;       // noise scale
}
model {
  real err;
  mu ~ normal(0, 10);
  phi ~ normal(0, 2);
  theta ~ normal(0, 2);
  sigma ~ cauchy(0, 5);
  err = y[1] - mu + phi * mu;
  err ~ normal(0, sigma);
  for (t in 2:N) {
    err = y[t] - (mu  + phi * y[t-1] + theta * err);
    err ~ normal(0, sigma);
  }
}
generated quantities {
  vector[M+N] r_hat;
  vector[M+N] log_lik;

  r_hat[1] = normal_rng(mu, sigma);
  log_lik[1] = normal_lpdf(y[1] | mu, sigma);
  
  {
      vector[M+N] z = append_row((y), (m));
      real err;
      err = y[1] - mu;
      for (t in 2:(M+N)) {
        real nu = (mu + phi * z[t-1] + theta * err);
        r_hat[t] = normal_rng(nu, sigma);
        log_lik[t] = normal_lpdf(z[t] | nu, sigma);
        err = z[t] - nu;
      }
  }

}