data {
  int<lower=0> T;                // number of time points
  vector[T] r;                     // return at time t
  int<lower=0> m;
  real sigma1;
}

parameters {
  real mu;                       // average return
  real<lower=0> sigma_eta;          // noise intercept
  cov_matrix[m] Cov;
}
transformed parameters {
  vector[T] sigma_squared;
  vector[m] mus;
  for(i in 1:m)
    mus[i] = mu;
  for(t in 1:m)
    sigma_squared[t] = sigma1;
  for (t in (m+1):T) {
    sigma_squared[t] = pow(sigma_eta, 2) + to_row_vector((r[(t-m):(t-1)]-mus))  *  Cov  *  (r[(t-m):(t-1)]-mus);
  }

}
model {
  mu ~ normal(0, 1);
  sigma_eta ~ normal(0, 1);
  //Cov ~ lkj_corr(1);
  r[(m+1):T] ~ normal(mu, sqrt(sigma_squared[(m+1):T]));

}

generated quantities {


}


