
data {
  int<lower=0> N;
  vector[N] y;
  real psi1;
  real threshold;
}
parameters {
  real<lower=0> omega0_l;
  real<lower=0> omega0_g;
  real<upper=1> omega1_l;
  real<upper=1> omega1_g;
  real<upper=(1-omega1_l)> epsilon1_l;
  real<upper=(1-omega1_g)> epsilon1_g;
  real<lower=0> alpha_l;
  real<lower=0> alpha_g;
}
transformed parameters {
  vector[N] psi;
  psi[1] = psi1;
  for(t in 2:N) {
    if(y[t-1] <= threshold)
      psi[t] = omega0_l + epsilon1_l * y[t-1] + omega1_l * psi[t-1];
    else
      psi[t] = omega0_g + epsilon1_g * y[t-1] + omega1_g * psi[t-1]; 
  }
}

model {
  omega0_l ~ normal(0.5,0.5);
  omega0_g ~ normal(0.5,0.5);
  omega1_l ~ normal(0.5,0.5);
  omega0_g ~ normal(0.5,0.5);
  epsilon1_l ~ normal(0.5,0.5);
  epsilon1_g ~ normal(0.5,0.5);
  alpha_l ~ normal(1,1);
  alpha_g ~ normal(1,1);
  for(t in 2:N) {
    if(y[t-1] <= threshold)
      y[t] ~ weibull(alpha_l, psi[t]);
    else
      y[t] ~ weibull(alpha_g, psi[t]);
  }

}
generated quantities {
   vector[N] log_lik;
   {
     log_lik[1] = weibull_lpdf(y[1] | alpha_l, psi[1]);

     for (t in 2:(N)) {
       if(psi[t] < 0) {
        log_lik[t] = 0;
        continue;
       }
       if(y[t-1] <= threshold)
          log_lik[t] = weibull_lpdf(y[t] |alpha_l, psi[t]);
        else
          log_lik[t] = weibull_lpdf(y[t] |alpha_g, psi[t]);
     }
   }
}
