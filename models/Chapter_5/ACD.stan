
data {
  int<lower=0> N;
  vector[N] y;
  int family;
  real psi1;
}
parameters {
  real<lower=0> omega0;
  real<upper=1> omega1;
  real<upper=(1-omega1)> epsilon1;
  real<lower=0> alpha;
  vector<lower=0>[family==1] ka;
}
transformed parameters {
  vector[N] psi;
  psi[1] = psi1;
  for(t in 2:N)
    psi[t] = omega0 + epsilon1 * y[t-1] + omega1 * psi[t-1];
}

model {
  omega0 ~ normal(0.5,0.5);
  omega1 ~ normal(0.5,0.5);
  epsilon1 ~ normal(0.5,0.5);
  alpha ~ normal(1,1);
  ka ~ normal(1,1);
  if(family==0)
    y ~ weibull(alpha, psi);
  else
    y ~ gamma(alpha, 1. ./ (ka[1]*psi));
}
generated quantities {
   vector[N] log_lik;
   {
    if(family==0)
      log_lik[1] = weibull_lpdf(y[1] | alpha, psi[1]);
    else
      log_lik[1] = gamma_lpdf(y[1] | alpha, 1. ./ (ka[1]*psi[1]));
    
     for (t in 2:(N)) {
       if(psi[t] < 0) {
        log_lik[t] = 0;
        continue;
       }
       if(family==0)
        log_lik[t] = weibull_lpdf(y[t] | alpha, psi[t]);
       else
        log_lik[t] = gamma_lpdf(y[t] | alpha, 1. ./ (ka[1]*psi[t]));
     }
   }
}
