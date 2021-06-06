functions {
  real geometric(real y, real phi){
    return(log(phi * (1 - phi)^y));
  }

}

data {
  int<lower=0> T; // length of data
  vector[T] log_dt; // log of time durations
  vector[T] D; // direction of movement
  int S[T]; // size of price change
  vector[T] N; // number of transactions
}

parameters {
  real  beta0 ;
  real  beta1 ;
  real  beta2 ;
  real<lower=0>  sigma ;
  real  alpha0 ;
  real  alpha1 ;
  real  gamma0 ;
  real  gamma1 ;
  real  omega0 ;
  real  omega1 ;
  real  omega2 ;
  real  Dbeta ;
  real  eta0_d ;
  real  eta1_d ;
  real  eta2_d ;
  real  eta3_d ;
  real  eta0_u ;
  real  eta1_u ;
  real  eta2_u ;
  real  eta3_u ;
}

model {
  sigma ~ cauchy(0,1);
  beta0 ~ normal(0,5);
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  gamma0 ~ normal(0,5);
  gamma1 ~ normal(0,5);
  alpha0 ~ normal(0,5);
  alpha1 ~ normal(0,5);
  omega0 ~ normal(0,5);
  omega1 ~ normal(0,5);
  omega2 ~ normal(0,5);
  Dbeta ~ normal(0,1);
  eta0_d ~ normal(0,5);
  eta1_d ~ normal(0,5);
  eta2_d ~ normal(0,5);
  eta3_d ~ normal(0,5);
  eta0_u ~ normal(0,5);
  eta1_u ~ normal(0,5);
  eta2_u ~ normal(0,5);
  eta3_u ~ normal(0,5);
  
  for(t in 5:T) {
    // formula (5.48)
    log_dt[t] ~ normal(beta0 + beta1*log_dt[t-1] + beta2*S[t-1], sigma);
    // formula (5.49) and (5.50)
    if (N[t] == 0)
      target += bernoulli_logit_lpmf(0 | alpha0 + alpha1*log_dt[t]);
    else
      target += bernoulli_logit_lpmf(1 | alpha0 + alpha1*log_dt[t])
                              + geometric(N[t]-1, inv_logit(gamma0 + gamma1*log_dt[t]));
    // we use the comulative normal to represent formula (5.51), nice trick!
    (D[t]==1) ~ bernoulli(normal_cdf(0, omega0 + omega1*D[t-1] + omega2*log_dt[t], exp(Dbeta * fabs(sum(D[(t-4):(t-1)])))));  
    // formulas (5.52) and (5.53)
    if(D[t] == -1)
      (S[t]-1) ~ poisson(exp( eta0_d + eta1_d*log(N[t]+1) + eta2_d*log_dt[t] + eta3_d*S[t-1]));
    else
      (S[t]-1) ~ poisson(exp( eta0_u + eta1_u*log(N[t]+1) + eta2_u*log_dt[t] + eta3_u*S[t-1]));
  }
}

