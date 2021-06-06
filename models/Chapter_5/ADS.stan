functions {
  real geometric(int y, real phi){
    return(log(phi * (1 - phi)^y));
  }
  real ADS_loglikelihood(vector I, int A, int D, int S, real peta, real eta, real lambda_u, real lambda_d){
    return(I[1] * log(1 - peta) + 
           I[2] * (log(peta) + log(eta) + log(lambda_u) + (S-1) * log(1 - lambda_u)) +
           I[3] * (log(peta) + log(1-eta) + log(lambda_d) + (S-1) * log(1 - lambda_d))
           );
  }
}

data {
  int<lower=0> N;
  int A[N];
  int D[N];
  int S[N];
}
transformed data {
  vector[N] A_real;
  vector[N] D_real;
  vector[3] I[N] ;
  int D_count = 0;
  for(t in 1:N) {
    A_real[t] = A[t];
    D_real[t] = D[t];
    D_count = D_count + 1;
    I[t][1] = 0;
    I[t][2] = 0;
    I[t][3] = 0;
    if(A[t] == 0)
      I[t][1] = 1;
    if(A[t] == 1 && D[t] == 1)
      I[t][2] = 1;
    if(A[t] == 1 && D[t] == -1)
      I[t][3] = 1;
      
  }

}

parameters {
  real beta0 ;
  real beta1 ;
  real gamma0 ;
  real gamma1 ;
  real phi0_d ;
  real phi0_u ;
  real phi1_d ;
  real phi1_u ;
}

model {
  beta0 ~ normal(0,1);
  beta1 ~ normal(0,1);
  gamma0 ~ normal(0,1);
  gamma1 ~ normal(0,1);
  phi0_d ~ normal(0,1);
  phi0_u ~ normal(0,1);
  phi1_d ~ normal(0,1);
  phi1_u ~ normal(0,1);
  /* formulas 5.25 to 5.28
  for(t in 2:N) {
    A[t] ~ bernoulli_logit(beta0 + beta1 * A[t-1]);
    if(A[t] == 1) {
      (D[t]==1) ~ bernoulli_logit(gamma0 + gamma1 * D[t-1]);
      if(D[t] == 1)
        target += geometric(S[t]-1, inv_logit(phi0_u + phi1_u * S[t-1]));
      else
        target += geometric(S[t]-1, inv_logit(phi0_d + phi1_d * S[t-1]));
    }
  }*/
  // formula 5.29
  for(t in 2:N)
    target += ADS_loglikelihood(I[t], A[t], D[t], S[t],
        inv_logit(beta0 + beta1 * A[t-1]),  
        inv_logit(gamma0 + gamma1 * D[t-1]), 
        inv_logit(phi0_u + phi1_u * S[t-1]), 
        inv_logit(phi0_d + phi1_d * S[t-1])) ;
}
