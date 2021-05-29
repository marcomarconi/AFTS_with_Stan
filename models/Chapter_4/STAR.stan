// PSIS-LOO-LFO fails here
data {
  int<lower=0> N;
  vector[N] y;
  int K;
  real sigma1;
  real gamma2_mean;
  real gamma2_sd;
}
parameters {
  real mu;
  real<lower=0> beta0;
  real<lower=0,upper=1> beta1;
  real<lower=0,upper=1> beta2;
  real gamma0;
  real gamma1;
  real<lower=0> gamma2;
}
transformed parameters {
  vector<lower=0>[N] sigma;
  for(i in 1:K)
    sigma[i] = sigma1; //beta0 + (gamma0/(1 + exp(0)));
  for (t in (K+1):N)
    sigma[t] = sqrt(
                    beta0 
                    + beta1*(y[t-1]-mu)^2 
                    + beta2*(y[t-2]-mu)^2 
                    + (gamma0 + gamma1*(y[t-1]-mu)^2)/(1 + exp(-gamma2*(y[t-1]-mu)))
                    );
}
model {
  mu ~ normal(0, 1);
  beta0 ~ normal(0.5, 0.5);
  beta1 ~ normal(0.5, 0.5);
  gamma0 ~ normal(0, 0.01); 
  gamma1 ~ normal(0, 0.1);
  gamma2 ~ normal(gamma2_mean, gamma2_sd);
  y ~ normal(mu, sigma);
}
generated quantities {
   vector[N] log_lik;

   { 
     for (t in 1:N) 
       log_lik[t] = normal_lpdf(y[t] | mu, sigma[t] );
     
   }
}
