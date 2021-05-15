data {
  int<lower=0> N;
  vector[N] y;
  vector[N] e;
  real sigma1;
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
}
transformed parameters {
  real<lower=0> sigma[N];
  sigma[1] = sigma1;
  for (t in 2:N)
    sigma[t] = sqrt(alpha0
                    + alpha1 * e[t-1]
                    + beta1 * pow(sigma[t-1], 2));
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
    vector[N] log_lik;
   {
     for (t in 1:N) 
       log_lik[t] = normal_lpdf(y[t] | mu, sigma[t] );
     
   }
}
