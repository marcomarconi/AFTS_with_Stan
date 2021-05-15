data {
  int<lower=0> N;
  vector[N] y;
  real sigma1;
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1-alpha1> beta1;
}
transformed parameters {
  vector<lower=0>[N] sigma;
  sigma[1] = sigma1;
  for (t in 2:N)
    sigma[t] = sqrt(alpha0
                    + alpha1 * square(y[t-1] - mu)
                    + beta1 * square(sigma[t-1]));
  
}
model {
  alpha0 ~ cauchy(0, 1);
  y ~ normal(mu, sigma);
}
generated quantities {
   vector[N] log_lik;
   {
     for (t in 1:N) 
       log_lik[t] = normal_lpdf(y[t] | mu, sigma[t] );
     
   }
}
