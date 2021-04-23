data {
  int<lower=0> N;
  vector[N] y;
  real sigma1;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
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
   vector[M+N] log_lik;
   vector[M+N] y_hat;
   {
     vector[M+N] z = append_row(y, m);
     vector[N+M] sigma_hat = append_row(sigma, rep_vector(0, M));
     y_hat[1] = normal_rng(mu, sigma1 );
     log_lik[1] = normal_lpdf(z[1] | mu, sigma1 );
     for (t in 2:(N+M)) {
        sigma_hat[t] = sqrt(alpha0
                    + alpha1 * pow(z[t-1] - mu, 2)
                    + beta1 * pow(sigma_hat[t-1], 2));
       y_hat[t] = normal_rng(mu, sigma_hat[t] );
       log_lik[t] = normal_lpdf(z[t] | mu, sigma_hat[t] );
     }
   }
}
