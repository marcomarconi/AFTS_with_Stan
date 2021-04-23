data {
  int<lower=0> N;
  vector[N] y;
  int Kar;
  real sigma1;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
}
transformed parameters {
  real<lower=0> sigma[N];
  for(i in 1:Kar)
    sigma[i] = sigma1;
  for (t in (Kar+1):N)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-Kar] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));
}
model {
  y ~ normal(mu, sigma);
}
generated quantities {
   vector[M+N] log_lik;
   vector[M+N] y_hat;
   {
     vector[M+N] z = append_row(y, m);
     vector[N+M] sigma_hat = append_row(to_vector(sigma), rep_vector(0, M));
     for (t in 1:N) {
        y_hat[t] = normal_rng(mu, sigma[t]);
        log_lik[t] = normal_lpdf(z[t] | mu, sigma[t]);
     }
     for (t in (N+1):(N+M)) {
        sigma_hat[t] = sqrt(alpha0
                    + alpha1 * pow(z[t-Kar] - mu, 2)
                    + beta1 * pow(sigma_hat[t-1], 2));
       y_hat[t] = normal_rng(mu, sigma_hat[t] );
       log_lik[t] = normal_lpdf(z[t] | mu, sigma_hat[t] );
     }
   }
}
