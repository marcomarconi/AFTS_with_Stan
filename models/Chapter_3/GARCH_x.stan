data {
  int<lower=0> N;
  real y[N];
  real x[N];
  int K;
}
transformed data {
  vector[N] x_std;
  real m = mean(x);
  real stdev = sd(x);
  x_std = (to_vector(x) - m) / stdev; // dividing by stdev avoids divergencies
}
parameters {
  real mu;
  real<lower=0> sigma1;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> beta1;
  real gamma;
}
transformed parameters {
  real<lower=0> sigma[N];
  for(i in 1:K)
    sigma[i] = sigma1;
  for (t in (K+1):N)
    sigma[t] = sqrt(
                      alpha0
                    + alpha1 * pow(y[t-K] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2)
                    + gamma * pow(x_std[t-1], 2)
                    )                    ;
}
model {
  mu ~ normal(0, 1);
  sigma1 ~ cauchy(0, 1);
  gamma ~ normal(0, 0.01);
  y ~ normal(mu, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (t in 1:N)
      log_lik[t] = normal_lpdf(y[t] | mu, sigma[t]);
  
}
