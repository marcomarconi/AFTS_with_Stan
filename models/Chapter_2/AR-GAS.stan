data {
  int<lower=0> N;
  vector[N] y;

}
parameters {
  real omega;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] ar;
  real err[N];
  ar[1] = 1;
  ar[2] = 1;
  err[1] = y[1];
  err[2] = y[2] - 2*(inv_logit(ar[2])-0.5) * y[1];
  for (t in 3:N) {
    ar[t] = (omega
                    + alpha * (err[t-1] * y[t-2]) / sigma^2
                    + beta * ar[t-1]);
 
    err[t] = y[t] - 2*(inv_logit(ar[t])-0.5) * y[t-1];
  }
}
model {
  omega ~ normal(0, 0.1);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  y[3:N] ~ normal(2*(inv_logit(ar[3:N])-0.5) .* y[2:(N-1)], sigma);
}

