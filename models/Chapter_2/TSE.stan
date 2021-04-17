data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real beta;
  real<lower=-1, upper=1> theta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] a;    // error term at time t
  a[1] = (y[1] - beta*x[1]);
  for (t in 2:N) 
    a[t] = (y[t] - beta*x[t]) - theta * a[t - 1];

}
model {
  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  target += normal_lpdf(a | 0, sigma);
}

