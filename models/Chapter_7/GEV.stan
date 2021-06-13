functions {
  real GEVd_lpdf(real y, real mu, real sigma, real xi) {
    //  GEV log pdf 
    real t;
    real log_t;
    if(xi == 0) 
      t = exp(-(y-mu)/sigma);
    else 
      t = pow((1 + xi*((y-mu)/sigma)), -1/xi);
    return log(1/sigma * pow(t, xi+1) * exp(-t) );
  }
}

data {
  int N;
  vector[N] y;
}

parameters {
  real mu;
  real xi;
  real<lower=0> sigma;
}

model {
  mu ~ normal(0,10);
  xi ~ normal(0,1);
  sigma ~ cauchy(0,1);
  for(i in 1:N)
    y[i] ~ GEVd(mu, sigma, xi);
}

