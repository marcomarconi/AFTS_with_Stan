data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> F;
  matrix[N,K] y;
  matrix[K,F] betas;
}


parameters {
  vector<lower=0>[K] sigma;
  matrix[N,F] eff;
}

model {
  sigma ~ exponential(1);
  for(t in 1:N) {
    vector[K] mu;  
    eff[t,] ~ normal(0,1);
    for(k in 1:K){
      mu[k] = 0;
      for(f in 1:F)
        mu[k] += betas[k,f] * eff[t,f];
    }
    y[t,] ~ normal(mu, sigma);    
  }
}

