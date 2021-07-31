data {
  int<lower=0> N;
  int<lower=0> K;
  vector[K] y[N];
  matrix[K,K] Sigma1;
  int approx; // use the approximated formula (if N is large)
}
parameters {
  real<lower=0,upper=1> lambda;
  vector[K] mu;
}
transformed parameters {
  matrix[K,K] Sigma[N];
  Sigma[1] = Sigma1;
  if(approx) { // approximated formula (if N is large)
    for (t in 2:N)
      Sigma[t] = (1-lambda) * ((y[t-1]-mu) * to_row_vector(y[t-1]-mu)) + lambda * Sigma[t-1];
  }else{ // recursive formula
    for (t in 2:N) {
      Sigma[t] = Sigma1 ;
      for(j in 1:(t-1)) 
        Sigma[t] = Sigma[t] + pow(lambda,(j-1)) * ((y[t-1]-mu) * to_row_vector((y[t-1]-mu))) ;
       Sigma[t] =  ((1-lambda) / (1-pow(lambda,(t-1)))) * Sigma[t];
    }
  }
      
}
model {
  lambda ~ beta(1, 1);
  for (t in 2:N)
    y[t] ~ multi_normal(mu, Sigma[t]);
}

