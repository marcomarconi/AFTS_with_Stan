data {
  int<lower=0> N;
  vector[N] y;
  int Kar;
  int family;
}

parameters {
  real mu;                       
  real<lower=0> nu[family==1];
  real<lower=0> alpha0;          
  real<lower=0,upper=1> alpha[Kar];  
}
transformed parameters {
  vector[N] sigma;
  for(i in 1:Kar)
    sigma[i] = alpha0;
  for (i in (Kar+1):N) {
    sigma[i] = alpha0;
    for(k in 1:Kar)
      sigma[i] += alpha[k] * pow(y[i-k] - mu, 2);
  }
  sigma = sqrt(sigma);
}
model {
  mu ~ normal(0, 1);
  nu ~ normal(0, 1);
  alpha0 ~ normal(1, 1);
  alpha ~ normal(0.5, 0.5);
  if(family == 0)
    y ~ normal(mu, sigma);
  else if(family == 1)  
    y ~ student_t(nu[1], mu, sigma);
}

generated quantities {
   vector[N] log_lik;
   {
     for(i in 1:N) 
      if(family == 0)
        log_lik[i] = normal_lpdf(y[i] | mu, sigma[i] );  
      else if(family == 1)    
        log_lik[i] = student_t_lpdf(y[i] | nu[1], mu, sigma[i] );  
   }
}


