data {
    int<lower=1> N;            // num observations
    vector[N] y;                 // observed outputs
    int<lower=0> E; 
    real e[E];
    int P;
    int Q;
    real sigma1;
    int distribution;
    int<lower = 0, upper = 1> garch_m;
    int<lower = 0, upper = 1> egarch;
    int<lower = 0, upper = 1> tgarch;
    int M;
    vector[M] m;
}
transformed data {
  int s = 0;
  if(P > Q)
    s = P;
  else
    s = Q;
  s += 1;  

}
parameters {
  real mu;                   // mean coeff
  real<lower = -1, upper = 1> phi[P];                  // autoregression coeff
  real<lower = -1, upper = 1> theta[Q];                // moving avg coeff
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=1> alpha2[E>0];
  real<lower=0,upper=(1-alpha1)> beta1;
  real premium[garch_m];
  real gamma[egarch];
  real<lower=0> lambda[tgarch];
  real<lower=0> nu[distribution==1];
  real<lower=0> rate[distribution==3];
}
transformed parameters {
  vector[N] err;
  real<lower=0> sigma[N];

  for(t in 1:(s-1)) {
    err[t] = y[t] - mu;
    if(garch_m == 1)
      err[t] -= premium[1]*(sigma[t])^2;
    for(p in 1:P)
      err[t] += phi[t] * mu;
  }
  for(t in s:N) {
    err[t] = y[t] - mu;
    if(garch_m == 1)
      err[t] -= premium[1]*(sigma[t])^2;
    for(p in 1:P)
      err[t] -= phi[p] * y[t-p];
    for(q in 1:Q)  
      err[t] -= theta[q] * err[t-q];
  }

  sigma[1] = sigma1;
  for (t in 2:N) {
    if(egarch)
      sigma[t] = sqrt(exp(alpha0
                    + alpha1 * ((fabs(err[t-1]) + gamma[1]*err[t-1]) / sigma[t-1] )
                    + beta1 * log(sigma[t-1]^2)));
    else if(tgarch)
      sigma[t] = sqrt(alpha0
                    + ((alpha1 + lambda[1]*(err[t-1] < 0)) * pow(err[t-1], 2))
                    + beta1 * pow(sigma[t-1], 2));
    else if(E > 0)
      sigma[t] = sqrt(alpha0
                    + alpha1 * err[t-1]^2
                    + alpha2[1] * e[t-1]^2
                    + beta1 * sigma[t-1]^2);
    else
      sigma[t] = sqrt(alpha0
                    + alpha1 * err[t-1]^2
                    + beta1 * sigma[t-1]^2);
  }
}
model {
  mu ~ normal(0, 10);
  nu ~ gamma(1, 0.1);
  rate ~ exponential(1);
  gamma ~ normal(0, 1);
  lambda ~ normal(0, 1);
  premium ~ normal(0, 2);
  phi ~ uniform(-1, 1);
  theta ~ uniform(-1, 1);
  alpha0 ~ normal(0, 1);
  
  if(distribution == 0)  
      target += normal_lpdf(err | 0, sigma);
  else if(distribution == 1)  
      target += student_t_lpdf(err | nu[1], 0, sigma);
  else if(distribution == 2)  
      target += double_exponential_lpdf(err | 0, sigma);
  else if(distribution == 3)  
      target += exp_mod_normal_lpdf(err |  0, sigma, rate[1]);
}
generated quantities {
  vector[M+N] log_lik;
  {
    vector[M+N] err_hat;
    real sigma_hat[N+M];
    vector[M+N] z = append_row((y), (m));

    for(t in 1:(N+M)){
      if(t <= N) {
        err_hat[t] = err[t];
        continue;
      }
      err_hat[t] = z[t] - mu;
      if(garch_m == 1)
        err_hat[t] -= premium[1]*(sigma_hat[t])^2;
      for(p in 1:P)
        err_hat[t] -= phi[p] * z[t-p];
      for(q in 1:Q)  
        err_hat[t] -= theta[q] * err_hat[t-q];
    }
      
    for (t in 1:(N+M)){
      if(t <= N) {
        sigma_hat[t] = sigma[t];
        continue;
      }
      if(egarch)
        sigma_hat[t] = sqrt(exp(alpha0
                    + alpha1 * ((fabs(err_hat[t-1]) + gamma[1]*err_hat[t-1]) / sigma_hat[t-1] )
                    + beta1 * log(sigma_hat[t-1]^2)));
      else if(tgarch)
        sigma_hat[t] = sqrt(alpha0
                    + ((alpha1 + lambda[1]*(err_hat[t-1] < 0)) * err_hat[t-1]^2)
                    + beta1 * sigma_hat[t-1]^2);    
      else if(E > 0)
        sigma_hat[t] = sqrt(alpha0
                    + alpha1 * err_hat[t-1]^2
                    + alpha2[1] * e[t-1]^2
                    + beta1 * sigma_hat[t-1]^2);                 
      else
        sigma_hat[t] = sqrt(alpha0
                    + alpha1 * err_hat[t-1]^2
                    + beta1 * sigma_hat[t-1]^2);  
    }
    

    for (t in 1:(N+M)) {
      if(distribution == 0)  
        log_lik[t] = normal_lpdf(err_hat[t] | 0, sigma_hat[t]); 
      else if(distribution == 1)  
        log_lik[t] = student_t_lpdf(err_hat[t] | nu[1], 0, sigma_hat[t]);
      else if(distribution == 2)    
        log_lik[t] = double_exponential_lpdf(err_hat[t] | 0, sigma_hat[t]);
      else if(distribution == 3)  
        log_lik[t] = exp_mod_normal_lpdf(err_hat[t] | 0, sigma_hat[t], rate[1]);
    }
  }
}
