data {
  int<lower=0> N;   // # time points (equally spaced)
  vector[N] y;      // mean corrected return at time t
  int distribution;
  int P;
  int Q;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
transformed data {
  int s;
  if(P > Q)
    s = P;
  else
    s = Q;
  s += 1;  
}
parameters {
  real mu; 
  real<lower=1> nu[distribution==1];  // std log volatility time t
  real<lower = -1, upper = 1> phi[P];                  // autoregression coeff
  real<lower = -1, upper = 1> theta[Q];                // moving avg coeff
  real alpha0;                     // mean log volatility
  real<lower=-1,upper=1> alpha1;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[N] h_std;  // std log volatility time t

}
transformed parameters {
  vector[N] h = h_std * sigma;  // now h ~ normal(0, sigma)
  h[1] /= sqrt(1 - alpha1 * alpha1);  // rescale h[1]
  h += alpha0;
  for (t in 2:N)
    h[t] += alpha1 * (h[t-1] - alpha0);
}

model {
  vector[N] err;
  mu ~ normal(0, 10);
  nu ~ gamma(1, 0.1);
  sigma ~ cauchy(0, 2);
  alpha1 ~ uniform(-1, 1);
  alpha0 ~ normal(0, 2);
  h_std ~ std_normal();
  for(t in 1:(s-1)) {
    err[t] = y[t] - mu;
    for(p in 1:P)
      err[t] += phi[t] * mu;
    if(distribution == 0)  
      err[t] ~ normal(0, exp(h[t] / 2));
    else if(distribution == 1)  
      err[t] ~ student_t(nu[1], 0, exp(h[t] / 2));
  }
  for (t in s:N) {
    err[t] = y[t] - mu;
    for(p in 1:P)
      err[t] -= phi[p] * y[t-p];
    for(q in 1:Q)  
      err[t] -= theta[q] * err[t-q];
  }
  if(distribution == 0)
    err ~ normal(0, exp(h / 2));
  else
    err ~ student_t(nu[1], 0, exp(h / 2));
}
// for PSIS-LOO-LFO
generated quantities {
  vector[M+N] log_lik;
  {
    vector[N+M] h_hat;
    vector[M+N] z = append_row((y), (m));
    vector[M+N] err;
    
    for (t in 1:(s-1)) {
      if(distribution == 0)  
        log_lik[t] = normal_lpdf(z[t] | mu, exp(h[t] / 2));
      else if(distribution == 1)  
        log_lik[t] = student_t_lpdf(z[t] | nu[1], mu, exp(h[t] / 2));
      err[t] = z[t] - mu;
      for(p in 1:P)
        err[t] += phi[t] * mu;
    }
    h_hat[1] = h[1];
    for (t in 2:(N+M)) {
      if(t <= N)
        h_hat[t] = h[t];
      else  
        h_hat[t] = alpha0 + alpha1 * (h_hat[t-1] - alpha0);
    }
    for (t in s:(N+M)) {
      err[t] = z[t] - mu;
      for(p in 1:P)
        err[t] -= phi[p] * z[t-p];
      for(q in 1:Q)  
        err[t] -= theta[q] * err[t-q];
      if(distribution == 0)
        log_lik[t] = normal_lpdf(err[t] |  0, exp(h_hat[t] / 2));
      else
        log_lik[t] = student_t_lpdf(err[t] |  nu[1], 0, exp(h_hat[t] / 2));
    }
  }
}
