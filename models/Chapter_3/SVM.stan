data {
  int<lower=0> N;   // # time points (equally spaced)
  vector[N] y;      // mean corrected return at time t
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
    vector[N] h_std;  // std log volatility time t

}
transformed parameters {
  vector[N] h = h_std * sigma;  // now h ~ normal(0, sigma)
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  for (t in 2:N)
    h[t] += phi * (h[t-1] - mu);
}

model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  mu ~ cauchy(0, 10);
  h_std ~ std_normal();
  y ~ normal(0, exp(h / 2));
}
// for PSIS-LOO-LFO
generated quantities {
  vector[M+N] log_lik;
  {
    vector[N+M] h_hat;
    h_hat = append_row(h, rep_vector(0, M));
    for (t in 1:N) 
        log_lik[t] = normal_lpdf(y[t] | 0, exp(h_hat[t] / 2));
    for (t in (N+1):(N+M)) {
        h_hat[t] = mu + phi * (h_hat[t-1] - mu);
        log_lik[t] = normal_lpdf(m[t-N] |  0, exp(h_hat[t] / 2));
    }
  }
}
