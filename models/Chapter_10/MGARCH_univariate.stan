
// first example in 10.7 APPLICATION,  Univariate Models
data {
  int<lower=0> N;
  int P;
  int K;
  vector[K] y[N];
  vector[K] sigma1;
  int<lower=0> M; 
  vector[M] m; // for PSIS-LOO-LFO
}
parameters {
  vector[K] mu;
  vector<lower=-1,upper=1>[K] ar[P];
  vector<lower=0>[K] alpha0;
  vector<lower=0,upper=1>[K] alpha1;
  vector<lower=0,upper=1>[K] beta1;
}
transformed parameters {
  vector<lower=0>[K] sigma[N];
  for(k in 1:K)
    sigma[1][k] = sigma1[k];
  for (t in 2:N)
    for(k in 1:K)
      sigma[t][k] = sqrt(alpha0[k]
                    + alpha1[k] * pow(y[t-1][k] - mu[k], 2)
                    + beta1[k] * pow(sigma[t-1][k], 2));
}
model {
  mu ~ normal(0, 2);
  //ar ~ normal(0, 1);
  alpha0 ~ cauchy(0, 1);
  alpha1 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  for (n in (P+1):N) {
    vector[K] tmp = mu;
    for(k in 1:K){
      for (p in 1:P)
        tmp[k] += ar[p][k] * y[n-p][k];
    }
    y[n] ~ normal(mu, sigma[n]);  
  }
}
generated quantities {
  vector[N] log_lik;
  for (t in 1:N) 
      log_lik[t] = normal_lpdf(y[t] | mu, sigma[t]);
  
}

