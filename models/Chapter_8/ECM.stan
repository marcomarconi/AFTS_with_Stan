data {
  int<lower=1> N;            
  int K;
  vector[K] y[N];
  int P;
}

parameters {
  vector[K] alpha;                   
  real beta;                  
  vector[K] c0;                  
  matrix[K,K] phi[P];                  
  cov_matrix[K] Sigma;
}
model {
  vector[K] err[N];
  c0 ~ normal(0, 10);
  beta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  for(i in 1:K)
    for(j in 1:K)
      for(p in 1:P)
        phi[p][i,j] ~ normal(0, 1);
    
  for (t in (P+2):N) {
    err[t,] = to_vector(y[t-1,]) + c0 + alpha * (to_row_vector({1,-beta}) * to_vector(y[t-1,]));
    for(p in 1:P)
      err[t,] += (phi[p] * to_vector(y[t-p,] - y[t-p-1,]));
  }
  y[(P+2):N] ~ multi_normal(err[(P+2):N], Sigma);
  
}





