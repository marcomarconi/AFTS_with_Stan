data {
  int<lower=0> N;
  int K;
  vector[K] y[N];
  real Q0;
  real P0;
  int likelihood;
  int prior;
}
transformed data {
  int J = floor(K*(K-1)/2);
}
parameters {
  vector[K] mu;
  real<lower=0> alpha10;
  real<lower=0> alpha11;
  real<lower=0> alpha21;
  real<lower=0> alpha22;
  real<lower=0> alpha20;
  real<lower=0> beta11;
  real<lower=0> beta21;
  real<lower=0> beta22;
  real lambda0;
  real lambda1;
  real lambda2;
}
transformed parameters {
  vector[K] Q[N];
  vector[J] P[N];
  vector[K] a[N];
  vector[K] b[N];
  Q[1] = Q0;
  P[1] = P0;
  a[1] = y[1] - mu;
  b[1][1] = a[1][1];
  b[1][2] = a[1][2] - Xi[1,3]*a[1][1] ;
  for (t in 2:N) {
    Q[t] =
    a[t] = y[t] - mu;
    b[t][1] = a[t][1];
    b[t][2] = a[t][2] - Xi[t,3]*a[t][1] ;
  }
  
}
model {
  matrix[K,K] G;
  matrix[K,K] L = diag_matrix(rep_vector(1,2));
  matrix [K,K] Sigma;
  real tmp;
  
  mu ~ normal(0, 2);
  [lambda0, lambda1, lambda2] ~ normal(0, 1);
  [alpha10, alpha11, alpha20, alpha21, alpha22] ~ normal(0, 2);
  [beta11, beta21, beta22] ~ normal(0, 1);

  for (t in 2:N) {
    G[1,1] = (Xi[t,1]);
    G[2,2] = (Xi[t,2]);
    G[2,1] = 0;
    G[1,2] = 0;
    L[2,1] = Xi[t,3];
    if(likelihood == 1)
      Sigma = L * G * L'; // multi_normal
    else if(likelihood == 2)  
      Sigma = L * sqrt(G); // multi_normal_cholesky
    if(!prior){
       if(likelihood == 1)
          y[t] ~ multi_normal(mu, Sigma);
       else if(likelihood == 2)     
          y[t] ~ multi_normal_cholesky(mu, Sigma);
       else {   // book's loglik
         tmp = 0;
         for(k in 1:K)
            tmp += log(Xi[t][k]) + pow(b[t][k],2) / Xi[t][k];
         target += -0.5 * tmp;
       }
    }
  }
}
generated quantities {/*
  vector[N] phi;
  vector[N] log_lik;
  matrix[K,K] G;
  matrix[K,K] L = diag_matrix(rep_vector(1,2));
  matrix [K,K] Sigma;
  phi[1]=0.5;
  for (t in 1:N) {
    phi[t] = (Xi[t,3] * sqrt((Xi[t,1]))) / sqrt((Xi[t,2]) +  pow(Xi[t,3],2) *  (Xi[t,1]));
    G[1,1] = (Xi[t,1]);
    G[2,2] = (Xi[t,2]);
    G[2,1] = 0;
    G[1,2] = 0;
    L[2,1] = Xi[t,3];
    //Sigma = L * G * L';
    Sigma = L * sqrt(G);
    log_lik[t] = multi_normal_cholesky_lpdf(y[t] | mu, Sigma);
  }*/
}
