data {
  int<lower=0> N;
  vector[2] y[N];
  vector[2+1] Xi0;
  int likelihood;
  int prior;
}
parameters {
  vector[2] mu;
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
  vector[2+1] Xi[N];
  vector[2] a[N];
  vector[2] b[N];
  Xi[1,] = Xi0;
  a[1] = y[1] - mu;
  b[1][1] = a[1][1];
  b[1][2] = a[1][2] - Xi[1,3]*a[1][1] ;
  // equation (10.28)
  for (t in 2:N) {
    Xi[t,1] = (alpha10 + alpha11 * pow(b[t-1][1], 2) + beta11 * Xi[t-1,1]); // g11
    Xi[t,2] = (alpha20 + alpha21 * pow(b[t-1][1], 2) + alpha22 * pow(b[t-1][2], 2) + beta21 * Xi[t-1,1] + beta22 * Xi[t-1,2]); // g22
    Xi[t,3] = lambda0 + lambda1 *  Xi[t-1,3] + lambda2 * a[t-1][2]; // q21
    a[t] = y[t] - mu;
    b[t][1] = a[t][1];
    b[t][2] = a[t][2] - Xi[t,3]*a[t][1] ;
  }
  
}
model {
  matrix[2,2] G;
  matrix[2,2] L = diag_matrix(rep_vector(1,2));
  matrix [2,2] Sigma;

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
    if(likelihood == 1) // multi_normal
      Sigma = L * G * L'; 
    else if(likelihood == 2)  // multi_normal_cholesky
      Sigma = L * sqrt(G); 
    if(!prior){
       if(likelihood == 1) // multi_normal
          y[t] ~ multi_normal(mu, Sigma);
       else if(likelihood == 2)     // multi_normal_cholesky
          y[t] ~ multi_normal_cholesky(mu, Sigma); 
       else if(likelihood == 0)  {   // book's loglik
         real tmp = 0;
         for(k in 1:2)
            tmp += log(Xi[t][k]) + pow(b[t][k],2) / Xi[t][k];
         target += -0.5 * tmp;
       }
    }
  }
}
generated quantities {
  vector[N] rho;
  vector[N] log_lik;
  matrix[2,2] G;
  matrix[2,2] L = diag_matrix(rep_vector(1,2));
  matrix [2,2] Sigma;
  rho[1]=0.5;
  for (t in 1:N) {
    rho[t] = (Xi[t,3] * sqrt((Xi[t,1]))) / sqrt((Xi[t,2]) +  pow(Xi[t,3],2) *  (Xi[t,1])); // equation (10.30)
    G[1,1] = (Xi[t,1]);
    G[2,2] = (Xi[t,2]);
    G[2,1] = 0;
    G[1,2] = 0;
    L[2,1] = Xi[t,3];
    // we always use multi_normal_cholesky here
    Sigma = L * sqrt(G);
    log_lik[t] = multi_normal_cholesky_lpdf(y[t] | mu, Sigma);
  }
}
