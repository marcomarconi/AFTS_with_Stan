data {
  int N;
  vector[N] y;
  int Kar ;
}
parameters {
  vector<lower = -1, upper = 1>[Kar] phi[2];
  real<lower=0> c1;
  real<upper=0> c2;
  vector<lower = 0>[2] sigma;
  vector<lower = 0, upper = 1>[2] p;
  real<lower = 0, upper = 1> xi1_init; 
}
transformed parameters {
  vector[2] c;
  matrix[N, 2] eta;
  matrix[N, 2] xi;
  vector[N] f;
  c[1] = c1;
  c[2] = c2;
  // fill in etas
  for(t in 1:N) {
    if(t > Kar) {
      for(k in 1:2)
        eta[t,k] = exp(normal_lpdf(y[t]| c[k] + phi[k] .* y[(t-Kar):(t-1)], sigma[k]));
    } else {
      for(k in 1:2)
        eta[t,k] = exp(normal_lpdf(y[t]| c[k] , sigma[k]));
    }
  }
  
  // work out likelihood contributions
  
  for(t in 1:N) {
    // for the first observation
    if(t==1) {
      f[t] = p[1]*xi1_init*eta[t,1] + // stay in state 1
             (1 - p[1])*xi1_init*eta[t,2] + // transition from 1 to 2
             p[2]*(1 - xi1_init)*eta[t,2] + // stay in state 2 
             (1 - p[2])*(1 - xi1_init)*eta[t,1]; // transition from 2 to 1
      
      xi[t,1] = (p[1]*xi1_init*eta[t,1] +(1 - p[2])*(1 - xi1_init)*eta[t,1])/f[t];
      xi[t,2] = 1.0 - xi[t,1];
    
    } else {
    // and for the rest
      
      f[t] = p[1]*xi[t-1,1]*eta[t,1] + // stay in state 1
             (1 - p[1])*xi[t-1,1]*eta[t,2] + // transition from 1 to 2
             p[2]*xi[t-1,2]*eta[t,2] + // stay in state 2 
             (1 - p[2])*xi[t-1,2]*eta[t,1]; // transition from 2 to 1
      
      // work out xi
      
      xi[t,1] = (p[1]*xi[t-1,1]*eta[t,1] +(1 - p[2])*xi[t-1,2]*eta[t,1])/f[t];
      
      // there are only two states so the probability of the other state is 1 - prob of the first
      xi[t,2] = 1.0 - xi[t,1];
    }
  }
  
}
model {
  // priors
  p ~ beta(2, 2);
  [c1, c2] ~ cauchy(0, 1);
  for(k in 1:2)
    phi[k] ~ normal(0,1);
  sigma ~ cauchy(0, 1);
  xi1_init ~ beta(2, 2);

  target += sum(log(f));
}
