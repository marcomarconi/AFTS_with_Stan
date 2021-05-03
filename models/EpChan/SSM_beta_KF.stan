data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[2] m0;
  matrix[2,2] P0;
}
transformed data {
  matrix[2,2] F;
  row_vector[2] H[N]; 
  matrix[2,2] I;
  F = [[1,0],
        [0,1]];
  for(t in 1:N)      
    H[t] = [1, x[t]];      
}
parameters {
  cov_matrix[2] Q;
  real<lower=0> R;
}
transformed parameters {
  vector[2] m[N];
  real S[N];
  real mu[N];
  vector[2] m_pred[N]; 
  matrix[2,2] P[N]; 
  matrix[2,2] P_pred[N];
  real v[N];
  vector[2] K[N]; 
  {
    m_pred[1] = m0;
    P_pred[1] = P0;
    for (t in 1:N) {
      if (t>1) {
        m_pred[t] = F * m[t-1];
        P_pred[t] = F * P[t-1] * F' + Q;
      }
      v[t] = y[t] - H[t] * m_pred[t];
      S[t] = (H[t] * P_pred[t] * H[t]' + R);
      K[t] = (P_pred[t] * H[t]') / S[t];
      m[t] = m_pred[t] + K[t] * v[t];   
      P[t] = P_pred[t] - P_pred[t] * H[t]' * K[t]';
      mu[t] = H[t] * m_pred[t];
    }
  }
}
model {
  Q ~ inv_wishart(2, diag_matrix(rep_vector(1,2)));
  R ~ student_t(4, 0, 1);
  y ~ normal(mu, sqrt(S));
}
generated quantities {

}
