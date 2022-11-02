data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma;
  real alpha;
}
model {
  sigma ~ inv_gamma(1,1);
  y ~ normal(alpha+X*beta, sigma);
}

generated quantities {
  vector[N] ystar;
  vector[N] log_lik;
  for (i in 1:N) {
      ystar[i] = normal_rng(alpha+X[i,]*beta,sigma);
      log_lik[i] = normal_lpdf(y[i]|alpha+X[i,]*beta,sigma);
  }
}
