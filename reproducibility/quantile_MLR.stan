data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
}

transformed data {
  real p = 0.5; //median regression
}

parameters {
  vector[P] beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  sigma ~ inv_gamma(1,1);
  target += skew_double_exponential_lpdf(y-alpha-X*beta|0.0,2.0*sigma,p);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  for (i in 1:N) {
      ystar[i] = skew_double_exponential_rng(0.0,2.0*sigma,p);
      ystar[i] = ystar[i] + alpha + X[i,]*beta;
      log_lik[i] = skew_double_exponential_lpdf(y[i]-alpha-X[i,]*beta|
      0.0,
      2.0*sigma,
      p);
  }
}
