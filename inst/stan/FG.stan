functions {
  real FG_lpdf(vector y, real w1, real sigma1, real sigma2, int N) {
    vector[N] z1;
    vector[N] z2;
    real p1;
    vector[N] p2;
    real p3;
    vector[N] p4;
    real output;
    real w2;
    w2 = 1.0 - w1;
    z1 = y/sigma1;
    z2 = y/sigma2;
    p1 = w1/sigma1;
    p2 = exp(z1-exp(z1));
    p3 = w2/sigma2;
    p4 = exp(-z2-exp(-z2));
    output = sum(log(p1*p2 + p3*p4));
    return output;
  }
}

data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
}

parameters {
  real<lower=0, upper=1> w;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector[P] beta;
}

model {
  sigma1 ~ inv_gamma(5,5);
  sigma2 ~ inv_gamma(5,5);
  target += FG_lpdf(y-X*beta | w,
  sigma1,
  sigma2,
  N);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  int z;
  real x1;
  real x2;
  for (i in 1:N) {
      z = bernoulli_rng(w);
      x1 = -gumbel_rng(0,sigma1);
      x2 = gumbel_rng(0,sigma2);
      ystar[i] = x1*z + x2*(1-z);
      ystar[i] = ystar[i] + X[i,]*beta;
      log_lik[i] = FG_lpdf(rep_vector(y[i] - X[i,]*beta,1) | w,
      sigma1,
      sigma2,
      1);
  }
}
