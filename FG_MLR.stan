functions {
  real FG_lpdf(vector y, real w1, real scale1, real scale2, int N) {
    vector[N] z1;
    vector[N] z2;
    real p1;
    vector[N] p2;
    real p3;
    vector[N] p4;
    real output;
    real w2;
    w2 = 1.0 - w1;
    z1 = y/scale1;
    z2 = y/scale2;
    p1 = w1/scale1;
    p2 = exp(z1-exp(z1));
    p3 = w2/scale2;
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
  real<lower=0, upper=1> w1;
  real<lower=0> scale1;
  real<lower=0> scale2;
  vector[P] beta;
  real alpha;
}

model {
  scale1 ~ inv_gamma(1,1);
  scale2 ~ inv_gamma(1,1);
  target += FG_lpdf(y-alpha-X*beta | w1,
  scale1,
  scale2,
  N);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  int z;
  real x1;
  real x2;
  for (i in 1:N) {
      z = bernoulli_rng(w1);
      x1 = -gumbel_rng(0,scale1);
      x2 = gumbel_rng(0,scale2);
      ystar[i] = x1*z + x2*(1-z);
      ystar[i] = ystar[i] + alpha + X[i,]*beta;
      log_lik[i] = FG_lpdf(rep_vector(y[i]-alpha-X[i,]*beta,1) | w1,
      scale1,
      scale2,
      1);
  }
}
