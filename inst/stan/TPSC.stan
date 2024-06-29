functions {
  vector vec_student_t_pdf(vector y, real delta,real sigma, int N) {
    return exp(lgamma((delta+1.0)/2.0) -
           lgamma(delta/2.0) -
           0.5*log(delta*pi()) -
           log(sigma) -
           (delta+1.0)/2.0 * log(1.0 + 1.0/delta * square(y/sigma)));
  }

  real skewed_student_t_lpdf(vector y, real w, real delta,
                             real sigma, int N) {
    real p1;
    vector[N] v1;
    vector[N] v2;
    vector[N] prob;
    v1 = 1.0 - (y./abs(y) + 1.0)/2.0;
    v2 = 1.0 - v1;
    prob = (2.0*w*vec_student_t_pdf(y,delta,sigma*sqrt(w/(1.0-w)),N)).*v1 +
    (2.0*(1-w)*vec_student_t_pdf(y,delta,sigma*sqrt((1.0-w)/w),N)).*v2;
    return sum(log(prob));
  }
}

data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
}

parameters {
  real<lower=0.0001,upper=0.9999> w;
  real<lower=0> delta;
  real<lower=0> sigma;
  vector[P] beta;
}

model {
  w ~ uniform(0,1);
  delta ~ inv_gamma(5,5);
  sigma ~ inv_gamma(5,5);
  target += skewed_student_t_lpdf(y-X*beta | w, delta, sigma, N);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  int z;
  real x1;
  real x2;
  for (i in 1:N) {
      z = bernoulli_rng(w);
      x1 = -abs(student_t_rng(delta,0,sigma*sqrt(w/(1.0-w))));
      x2 = abs(student_t_rng(delta,0,sigma*sqrt((1.0-w)/w)));
      ystar[i] = x1*z + x2*(1-z);
      ystar[i] = ystar[i] + X[i,]*beta;
      log_lik[i] = skewed_student_t_lpdf(rep_vector(y[i]-X[i,]*beta,1)|
      w,
      delta,
      sigma,
      1);
  }
}
