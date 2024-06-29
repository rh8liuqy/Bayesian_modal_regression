functions {
  vector vec_student_t_pdf(vector y, real delta, real sigma, int N) {
    return exp(lgamma((delta+1.0)/2.0) -
           lgamma(delta/2.0) -
           0.5*log(delta*pi()) -
           log(sigma) -
           (delta+1.0)/2.0 * log(1.0 + 1.0/delta * square(y/sigma)));
  }

  real skewed_student_t_lpdf(vector y, real sigma1, real sigma2,
                             real delta1, real delta2, int N) {
    real p1;
    vector[N] v1;
    vector[N] v2;
    vector[N] prob;
    real w;
    w = sigma1*exp(student_t_lpdf(0.0|delta2,0.0,1.0))/
    (sigma1*exp(student_t_lpdf(0.0|delta2,0.0,1.0))+
    sigma2*exp(student_t_lpdf(0.0|delta1,0.0,1.0)));
    v1 = 1.0 - ((y)./abs(y) + 1.0)/2.0;
    v2  = 1.0 - v1;
    prob = 2.0*w*vec_student_t_pdf(y,delta1,sigma1,N).*v1+
    2.0*(1.0-w)*vec_student_t_pdf(y,delta2,sigma2,N).*v2;
    return sum(log(prob));
  }
}

data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
}

parameters{
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> delta1;
  real<lower=0> delta2;
  vector[P] beta;
}

model {
  sigma1 ~ inv_gamma(5,5);
  sigma2 ~ inv_gamma(5,5);
  delta1 ~ inv_gamma(5,5);
  delta2 ~ inv_gamma(5,5);
  target += skewed_student_t_lpdf(y-X*beta|sigma1,sigma2,delta1,delta2,N);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  int z;
  real x1;
  real x2;
  for (i in 1:N) {
      z = bernoulli_rng(sigma1*exp(student_t_lpdf(0.0|delta2,0.0,1.0))/
      (sigma1*exp(student_t_lpdf(0.0|delta2,0.0,1.0))+
      sigma2*exp(student_t_lpdf(0.0|delta1,0.0,1.0))));
      x1 = -abs(student_t_rng(delta1,0,sigma1));
      x2 = abs(student_t_rng(delta2,0,sigma2));
      ystar[i] = x1*z + x2*(1-z);
      ystar[i] = ystar[i] + X[i,]*beta;
      log_lik[i] = skewed_student_t_lpdf(rep_vector(y[i] - X[i,]*beta,1)|
      sigma1,
      sigma2,
      delta1,
      delta2,
      1);
  }
}
