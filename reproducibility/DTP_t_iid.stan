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
  vector[N] y;
}

parameters{
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> delta1;
  real<lower=0> delta2;
  real theta;
}

model {
  sigma1 ~ inv_gamma(5,5);
  sigma2 ~ inv_gamma(5,5);
  delta1 ~ inv_gamma(5,5);
  delta2 ~ inv_gamma(5,5);
  target += skewed_student_t_lpdf(y-theta|sigma1,sigma2,delta1,delta2,N);
}
