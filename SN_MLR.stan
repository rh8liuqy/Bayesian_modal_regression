functions {
  real vec_SN_DP_lpdf(vector y, real xi, real omega, real alpha) {
    int N = num_elements(y);
    vector[N] z = (y - xi) / omega;
    real p1 = N*log(2.0);
    real p2 = -N*log(omega);
    real p3 = std_normal_lpdf(z);
    real p4 = std_normal_lcdf(alpha * z);
    
    return p1 + p2 + p3 + p4;
}
  
  real vec_SN_CP_lpdf(vector y, real mu, real sigma, real gamma1) {
    real output;
    real constant_c;
    real mu_z;
    real omega;
    real xi;
    real delta;
    real alpha;
    constant_c = pow(2.0*gamma1/(4.0 - pi()), 1.0/3.0);
    mu_z = constant_c/sqrt(1.0+pow(constant_c,2.0));
    omega = sigma/sqrt(1-pow(mu_z,2.0));
    xi = mu - omega * mu_z;
    delta = mu_z / sqrt(2.0/pi());
    if (delta > 0.0) {
      alpha = 1.0 / sqrt(1.0 / pow(delta,2.0) - 1.0);
    } else {
      alpha = -1.0 / sqrt(1.0 / pow(delta,2.0) - 1.0);
    }
    output = vec_SN_DP_lpdf(y | xi, omega, alpha);
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
  vector[P] beta;
  real alpha;
  real<lower = 0> sigma;
  real<lower = -1, upper = 1> gamma1;
}

model {
  sigma ~ inv_gamma(1,1);
  target += vec_SN_CP_lpdf(y - alpha - X*beta | 0.0, sigma, gamma1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] ystar;
  real constant_c;
  real mu_z;
  real omega;
  real xi;
  real delta;
  real alpha_SN;
  real U1;
  real U2;
  real Z;
  for (i in 1:N) {
    constant_c = pow(2.0*gamma1/(4.0 - pi()), 1.0/3.0);
    mu_z = constant_c/sqrt(1.0+pow(constant_c,2.0));
    omega = sigma/sqrt(1-pow(mu_z,2.0));
    xi = alpha + X[i,]*beta - omega * mu_z;
    delta = mu_z / sqrt(2.0/pi());
    if (delta > 0.0) {
      alpha_SN = 1.0 / sqrt(1.0 / pow(delta,2.0) - 1.0);
    } else {
      alpha_SN = -1.0 / sqrt(1.0 / pow(delta,2.0) - 1.0);
    }
    U1 = abs(std_normal_rng());
    U2 = std_normal_rng();
    Z = (alpha_SN * U1 + U2)/sqrt(1 + pow(alpha_SN, 2.0));
    ystar[i] = xi + omega * Z;
    log_lik[i] = vec_SN_CP_lpdf(rep_vector(y[i] - alpha - X[i,]*beta, 1) | 
    0.0, 
    sigma, 
    gamma1);
  }
}
