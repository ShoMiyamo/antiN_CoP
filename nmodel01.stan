data {
  int<lower=0> N;          // number of measurements
  vector[N] y;             // measurements
  vector[N] t;             // time of the measurement
}

parameters {
  real h;                  // max response if lambda = 0
  real<lower=0> alpha;              // shape
  real<lower=0> beta;               // inverse scale
  real<lower=0> lambda;             // decay rate
  real<lower=0> sigma;     // measurement error
}

model {
  // Priors
  h ~ normal(0, 5);
  alpha ~ student_t(4, 0, 0.5);
  beta ~ student_t(4, 0, 0.5);
  lambda ~ student_t(4, 0, 0.5);
  sigma ~ student_t(4, 0, 0.5);
  
  // Likelihood
  for(i in 1:N) {
    y[i] ~ normal(h * gamma_cdf(t[i], alpha, beta) * exp(-lambda * t[i]), sigma);
  }
}

generated quantities {
  vector[N] y_pred;  // predicted values of y

  for(i in 1:N) {
    y_pred[i] = normal_rng(h * gamma_cdf(t[i], alpha, beta) * exp(-lambda * t[i]), sigma);
  }
}
