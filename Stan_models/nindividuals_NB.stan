data {
  int<lower=0> N;
  int<lower=1> K; //number of covariates
  int<lower=0> y[N]; //density
  matrix[N,K] X; //covariate matrix
  int<lower=0> Nsim;
  matrix[Nsim, K] Xsim; //data to predict
}
parameters {
  real a0;
  vector[K] B;
  real<lower=0> phi;
}
transformed parameters{
  vector[N] mu;
  mu = a0 + X*B;
}
model {
  a0 ~ normal(3, 1);
  B ~ normal(0, 1);
  phi ~ normal(0, 1);
  target += neg_binomial_2_log_lpmf(y |mu, phi);
}
generated quantities{
  vector[N] y_rep; //for goodness of fit (same data)
  vector[Nsim] mu_sim; //predicting to new data
  for(i in 1:N){
    y_rep[i] = neg_binomial_2_log_rng(mu[i], phi);    
  }
  mu_sim = exp(a0 + Xsim*B);
}
