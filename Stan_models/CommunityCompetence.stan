data {
  int<lower=0> N;
  int<lower=1> K; //number of covariates
  vector[N] y; //logged community competence
  matrix[N,K] X; //covariate matrix
  int<lower=0> Nsim;
  matrix[Nsim, K] Xsim; //data to predict
}
parameters {
  real a0;
  vector[K] B;
  real<lower=0> sigma;
}
transformed parameters{
  vector[N] mu;
  mu = a0 + X*B;
}
model {
  a0 ~ normal(6, 2);
  B ~ normal(0, 1);
  sigma ~ exponential(1);
  target += normal_lpdf(y| mu, sigma);
}
generated quantities{
  vector[N] y_rep; //for goodness of fit (same data)
  vector[Nsim] mu_sim; //predicting to new data, not-log scale
  for(i in 1:N){
    y_rep[i] = normal_rng(mu[i], sigma);    
  }
  mu_sim = exp(a0 + Xsim*B);
}
