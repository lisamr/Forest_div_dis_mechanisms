data {
  int<lower=0> N;
  int<lower=1> K; //number of covariates
  real<lower=0> y[N]; //basal area density
  matrix[N,K] X; //covariate matrix
  int<lower=0> Nsim;
  matrix[Nsim, K] Xsim; //data to predict
}
parameters {
  real a0;
  vector[K] B;
  real<lower=0> va;
}
transformed parameters{
  vector[N] mu;
  mu = exp(a0 + X*B);
}
model {
  a0 ~ normal(-1, 1);
  B ~ normal(0, 1);
  va ~ normal(0, 1);
  
  for (i in 1:N){
    target += gamma_lpdf(y[i] | mu[i]*mu[i]/va, mu[i]/va);
  }
  
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep; //for goodness of fit (same data)
  vector[Nsim] mu_sim; //predicting to new data
  vector[Nsim] y_sim; //predicting to new data
  for(i in 1:N){
      y_rep[i] = gamma_rng(mu[i]*mu[i]/va, mu[i]/va);
      log_lik[i] = gamma_lpdf(y[i] | mu[i]*mu[i]/va, mu[i]/va);
  }
  mu_sim = exp(a0 + Xsim*B);
  for(i in 1:Nsim) 
    y_sim[i] = gamma_rng(mu_sim[i]*mu_sim[i]/va, mu_sim[i]/va);
}
