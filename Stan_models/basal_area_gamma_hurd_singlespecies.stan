data {
  int<lower=0> N;
  int<lower=1> K; //number of covariates
  real<lower=0> y[N]; //basal area density
  matrix[N,K] X; //covariate matrix
  int<lower=0> Nsim;
  matrix[Nsim,K] Xsim; //data to predict
}
parameters {
  real theta;
  real a0;
  vector[K] B;
  vector[K] Bz;
  real<lower=0> va;
}
transformed parameters{
  vector[N] mu;
  mu = exp(a0 + X*B);
}
model {
  theta ~ normal(0, 1);
  a0 ~ normal(-1, 1);
  B ~ normal(0, 1);
  Bz ~ normal(0, 1);
  va ~ normal(0, 1);
  
  for (i in 1:N){
    if(y[i] == 0){
      target += log(inv_logit(theta + X[i]*Bz)); 
    } else{
      target += log1m(inv_logit(theta + X[i]*Bz)); 
      target += gamma_lpdf(y[i] | mu[i]*mu[i]/va, mu[i]/va);
    }
  }
}
generated quantities{
  vector[N] y_rep;
  vector[Nsim] mu_sim;
  vector[Nsim] p_sim;
  for(i in 1:N){
    if(bernoulli_logit_rng(theta + X[i]*Bz)){
      y_rep[i] = 0;
    }else{
      y_rep[i] = gamma_rng(mu[i]*mu[i]/va, mu[i]/va);
    }
  }
  mu_sim = exp(a0 + Xsim*B);
  p_sim = 1 - inv_logit(theta + Xsim*Bz); //probabilyt of occurence
}
