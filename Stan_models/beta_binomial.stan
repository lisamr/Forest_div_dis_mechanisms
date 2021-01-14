data{
  int<lower=1> N; //number of observations
  int<lower=1> K; //number of predictor variables
  matrix[N,K] X; //the model matrix 
  int I[N];//the response variable
  int n[N]; 

}
parameters{
  real a0; 
  vector[K] beta; //the regression parameters
  real<lower=0> phi; //dispersion parameter

}
transformed parameters{ 
  vector[N] p; //mean of the linear model
  vector[N] Alpha; //shape parameters for the beta distribtuion
  vector[N] Beta;
  real theta;
  
  //linear models
  theta = phi + 2; //want the mass near 2 to get an even spread between 0 and 1
  p = inv_logit(a0 + X*beta);
  Alpha = p*theta;
  Beta = (1-p)*theta;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, .5);
  phi ~ exponential(.5); //prior from McElreath ed. 2, pg. 371. I made it wider by a little.
  
  target += beta_binomial_lpmf(I | n, Alpha, Beta);
}

generated quantities{
  real log_lik[N];
  int<lower=0> y_rep[N];
  for(i in 1:N) y_rep[i] = beta_binomial_rng(n[i], Alpha[i], Beta[i]);

  for(i in 1:N){
      log_lik[i] = beta_binomial_lpmf(I[i] | n[i], Alpha[i], Beta[i]);
  }
}
