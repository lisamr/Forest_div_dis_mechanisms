data{
  int<lower=1> N; //number of observations
  int<lower=1> K; //number of predictor variables
  matrix[N,K] X; //the model matrix 
  int I[N];//the response variable
  int n[N]; 
  matrix[N,N] D; //pairwise distance matrix
}
parameters{
  real a0; 
  vector[K] beta; //the regression parameters
  real<lower=0> phi; //dispersion parameter
  //GP pars
  vector[N] z;
  real<lower=0> eta;
  real<lower=0> rho;
  //real<lower=0> sigma; //parameter fixed here because no repeat plots in this model

}
transformed parameters{ 
  vector[N] p; //mean of the linear model
  vector[N] Alpha; //shape parameters for the beta distribtuion
  vector[N] Beta;
  real theta;
  cov_matrix[N] Sigma;
  vector[N] w;
  real<lower=0> eta_sq;
  real<lower=0> rho_sq;
  //real<lower=0> sig_sq;
  
  //Gaussian process intercept
  eta_sq = pow(eta, 2);
  rho_sq = pow(rho, 2);
  //sig_sq = pow(sigma, 2);
    for (i in 1:(N-1)) {
    for (j in (i + 1):N) {
      Sigma[i, j] = eta_sq * exp(-.5/rho_sq * pow(D[i, j], 2));
      Sigma[j, i] = Sigma[i, j];
    }
  }
  for (k in 1:N) Sigma[k, k] = eta_sq + .000001; //sig_sq;
  w = cholesky_decompose(Sigma) * z;
  
  
  //linear models
  theta = phi + 2; //want the mass near 2 to get an even spread between 0 and 1
  p = inv_logit(a0 + X*beta + w);
  Alpha = p*theta;
  Beta = (1-p)*theta;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, .5);
  phi ~ exponential(.5); //prior from McElreath ed. 2, pg. 371. I made it wider by a little.
  eta ~ normal(0, 1);
  //sigma ~ normal(0, 1);
  rho ~ normal(0, 1);
  z ~ normal(0, 1);
  
  //likelihood, holding out some data
  for(i in 1:N){
        target += beta_binomial_lpmf(I | n, Alpha, Beta);
  }
}

generated quantities{
  real log_lik[N];
  int<lower=0> y_rep[N];
  for(i in 1:N) y_rep[i] = beta_binomial_rng(n[i], Alpha[i], Beta[i]);

  for(i in 1:N){
      log_lik[i] = beta_binomial_lpmf(I[i] | n[i], Alpha[i], Beta[i]);
  }
}
