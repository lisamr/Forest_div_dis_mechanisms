data{
  int<lower=1> N; //number of observations
  int<lower=1> S; //number of species
  int<lower=1> J; //number of plots
  int<lower=1> Kp; //number of predictor variables, plot
  matrix[J,Kp] Xp; //the model matrix, plot
  int<lower=0, upper=1> I[N];//infected or not
  int<lower=1> SpID[N]; //indices
  int<lower=1> PlotID[N];
  matrix[N, 2] XS; //model matrix, varies by species (int, diversity)
  vector[N] BA;
  matrix[J,J] D; //pairwise distance matrix
}

parameters{
  real bBA;
  vector[Kp] betap;

  //parameters for intercept and slope that varies by species 
  matrix[2, S] zS; // matrix of intercepts and slope
  vector<lower=0>[2] sigmaS; // sd for intercept and slope
  vector[2] betaS; // intercept and slope hyper-priors
  cholesky_factor_corr[2] LS; // Cholesky correlation matrix
  
  //GP pars
  vector[J] zplot;
  real<lower=0> eta;
  real<lower=0> rhoGP;
  real<lower=0> sigma;
}

transformed parameters{ 
  matrix[2, S] z; // non-centered species-varying effects
  vector[N] p; //main model
  
  //Gaussian process intercept (varies by plot)
  cov_matrix[J] Sigma;
  vector[J] w;
  real<lower=0> eta_sq = pow(eta, 2);
  real<lower=0> rho_sq = pow(rhoGP, 2);
  real<lower=0> sig_sq = pow(sigma, 2);
    for (i in 1:(J-1)) {
    for (j in (i + 1):J) {
      Sigma[i, j] = eta_sq * exp(-.5/rho_sq * pow(D[i, j], 2));
      Sigma[j, i] = Sigma[i, j];
    }
  }
  for (k in 1:J) Sigma[k, k] = eta_sq + sig_sq;
  w =  Xp*betap + cholesky_decompose(Sigma) * zplot;
  
  //linear models
   z = diag_pre_multiply(sigmaS, LS) * zS; //varies by species
  for(i in 1:N){
    p[i] = w[PlotID[i]] + (z[1,SpID[i]] + betaS[1]) + (z[2,SpID[i]] + betaS[2])*XS[i, 2] + bBA*BA[i]; //main model
  }
}
model{
  //priors
  bBA ~ normal(0, 1);
  betap ~ normal(0, 1);
  //multinormal prior for varying int & slope
  to_vector(zS) ~ normal(0, 1);
  sigmaS ~ exponential(1);
  betaS ~ normal(0, 1); 
  LS ~ lkj_corr_cholesky(2);
  //priors for GP
  eta ~ normal(0, 1.5);
  sigma ~ normal(0, 1);
  rhoGP ~ normal(0, 3);
  zplot ~ normal(0, 1);

  //likelihood
  I ~ bernoulli_logit(p);
}
generated quantities{
  int<lower=0, upper=1> y_rep[N];
  vector[N] log_lik;
  
  for(i in 1:N){
    y_rep[i] = bernoulli_logit_rng(p[i]);
    log_lik[i] = bernoulli_logit_lpmf(I[i]| p[i]);
  }
}
