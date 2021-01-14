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
}

parameters{
  real bBA;
  vector[Kp] betap;
  real<lower=0> sdPlot;
  vector[J] zPlot;

  //parameters for intercept and slope that varies by species 
  matrix[2, S] zS; // matrix of intercepts and slope
  vector<lower=0>[2] sigmaS; // sd for intercept and slope
  vector[2] betaS; // intercept and slope hyper-priors
  cholesky_factor_corr[2] LS; // Cholesky correlation matrix
}

transformed parameters{ 
  vector[J] aPlot; //noncenter plot-varying effects
  matrix[2, S] z; // non-centered species-varying effects
  vector[N] p; //main model
  //linear models
   aPlot = Xp*betap + sdPlot*zPlot; //varies by plot
   z = diag_pre_multiply(sigmaS, LS) * zS; //varies by species
  for(i in 1:N){
    p[i] = aPlot[PlotID[i]] + (z[1,SpID[i]] + betaS[1]) + (z[2,SpID[i]] + betaS[2])*XS[i, 2] + bBA*BA[i]; //main model
  }
}
model{
  //priors
  bBA ~ normal(0, 1);
  betap ~ normal(0, 1);
  sdPlot ~ normal(0, 1);
  zPlot ~ normal(0,1);
  //multinormal prior for varying int & slope
  to_vector(zS) ~ normal(0, 1);
  sigmaS ~ exponential(1);
  betaS ~ normal(0, 1); 
  LS ~ lkj_corr_cholesky(2);

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
