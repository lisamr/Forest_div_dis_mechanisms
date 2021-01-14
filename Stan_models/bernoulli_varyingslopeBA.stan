data{
  int<lower=1> N; //number of observations
  int<lower=1> S; //number of species
  int<lower=1> J; //number of plots
  int<lower=1> Kp; //number of predictor variables, plot
  matrix[J,Kp] Xp; //the model matrix, plot
  int<lower=0, upper=1> I[N];//infected or not
  int<lower=1> SpID[N]; //indices
  int<lower=1> PlotID[N];
  vector[J] diversity;
  vector[N] BA;
}

parameters{
  vector[Kp] betap;
  real<lower=0> sdPlot;
  vector[J] zPlot;
  vector<lower=0>[2] sigma_pars;
  vector[S] aSp;
  vector[S] bd;
  real bBA;
  real bbar;
  real abar;
  corr_matrix[2] rho;
}

transformed parameters{ 
  //noncenter random effects
  vector[J] aPlot;
  vector[N] p;
  
  //linear models
   aPlot = Xp*betap + sdPlot*zPlot; //plot-level

  for(i in 1:N){
      p[i] = aPlot[PlotID[i]] + aSp[SpID[i]] + bBA*BA[i] + bd[SpID[i]]*diversity[PlotID[i]]; //observation-level
  }
  
}
model{
  //priors
  betap ~ normal(0, .5);
  bBA ~ normal(0, .5);
  zPlot ~ normal(0,1);
  sdPlot ~ normal(0, 1);
  //multinormal prior for varying int & slope
  abar ~ normal(0, 1);
  bbar ~ normal(0, .5);
  sigma_pars ~ exponential(1);
  rho ~ lkj_corr(2);
  {
    vector[2] YY[S];
    vector[2] MU;
    MU = [abar, bbar]';
    for ( s in 1:S ) YY[s] = [aSp[s] , bd[s]]';
    YY ~ multi_normal(MU, quad_form_diag(rho, sigma_pars) );
  }

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
