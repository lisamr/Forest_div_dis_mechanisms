data {
  int<lower=0> N;
  int<lower=1> K; //number of covariates + intercept
  int<lower=0> J; 
  int<lower=0> S;
  real<lower=0> y[N]; //basal area density
  matrix[N,K] X; //covariate matrix
  int<lower=0, upper=S> spID[N];
  int<lower=0, upper=J> plotID[N];
  int<lower=0> Nsim;
  matrix[Nsim,K] Xsim; //data to predict
  int<lower=0, upper=S> spID_sim[Nsim];
}
parameters {
  real<lower=0> nu;
  real<lower=0> nuz;
  real<lower=0> va;
  real<lower=0> sd_j;
  vector[J] z_j;
  vector[K] B_bar;
  vector[K] B[S];
  corr_matrix[K] rho;
  vector<lower=0>[K] sigma_p;
  real<lower=0> sd_jz;
  vector[J] z_jz;
  vector[K] Bz_bar;
  vector[K] Bz[S];
  corr_matrix[K] rhoz;
  vector<lower=0>[K] sigma_pz;
}
transformed parameters{
  vector[J] aj;
  vector[J] ajz;
  vector[N] mu;
  vector[N] pz;
  
  aj = sd_j*z_j; //noncentered plot intercept
  ajz = sd_jz*z_jz; //noncentered plot intercept, zero probability
  for(i in 1:N){
    mu[i] = exp(X[i]*B[spID[i]] + aj[plotID[i]]); //mean model for gamma process
    pz[i] = inv_logit(X[i]*Bz[spID[i]] + ajz[plotID[i]]);//mean model for bernoulli process
  }

}
model {
  //gamma process priors
  nu ~ gamma(2, .1); //recommended by Stan manual
  va ~ normal(0, 2);
  sd_j ~ normal(0, 1);
  z_j ~ normal(0, 1);
  B_bar ~ normal(0, 1); 
  rho ~ lkj_corr(2);
  sigma_p ~ normal(0, 5);
  B ~ multi_student_t(nu, B_bar, quad_form_diag(rho, sigma_p));
  //bernoulli process priors
  nuz ~ gamma(2, .1);
  sd_jz ~ normal(0, 1);
  z_jz ~ normal(0, 1);
  Bz_bar ~ normal(0, 1);
  rhoz ~ lkj_corr(2);
  sigma_pz ~ normal(0, 5);
  Bz ~ multi_student_t(nuz, Bz_bar, quad_form_diag(rhoz, sigma_pz));
  
  for (i in 1:N){
    if(y[i] == 0){
      target += log(pz[i]); 
    } else{
      target += log1m(pz[i]); 
      target += gamma_lpdf(y[i] | mu[i]*mu[i]/va, mu[i]/va);
    }
  }
}
generated quantities{
  vector[N] y_rep;
  vector[Nsim] mu_sim;
  vector[Nsim] pz_sim;
  for(i in 1:N){
    if(bernoulli_logit_rng(pz[i])){
      y_rep[i] = 0;
    }else{
      y_rep[i] = gamma_rng(mu[i]*mu[i]/va, mu[i]/va);
    }
  }
    for(i in 1:Nsim){
    mu_sim[i] = exp(Xsim[i]*B[spID_sim[i]]); //predictions at the 'average' plot
    pz_sim[i] = 1 - inv_logit(Xsim[i]*Bz[spID_sim[i]]);//occurence at the 'average' plot
  }
}
