/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter

author Carlos Traynor
*/

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> J1; //cohort
  int<lower=0> J2; //intClust
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int J1obs[Nobs]; //cohort
  int J1cen[Ncen]; 
  int J2obs[Nobs]; //intClust
  int J2cen[Ncen]; 
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  tau_al = 10.0;
  tau_mu = 10.0;
  N = Nobs + Ncen;
}

parameters {
  real alpha_raw;
  real mu;
  
  real zeta1;
  vector[J1] b1;
  real<lower=0> kappa1;
  
  real zeta2;
  vector[J2] b2;
  real<lower=0> kappa2;
}

transformed parameters {
  real<lower=0> alpha;
  alpha = exp(tau_al * alpha_raw);
  
}

model {
  yobs ~ weibull(alpha, exp(-( mu + zeta1 + zeta2 + b1[J1obs]+ b2[J2obs])/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-( mu + zeta1 + zeta2  + b1[J1cen]+ b2[J2cen])/alpha));
  
  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta1 ~ normal(0, 1);
  b1 ~ normal(0, kappa1);
  kappa1 ~ gamma(2, .1);
  
  zeta2 ~ normal(0, 1);
  b2 ~ normal(0, kappa2);
  kappa2 ~ gamma(2, .1);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu + zeta1 + zeta2 + b1[J1obs[n]]+ b2[J2obs[n]])/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu + zeta1 + zeta2 + b1[J1cen[n]]+ b2[J2cen[n]])/alpha));
  }
  
}


