/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter

author Carlos Traynor, heavily inspired and grateful to Tomi Peltola and Jacqueline Buros
*/
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }

    return res;
  }

  //Horseshoe Prior
  vector prior_hs_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, vector r1_localplus, vector r2_localplus, real nu){
    vector[num_elements(r1_local)] lambda;
    vector[num_elements(r1_local)] lambdaplus;
    real tau;
    
    //half-t prior for lambdas (nu = 1 corresponds to horseshoe+)
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);
    r1_localplus ~ normal(0.0, 1.0);
    r2_localplus ~ inv_gamma(0.5 * nu, 0.5 * nu);
    //half-cauchy prior for tau
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    lambda = r1_local .* sqrt_vec(r2_local);
    lambdaplus = r1_localplus .* sqrt_vec(r2_localplus);
    tau = r1_global * sqrt(r2_global);

    return  (lambda .* lambdaplus * tau);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  int<lower=0> J1; //cohort
  int<lower=0> J2; //intClust
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int Jobs1[Nobs]; //cohort
  int Jcen1[Ncen]; 
  int Jobs2[Nobs]; //intClust
  int Jcen2[Ncen]; 
  int<lower=0> M;
  matrix[Nobs, M] Zobs;
  matrix[Ncen, M] Zcen;
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
  
  real<lower=0> tau1_global;
  vector<lower=0>[M] tau2_global;
  
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;
  real sigma; //noise std
  
  sigma= exp(logsigma);
  beta_b = beta_b_raw  .* prior_hs_lp(tau1_global, tau2_global, tau1_local, tau2_local, nu_local, scale_global, nu_global) ;
  alpha = exp(tau_al * alpha_raw);

}

model {
  yobs ~ weibull(alpha, exp(-( mu + zeta1 + zeta2 + Zobs * beta_b + b1[Jobs1]+ b2[Jobs2])/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-( mu + zeta1 + zeta2 + Zcen * beta_b + b1[Jcen1]+ b2[Jcen2])/alpha));
  
  beta_b_raw ~ normal(0.0, 1.0);
  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
  
  zeta1 ~ normal(0, 1);
  b1 ~ normal(0, kappa1);
  kappa1 ~ gamma_lpdf(2, .1);
  
  zeta2 ~ normal(0, 1);
  b2 ~ normal(0, kappa2);
  kappa2 ~ gamma_lpdf(2, .1);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu + zeta1 + zeta2 + Zobs[n,] * beta_b + b1[Jobs1[n]]+ b2[Jobs2[n]])/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu + zeta1 + zeta2 +Zcen[n,] * beta_b + b1[Jcen1[n]]+ b2[Jcen2[n]])/alpha));
  }
  
}


