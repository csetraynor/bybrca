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
  vector prior_hs_lp(real r1_global, real r2_global, vector r1_local, vector r2_local,  real nu_local, real scale_global, real nu_global){
    vector[num_elements(r1_local)] lambda;
    real tau;
    
    //half-t prior for lambdas 
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
    
    //half-t prior for tau
    r1_global ~ normal(0.0, scale_global);
    r2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);

    lambda = r1_local .* sqrt_vec(r2_local);
    tau = r1_global * sqrt(r2_global);

    return  (lambda * tau);
  }
}

data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int<lower=0> M;
  matrix[Nobs, M] Zobs;
  matrix[Ncen, M] Zcen;
}

transformed data {
  real<lower=0> tau_al;
  real<lower=0> tau_mu;
  int<lower=0> N;
  real<lower=0> scale_global;
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  tau_al = 10.0;
  tau_mu = 10.0;
  N = Nobs + Ncen;
  nu_local = 1; //horshoe prior
  nu_global = 1; //half-cauchy
  scale_global = .007;
  
}

parameters {
  real alpha_raw;
  real mu;
  
  
  real<lower=0> tau1_global;
  real<lower=0> tau2_global;
  vector<lower=0>[M] tau1_local;
  vector<lower=0>[M] tau2_local;
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;
  
  beta_b = prior_hs_lp(tau1_global, tau2_global, tau1_local, tau2_local, nu_local, scale_global, nu_global) .* beta_b_raw; 
  
  alpha = exp(tau_al * alpha_raw);
}

model {
     yobs ~ weibull(alpha, exp(-( mu +Zobs * beta_b )/alpha)); 
     
     target += weibull_lccdf(ycen | alpha, exp(-( mu + Zcen * beta_b)/alpha)); 
  
  beta_b_raw ~ normal(0.0, 1.0);
  
  alpha_raw ~ normal(0.0, 1.0);
  mu ~ normal(0 , tau_mu);
}

generated quantities {
  vector[N] log_lik;
  
  for (n in 1:Nobs){
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu +  Zobs[n,] * beta_b)/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu  +Zcen[n,] * beta_b)/alpha));
  }
  
}


