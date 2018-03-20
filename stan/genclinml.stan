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
  //Laplace prior
    vector g_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
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
  int<lower=0> J1; //cohort
  int<lower=0> J2; //intClust
  vector[Nobs] yobs;
  vector[Ncen] ycen;
  int J1obs[Nobs]; //cohort
  int J1cen[Ncen]; 
  int J2obs[Nobs]; //intClust
  int J2cen[Ncen]; 
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
  
  real zeta1;
  vector[J1] b1;
  real<lower=0> kappa1;
  
  real zeta2;
  vector[J2] b2;
  real<lower=0> kappa2;
  
  vector<lower=0>[J2]tau1_global;
  vector<lower=0>[J2]tau2_global;
    
  vector<lower=0>[J2] tau1_local[M];
  vector<lower=0>[J2] tau2_local[M];
  
  vector[J2] beta_g_raw[M];
  
  real<lower=0> tau_s_b_raw;
  vector<lower=0>[M] tau_b_raw;
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;
  matrix[J2,M] beta_g;
  vector[Nobs] Z_g_obs;
  vector[Ncen] Z_g_cen;
  

  beta_b = g_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;

  for (j in 1:J2){
     beta_g[j,] = to_row_vector(prior_hs_lp(tau1_global[j], tau2_global[j],tau1_local[j],tau2_local[j], nu_local, scale_global, nu_global) .* beta_g_raw[j]); 
  }
  
  alpha = exp(tau_al * alpha_raw);
  
  for(n in 1:Nobs){
    Z_g_obs[n] = Zobs[n, ] * to_vector(beta_g[J2obs[n],]);
  }
  for(n in 1:Ncen){
    Z_g_cen[n] = Zcen[n, ] * to_vector(beta_g[J2cen[n],]); 
  }

}

model {
     yobs ~ weibull(alpha, exp(-( mu + zeta1 + zeta2 + Z_g_obs + b1[J1obs]+ Zobs * beta_b + b2[J2obs])/alpha)); 

    target += weibull_lccdf(ycen | alpha, exp(-( mu + zeta1 + zeta2 + Z_g_cen+ Zcen * beta_b + b1[J1cen]+ b2[J2cen])/alpha)); 
     
  beta_b_raw ~ normal(0.0, 1.0);
  
  for(j in 1:J2){
    beta_g_raw[j] ~ normal(0, 1); 
  }

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
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu + zeta1 + zeta2 + Z_g_obs[n] + Zobs[n,] * beta_b + b1[J1obs[n]]+ b2[J2obs[n]])/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu + zeta1 + zeta2 +Z_g_cen[n] + Zcen[n,] * beta_b + b1[J1cen[n]]+ b2[J2cen[n]])/alpha));
  }
  
}
