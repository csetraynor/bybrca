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

  vector g_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);

    return r_global * sqrt_vec(r_local);
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
  
  real<lower=0> tau_s_b_raw;
  vector<lower=0>[M] tau_b_raw;
  vector[M] beta_b_raw;
}

transformed parameters {
  real<lower=0> alpha;
  vector[M] beta_b;

  
  beta_b = g_prior_lp(tau_s_b_raw, tau_b_raw) .* beta_b_raw;
  alpha = exp(tau_al * alpha_raw);
  

}

model {
  yobs ~ weibull(alpha, exp(-( mu + zeta1 + zeta2 + Zobs * beta_b + b1[J1obs]+ b2[J2obs])/alpha));
  target += weibull_lccdf(ycen | alpha, exp(-( mu + zeta1 + zeta2 + Zcen * beta_b + b1[J1cen]+ b2[J2cen])/alpha));
  
  beta_b_raw ~ normal(0.0, 1.0);
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
    log_lik[n] = weibull_lpdf(yobs[n] | alpha, exp(-( mu + zeta1 + zeta2 + Zobs[n,] * beta_b + b1[J1obs[n]]+ b2[J2obs[n]])/alpha));
  }
  for (n in 1:Ncen){
    log_lik[Nobs + n] = weibull_lccdf(ycen[n]| alpha, exp(-( mu + zeta1 + zeta2 +Zcen[n,] * beta_b + b1[J1cen[n]]+ b2[J2cen[n]])/alpha));
  }
  
}


