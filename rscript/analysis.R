#Analysis
library(loo)
library(rstanarm)
library(bayesplot)
stannull <- read_rds("C:/RFactory/bymetabric_files/bysfit/stannul.rds")
log_liknull <- loo::extract_log_lik(stannull, parameter_name = "log_lik")
loonull <- loo::loo(log_liknull)
rm(list = c('stannull', 'log_liknull'))

stannullml <- read_rds("C:/RFactory/bymetabric_files/bysfit/nullml.rds")
log_liknullml <- loo::extract_log_lik(stannullml, parameter_name = "log_lik")
loonullml <- loo::loo(log_liknullml)
rm(list = c('stannullml', 'log_liknullml'))

stanclin <- read_rds("C:/RFactory/bymetabric_files/bysfit/clin.rds")
log_likclin <- loo::extract_log_lik(stanclin, parameter_name = "log_lik")
looclin <- loo::loo(log_likclin)
rm(list = c('stanclin', 'log_likclin'))

data.frame(Model = c("Null", "Varying Intercept", "Clinical"),
           elpd_loo = c(loonull$elpd_loo, loonullml$elpd_loo,  looclin$elpd_loo),
           se = c(loonull$se_elpd_loo, loonullml$se_elpd_loo, looclin$se_elpd_loo))
compare(loonull, loonullml)

yrep <- posterior_predict(stannull)
psis <- psislw(-log_lik(fit), cores = 2)
