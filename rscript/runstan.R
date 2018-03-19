##Run Stan
# !diagnostics off
library(dplyr)
library(readr)
library(survival)
library(rstan)
library(loo)
library(caret)
library(Biobase)
memory.limit(1e10)
md <- read_rds("C:/RFactory/bymetabric_files/rdsmetabric/Med_Data_Clean.rds")
gd <- read_rds("C:/RFactory/bymetabric_files/rdsmetabric/Gen_Data.rds")
load("Gen_data_fun.Rdata")
# Run null model
stan_file_null <- "bybrca/stan/null.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stannull <- rstan::stan(stan_file_null,
                        data = gen_stan_data0(md),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits0())
log_liknull <- loo::extract_log_lik(stannull, parameter_name = "log_lik")
loonull <- loo::loo(log_liknull)
print(loonull)
saveRDS(stannull, file = "C:/RFactory/bymetabric_files/bysfit/stannul.rds")
rm(list = c('stannull', 'log_liknull'))


#Run Null with Multilevel for cohort
stan_file_null_ml <- "bybrca/stan/nullml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stannullml <- rstan::stan(stan_file_null_ml,
                        data = gen_stan_data1(md),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits1(J1 = 5, J2 = 11))
log_liknullml <- loo::extract_log_lik(stannullml, parameter_name = "log_lik")
loonullml <- loo::loo(log_liknullml)
print(loonullml)
compare(loonull, loonullml) #preference for the second model!
saveRDS(stannullml, file = "bysfit/nullml.rds")
rm(list = c('stannullml', 'log_liknullml'))

#---------------------------------------
##Run Multilevel with classical clinical vars
stan_file_clin <- "bybrca/stan/clinml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanclinml <- rstan::stan(stan_file_clin,
                        data = gen_stan_data2(md,
                                              formula = "~ size + grade +tumor_stage "),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 1000,
                        init = gen_inits2(J1 = 5, J2 = 11, M = 7))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

log_likclin <- loo::extract_log_lik(stanclinml, parameter_name = "log_lik")
looclin <- loo(log_likclin)
print(looclin)
compare(loonullml, looclin)
saveRDS(stanclinml, file = "bysfit/clin.rds")
rm(list = c('stanclinml', 'log_likclin'))

#---------------------------------------
##Run Ml Stan with genomic
stan_file_gen <- "bybrca/stan/genml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 1
stangene <- rstan::stan(stan_file_gen,
                        data = gen_stan_data3(md,
                                              es = gd),
                        cores = min(nChain, parallel::detectCores()),
                        chains = nChain,
                        iter = 5,
                        init = gen_inits3(J1 = 5, J2 = 11, M = 24368))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

likgene <- loo::extract_log_lik(stangene, parameter_name = "log_lik")
loogene <- loo(likgene)
saveRDS(stangene, file = "bysfit/gene.rds")
rm(list = c('stangene', 'likgene'))


#---------------------------------------
##Run ClinicoGenomic Stan
stan_file_clingen <- "bys/bys/clingenml.stan"
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
nChain <- 4
stanclingene <- rstan::stan(stan_file_clingen,
                            data = gen_stan_data4(md,
                                                  es = brcaES, formula =  "~ stage + er + pr+ her2 "),
                            cores = min(nChain, parallel::detectCores()),
                            chains = nChain,
                            iter = 1000,
                            init = gen_inits4(J = 6, M = 12, M_g = 14666))
# if (interactive())
#   shinystan::launch_shinystan(stanfit)

likclingene <- loo::extract_log_lik(stanclingene, parameter_name = "log_lik")
looclingene <- loo(likclingene)
saveRDS(stanclingene, file = "bysfit/clingene.rds")
rm(list = c('stanclingene', 'likclingene'))
