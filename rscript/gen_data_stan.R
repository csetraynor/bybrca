# !diagnostics off
rm(list = ls())
library(readr)
library(dplyr)
library(caret)
library(Biobase)

##--- Gen Data Null
gen_stan_data0 <- function(data) {

  ind = (data$status == 1)
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind]
  )
}
gen_inits0 <- function() {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1)
    )
}
 # into_data <- gen_stan_data(md)
 # glimpse(into_data)
 # rstan::stan_rdump(ls(into_data), file = "checking.data.R",
 #                   envir = list2env(into_data))

##--- Gen Stan Data Null with Intercept ---#
gen_stan_data1 <- function(data) {
  # Subgroup indicator
  J1group = as.integer(data$cohort)
  
  J2group = as.integer(data$intclust)
  ind = data$status == 1
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    J1 = n_distinct(data$cohort),
    J1obs = as.numeric(J1group[ind]),
    J1cen = as.numeric(J1group[!ind]),
    J2 = n_distinct(data$intclust),
    J2obs = as.numeric(J2group[ind]),
    J2cen = as.numeric(J2group[!ind])
  )
}
 gen_stan_data1(md) %>% glimpse
 #glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits1 <- function(J1, J2) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      zeta1 = rnorm(1),
      b1 = rnorm(J1),
      kappa1 = abs(rcauchy(1, 0, 2)),
      
      zeta2 = rnorm(1),
      b2 = rnorm(J2),
      kappa2 = abs(rcauchy(1, 0, 2))
    )
}
# inits <- gen_inits(J = 6)
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))

##--- Gen Stan Data Clinical Model ---#
gen_stan_data2 <- function(data, formula = as.formula(~1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  #Covariates (clinical) Matrix
  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- as.matrix(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }
  tmp <- tbl_df(Z)
  if ("size" %in% colnames(tmp)){
    tmp <- tmp %>% select(-size)
  }
  prop <- apply(tmp, 2, sum) / apply(tmp, 2, length) 
  tmp <- sweep(tmp, 2, prop) #centering
  Z <- tbl_df(Z) %>% select(size) %>% cbind(tmp) %>% as.matrix
  
  # Subgroup indicator
  J1group = as.integer(data$cohort)
  J2group = as.integer(as.factor(data$intclust))
  
  #Censoring indicator
  ind = data$status == 1
  Zobs <- Z[ind,]
  Zcen <- Z[!ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[!ind,]), M)),
    J1 = n_distinct(data$cohort),
    J1obs = as.numeric(J1group[ind]),
    J1cen = as.numeric(J1group[!ind]),
    J2 = n_distinct(data$intclust),
    J2obs = as.numeric(J2group[ind]),
    J2cen = as.numeric(J2group[!ind])
  )
}
# gen_stan_data2(md, formula = "~ tumor_stage + size + grade") %>% 
#   glimpse

# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits2 <- function(J1, J2, M) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      beta_b_raw = rnorm(M),
      
      zeta1 = rnorm(1),
      b1 = rnorm(J1),
      kappa1 = abs(rcauchy(1, 0, 2)),
      
      zeta2 = rnorm(1),
      b2 = rnorm(J2),
      kappa2 = abs(rcauchy(1, 0, 2))
    )
}
# init
# rstan::stan_rdump(ls(inits), file = "checking.init.R",
#                   envir = list2env(inits))

##--- Gen Stan Data Genomic Model---#
gen_stan_data3 <- function(data, eset = NA) {
  
  Z <- t(exprs(eset))
  M <- ncol(Z)
  
  # Subgroup indicator
  J1group = as.integer(data$cohort)
  J2group = as.integer(as.factor(data$intclust))
  
  #Censoring indicator
  ind = data$status == 1
  Zobs <- Z[ind,]
  Zcen <- Z[!ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[!ind,]), M)),
    J1 = n_distinct(data$cohort),
    J1obs = as.numeric(J1group[ind]),
    J1cen = as.numeric(J1group[!ind]),
    J2 = n_distinct(data$intclust),
    J2obs = as.numeric(J2group[ind]),
    J2cen = as.numeric(J2group[!ind])
  )
}
# Z <- md %>%
#   model.matrix(~ stage + er + pr+ her2 + menopause , data = .)
# attr(Z, "dimnames")[[2]][-1]
gen_stan_data3(md, eset = gd) %>% glimpse
 into_data = gen_stan_data3(md, eset = gd)
  rstan::stan_rdump(ls(into_data), file = "gene.data.R",
                    envir = list2env(into_data))
#Inits for select
gen_inits3 <- function(J1, J2, M) {
 # function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      beta_g_raw = array(rnorm(M*J2), dim = c(J2, M)),
      
      tau1_global = 0.1*abs(rnorm(J2)),
      tau2_global = 0.1*abs(rnorm(J2)),
      tau1_local = abs(array(rnorm(M*J2), dim = c(J2, M))),
      tau2_local = abs(array(rnorm(M*J2), dim = c(J2, M))),
      
      zeta1 = rnorm(1),
      b1 = rnorm(J1),
      kappa1 = abs(rcauchy(1, 0, 2)),
      
      zeta2 = rnorm(1),
      b2 = rnorm(J2),
      kappa2 = abs(rcauchy(1, 0, 2))
    )
}
 init =  gen_inits3(J1 = 5, J2 = 11, M = 4715)
 rstan::stan_rdump(ls(init), file = "gen.init.R",
                   envir = list2env(init))
##--- Gen Stan Data Clinico Genomic Model---#
gen_stan_data4 <- function(data, eset = NA , formula = as.formula(~1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  # Subgroup indicator
  J1group = as.integer(data$cohort)
  J2group = as.integer(as.factor(data$intclust))
  
  #Covariates (clinical) Matrix
  Z <- data %>%
    model.matrix(formula, data = .)
  M <- ncol(Z)
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Z))
      Z <- as.matrix(Z[,-1], dim = c(nrow(data), M - 1))
    M <- ncol(Z)
  }
  tmp <- tbl_df(Z)
  if ("size" %in% colnames(tmp)){
    tmp <- tmp %>% select(-size)
  }
  #centering
  prop <- apply(tmp, 2, sum) / apply(tmp, 2, length) 
  tmp <- sweep(tmp, 2, prop)
  Z <- tbl_df(Z) %>% select(size) %>% cbind(tmp) %>% as.matrix
  
  Z_g <- t(exprs(eset))
  M_g <- ncol(Z_g)
  
  #Censoring indicator
  ind = data$status == 1
  Zobs <- Z[ind,]
  Zcen <- Z[!ind,]
  
  Zobs_g <- Z_g[ind,]
  Zcen_g <- Z_g[!ind,]
  
  stan_data <- list(
    Nobs = nrow(data[ind,]),
    Ncen = nrow(data[!ind,]),
    yobs = data$time[ind],
    ycen = data$time[!ind],
    M = M,
    Zobs = array(Zobs, dim = c(nrow(data[ind,]), M)),
    Zcen = array(Zcen, dim = c(nrow(data[!ind,]), M)),
    M_g = M_g,
    Zobs_g = array(Zobs_g, dim = c(nrow(data[ind,]), M_g)),
    Zcen_g = array(Zcen_g, dim = c(nrow(data[!ind,]), M_g)),
    J1 = n_distinct(data$cohort),
    J1obs = as.numeric(J1group[ind]),
    J1cen = as.numeric(J1group[!ind]),
    J2 = n_distinct(data$intclust),
    J2obs = as.numeric(J2group[ind]),
    J2cen = as.numeric(J2group[!ind])
  )
}

# attr(Z, "dimnames")[[2]][-1]
# into_data <- gen_stan_data4(md, eset = brcaES, formula =  "~ stage + er + pr+ her2 ")
# glimpse(into_data)
# rstan::stan_rdump(ls(into_data), file = "checking.data.R",
#                   envir = list2env(into_data))

gen_inits4 <- function(J1, J2, M, M_g) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      
      tau_s_b_raw = 0.1*abs(rnorm(1)),
      tau_b_raw = abs(rnorm(M)),
      
      beta_b_raw = rnorm(M),
      
      tau1_global = 0.1*abs(rnorm(J2)),
      tau2_global = 0.1*abs(rnorm(J2)),
      tau1_local = abs(array(rnorm(M*J2), dim = c(M, J2))),
      tau2_local = abs(array(rnorm(M*J2), dim = c(M, J2))),
      
      beta_g_raw = abs(rnorm(M_g)),
      
      zeta1 = rnorm(1),
      b1 = rnorm(J1),
      kappa1 = abs(rcauchy(1, 0, 2)),
      
      zeta2 = rnorm(1),
      b2 = rnorm(J2),
      kappa2 = abs(rcauchy(1, 0, 2))
    )
}

# md <- read_rds("Med_Data_Clean.rds")
# gd <- read_rds("Gen_Data.rds")
# gen_stan_data2(md, formula = "~ tumor_Stage +") %>% glimpse
save(list = ls(), file = "Gen_data_fun.Rdata")
rm(list = ls())
