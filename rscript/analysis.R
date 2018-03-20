#Analysis
library(loo)
library(rstanarm)
library(bayesplot)
library(rethinking)
library(bayesplot)
library(ggplot2)
library(rstanarm)
stannull <- read_rds("C:/RFactory/bymetabric_files/bysfit/stannul.rds")
stannull <- read_rds("bysfit/stannul.rds")
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

parsofinterest = c( "b1[1]", 
                   "b1[2]", "b1[3]", "b1[4]","b1[5]",
                   "b2[1]","b2[2]","b2[4]","b2[5]","b2[6]",
                   "b2[7]","b2[8]","b2[9]","b2[10]","b2[11]","b2[3]")

posterior <- as.array(stannullml, pars = parsofinterest)
dimnames(posterior)$parameters = c("cohort:1","cohort:2","cohort:3","cohort:4","cohort:5","IntClust:1", "IntClust:2","IntClust:3","IntClust:4ER+","IntClust:4ER-","IntClust:5","IntClust:6","IntClust:7","IntClust:8","IntClust:9", "IntClust:10")

theme_set(theme_bw())
#central posterior uncertainty intervals
mcmc_intervals(posterior) +
  ggplot2::labs(
    title = "Posterior distributions",
    subtitle = "with medians, 50% and 90% intervals",
    y = "Shared Frailty parameters"
  )+ vline_0()
#uncertainty intervals as shaded areas
mcmc_areas(
  posterior, 
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

help("MCMC-intervals")
color_scheme_set("red")
mcmc_intervals(posterior)
#Compute elpddiff
elpddiffml = compare(loonullml, loonull )
elpddiffclin = compare(looclin, loonull)
d = data.frame(Model = c("Multilevel null", "Multilevel clinical"),
               elpd_diff = c(elpddiffml[1],elpddiffclin[1]),
               elpd_diff_upper = c(elpddiffml[1] + 1.96*elpddiffml[2], elpddiffclin[1] + 1.96*elpddiffclin[2]),
               elpd_diff_lower = c(elpddiffml[1] - 1.96*elpddiffml[2], elpddiffclin[1] - 1.96*elpddiffclin[2]))
d$Model <- factor(d$Model, levels = c("Multilevel null", "Multilevel clinical"))
ggplot2::ggplot(data = d , aes(x = Model, y = elpd_diff, ymin = elpd_diff_lower, ymax = elpd_diff_upper)) + geom_errorbar(width=0.2, size=1, color="blue")  + geom_point( mapping=aes(x=Model, y=elpd_diff), size=4, shape=21, fill="white") + ylim(c( -200, 10)) + hline_0()+labs(title = "Leave One Out cross-validation error", subtitle = "Comparision with null model", y = "Leave One Out cross-validation error")
