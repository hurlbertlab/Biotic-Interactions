#### Bayesian of all matrices not just subset ####
library(brms)
library(rstudioapi)
mmslope_all <- brm(sp_success | trials(sp_success+sp_fail) ~ c_s +  abTemp + abElev + abPrecip + abNDVI + (c_s + abTemp + abElev + abPrecip + abNDVI|FocalAOU), family = binomial(link = logit), data = occumatrix_all , cores = 4, chains=4, iter=5000,warmup=2000,control = list(max_treedepth = 15),set_prior("lkj(1)", class = "cor"))

save(mmslope_all, filename ="mmslope_all.rda")
