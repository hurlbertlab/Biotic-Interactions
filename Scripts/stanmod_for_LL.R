#r mod for LL
library(rstanarm)

occumatrix = read.csv("occumatrix.csv", header = TRUE) 


mmslope <- stan_glmer(cbind(sp_success, sp_fail) ~ c_s + 
abTemp + abElev + abPrecip + abNDVI + (c_s + abTemp + abElev + abPrecip + abNDVI|FocalAOU), family = binomial(link = logit), data = occumatrix,warmup=5000,iter=10000, chains = 4, cores = 4, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1))
save(mmslope, filename ="mmslope.rda")
write.csv(mmslope, "mmslope.csv", row.names = FALSE)