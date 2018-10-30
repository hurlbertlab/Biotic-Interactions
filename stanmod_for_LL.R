#r mod for LL
library(rstanarm)

occuenv = read.csv("all_expected_pres.csv", header = TRUE) 
occumatrix$c_s = scale(occumatrix$comp_scaled, scale = T, center = T)
occumatrix$abTemp=abs(occumatrix$zTemp)
occumatrix$abElev=abs(occumatrix$zElev)
occumatrix$abPrecip=abs(occumatrix$zPrecip)
occumatrix$abNDVI=abs(occumatrix$zNDVI)

# using equation species sum*Focal occ to get success and failure for binomial anlaysis
occumatrix$sp_success = 15 * occumatrix$FocalOcc
occumatrix$sp_fail = 15 * (1 - occumatrix$FocalOcc)

occumatrix$sp_success = as.integer(occumatrix$sp_success)
occumatrix$sp_fail = as.integer(occumatrix$sp_fail)
mm <- stan_glmer(cbind(sp_success, sp_fail) ~ c_s + 
       abTemp + abElev + abPrecip + abNDVI + (c_s + abTemp + abElev + abPrecip + abNDVI|FocalAOU), family = binomial(link = logit), data = occumatrix, iter = 10000, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1))
write.csv(summary(mm), row.names= TRUE)
