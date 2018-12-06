library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggExtra)
library(rstanarm)
library(boot)

# read in files created in data cleaning script
temp_occ = read.csv("data/bbs_sub1.csv", header=TRUE) # BBS occ script
centroid=read.csv("data/centroid.csv", header=TRUE) # GIS script
occuenv= read.csv("data/all_expected_pres.csv", header = TRUE) # Data cleaning script
subsetocc = read.csv('data/subsetocc.csv', header = T) # Hurlbert Lab
tax_code = read.csv("data/Tax_AOU_Alpha.csv", header = TRUE) # Hurlbert Lab
bbs_abun = read.csv("data/bbs_abun.csv", header =TRUE)
nsw = read.csv("data/new_spec_weights.csv", header = TRUE)
shapefile_areas = read.csv("data/shapefile_areas.csv", header =TRUE)

#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
subsetocc$AOU[subsetocc$AOU == 4810] = 4812
temp_occ$Aou[temp_occ$Aou == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4123] = 4120

# have to omit focal/competitors with no overlap
shapefile_areas = na.omit(shapefile_areas)
AOU = read.csv("data/Bird_Taxonomy.csv", header = TRUE) # taxonomy data
bbs = read.csv('data/bbs_abun.csv', header = T) # BBS abundance data - from Hurlbert Lab 
maincomp = read.csv("data/shapefileareas_w_comp.csv", header = TRUE)
maincomp1 = subset(maincomp, mainCompetitor == 1)
maincomp1.5 = left_join(maincomp1, tax_code, by = c("focalAOU" = "AOU_OUT"))
maincomp1.5$Focal_Common_Name = maincomp1.5$PRIMARY_COM_NAME
maincomp2 = left_join(maincomp1.5[,c("Focal", "focalAOU","Competitor", "compAOU", "FocalArea", "CompArea", "area_overlap", "PropOverlap", "Focal_Common_Name")], tax_code, by = c("compAOU" = "AOU_OUT"))


# make Table S1 #
focalspp = read.csv("data/focal spp.csv", header = TRUE)
focalspp1 = left_join(nsw[,c("focalAOU","Focal", "FocalSciName","CompAOU","Competitor", "CompSciName", "Family")], maincomp[,c("focalAOU","compAOU", "mainCompetitor")], by = c("focalAOU"= "focalAOU", "CompAOU" =  "compAOU"))
# write.csv(focalspp1, "data/TableS1.csv", row.names = FALSE)

# rbind missing spp to tax code data frame
pacific = c("Pacific Wren", 7221,  "PAWR")
quail = c("California Quail", 2940,  "CAQU")
quail2 = c("Scaled Quail", 2930,  "SCQU")
lost_spp = t(data.frame(pacific, quail, quail2, stringsAsFactors = FALSE))
colnames(lost_spp) = c("PRIMARY_COM_NAME", "AOU_OUT" ,"ALPHA.CODE")
tax_code = rbind(tax_code, lost_spp)
tax_code$AOU_OUT = as.numeric(tax_code$AOU_OUT)

# rescaling all occupancy values  - odds ratio
# need to get rid of ones in order to not have infinity values 
edge_adjust = .005 
occuenv$FocalOcc_scale = (occuenv$FocalOcc * (1 - 2*edge_adjust)) + edge_adjust
# create logit transformation function, did on rescaled vals
occuenv$occ_logit =  log(occuenv$FocalOcc_scale/(1-occuenv$FocalOcc_scale)) 

##### LIN REG #######
# for loop subsetting env data to expected occurrence for focal species
envoutput = c()
envoutputa = c()
# create beta output data frame
beta_occ = c()
beta_abun = c()

l_focalspecies = inner_join(shapefile_areas, occuenv, c("focalAOU"="FocalAOU"))
subfocalspecies = unique(occuenv$FocalAOU) #[-99] excluding sage sparrow bc no good routes

for (sp in 1:length(subfocalspecies)){
  print(sp)
  temp = subset(occuenv, FocalAOU == subfocalspecies[sp])
  
  competition <- lm(temp$occ_logit ~  temp$comp_scaled)  # changes between main and all comps
  # z scores separated out for env effects (as opposed to multivariate variable)
  env_z = lm(temp$occ_logit ~ abs(zTemp) + abs(zElev) + abs(zPrecip) + abs(zNDVI), data = temp)
  # z scores separated out for env effects
  both_z = lm(temp$occ_logit ~  temp$comp_scaled + abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
  
  # abundance, not temp occ - same results?
  competition_abun <- lm(temp$FocalAbundance ~  temp$comp_scaled) 
  # z scores separated out for env effects - abundance
  env_abun = lm(temp$FocalAbundance ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  # z scores separated out for env effects - abundance
  both_abun = lm(temp$FocalAbundance ~  comp_scaled + abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  
  #variance_partitioning 
  ENV = summary(both_z)$r.squared - summary(competition)$r.squared #env only
  COMP = summary(both_z)$r.squared - summary(env_z)$r.squared #competition only
  SHARED = summary(competition)$r.squared - COMP #shared variance
  NONE = 1 - summary(both_z)$r.squared # neither variance
  sp1 = unique(temp$FocalAOU)
  sum = sum(ENV, COMP, SHARED)
  n = length(temp$stateroute)
  envoutput = rbind(envoutput, c(sp1, ENV, COMP, SHARED, NONE, n)) #, sum
  
  #variance_partitioning 
  ENVa = summary(both_abun)$r.squared - summary(competition_abun)$r.squared
  COMPa = summary(both_abun)$r.squared - summary(env_abun)$r.squared
  SHAREDa = summary(competition_abun)$r.squared - COMP
  NONEa = 1 - summary(both_abun)$r.squared
  
  sp1 = unique(temp$FocalAOU)
  
  
  envoutputa = rbind(envoutputa, c(sp1, ENVa, COMPa, SHAREDa, NONEa))
  
  if(length(unique(temp$comp_scaled[!is.na(temp$comp_scaled)])) > 2){
  # saving model output into separate data frames
  occ_comp_est = summary(competition)$coef[2,"Estimate"]
  occ_comp_p = summary(competition)$coef[2,"Pr(>|t|)"]
  occ_comp_r = summary(competition)$r.squared
  #occ_env_est = mean(summary(env_z)$coef[,"Estimate"])
  #occ_env_p = summary(env_z)$coef[2,"Pr(>|t|)"]
  occ_env_r = summary(env_z)$r.squared 
  #occ_b_est = summary(both_z)$coef[2,"Estimate"]
  occ_b_p = summary(both_z)$coef[2,"Pr(>|t|)"]
  occ_b_r = summary(both_z)$r.squared 

  abun_comp_est = summary(competition_abun)$coef[2,"Estimate"]
  abun_comp_p = summary(competition_abun)$coef[2,"Pr(>|t|)"]
  abun_comp_r = summary(competition_abun)$r.squared #using multiple rsquared
  #abun_env_est = summary(env_abun)$coef[2,"Estimate"]
  #abun_env_p = summary(env_abun)$coef[2,"Pr(>|t|)"]
  abun_env_r = summary(env_abun)$r.squared 
  #abun_both_est = summary(both_abun)$coef[2,"Estimate"]
  abun_both_p = summary(both_abun)$coef[2,"Pr(>|t|)"]
  abun_both_r = summary(both_abun)$r.squared
  
  beta_occ = rbind(beta_occ, c(sp1, occ_comp_est, occ_comp_p, occ_comp_r,occ_env_r, occ_b_p, occ_b_r))
  beta_abun = rbind(beta_abun, c(sp1, abun_comp_est, abun_comp_p, abun_comp_r, abun_env_r, abun_both_p , abun_both_r))
  } 
  else{print(FALSE)}
}         


envoutput = data.frame(envoutput)
envoutputa = data.frame(envoutputa)

names(envoutput) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE", "n")
names(envoutputa) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")

envoutput1 = merge(envoutput, tax_code, by.x = 'FocalAOU', by.y = "AOU_OUT", all.x = TRUE) 

envoutput2 = merge(envoutput, subsetocc[,c("AOU", "migclass", "Trophic.Group")], by.x='FocalAOU', by.y='AOU', all.x = TRUE)

envloc = merge(envoutput2, centroid[, c("FocalAOU", "Long", "Lat")], by = 'FocalAOU', all.x = TRUE)

#write.csv(envoutputa, "data/envoutputa.csv", row.names = FALSE)
beta_occ = data.frame(beta_occ)
names(beta_occ) = c("FocalAOU", "Competition_Est", "Competition_P", "Competition_R2", "EnvZ_R2", "BothZ_P", "BothZ_R2")
beta_occ_comp = left_join(beta_occ, maincomp2[,c("focalAOU","Focal_Common_Name", "PRIMARY_COM_NAME")], by = c("FocalAOU" = "focalAOU"))
beta_abun = data.frame(beta_abun)
names(beta_abun) = c("FocalAOU", "Competition_Est", "Competition_P", "Competition_R2", "EnvZ_R2", "BothZ_P", "BothZ_R2")
beta_abun_comp = left_join(beta_abun, maincomp2[,c("focalAOU","Focal_Common_Name", "PRIMARY_COM_NAME")], by = c("FocalAOU" = "focalAOU"))

beta_TableS4 = left_join(beta_occ, beta_abun_comp, by = "FocalAOU")

leftover_birds = left_join(envoutput, maincomp2[,c("focalAOU","Focal_Common_Name", "PRIMARY_COM_NAME")], by = c("FocalAOU" = "focalAOU"))

# write.csv(beta_TableS4, "data/beta_TableS4.csv", row.names = FALSE)

#### ---- GLM fitting  ---- ####
# add on success and failure columns by creating # of sites where birds were found
# and # of sites birds were not found from original bbs data
# occumatrix = merge(temp_occ, occuenv, by.x=c("Aou", "stateroute"),by.y=c("Species", "stateroute"))
occumatrix=occuenv
occumatrix$c_s = scale(occumatrix$comp_scaled, scale = T, center = T)
occumatrix$abTemp=abs(occumatrix$zTemp)
occumatrix$abElev=abs(occumatrix$zElev)
occumatrix$abPrecip=abs(occumatrix$zPrecip)
occumatrix$abNDVI=abs(occumatrix$zNDVI)

# using equation species sum*Focal occ to get success and failure for binomial anlaysis
occumatrix$sp_success = as.integer(15 * occumatrix$FocalOcc)
occumatrix$sp_fail = as.integer(15 * (1 - occumatrix$FocalOcc))

#### GLM of all matrices not just subset ####
# glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ c_s + 
 #   abTemp + abElev + abPrecip + abNDVI + (1|FocalAOU), family = binomial(link = logit), data = occumatrix)
# summary(glm_occ_rand_site)                                    


mmslope <- stan_glmer(cbind(sp_success, sp_fail) ~ c_s + 
       abTemp + abElev + abPrecip + abNDVI + (c_s + abTemp + abElev + abPrecip + abNDVI|FocalAOU), family = binomial(link = logit), data = occumatrix, iter = 40000, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1))
# write.csv(summary(mm), "mm_slope.csv", row.names= TRUE)

mmint <- stan_glmer(cbind(sp_success, sp_fail) ~ c_s + abTemp + abElev + abPrecip + abNDVI + (1|FocalAOU), family = binomial(link = logit), data = occumatrix, iter = 10000, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1))

# Default priors =  default_prior_test <- stan_glmer(cbind(sp_success, sp_fail) ~ c_s + 
#        abTemp + abElev + abPrecip + abNDVI + (1|FocalAOU), family = binomial(link = logit), data = occumatrix, chains =1, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1))
# The prior_summary function provides a concise summary of the priors used: prior_summary(default_prior_test)


mmslope = read.csv("data/mm_slope_full.csv", header = TRUE)
mm2 = subset(mmslope, mean > -2.52e+05)
mm2$X = gsub("b", "", mm2$X) 
mm2$X = gsub("c_s", "cs", mm2$X) 
mm2$X = gsub("FocalAOU", "", mm2$X) 
mm2$X = gsub(":", "_", mm2$X) 
mm2$X = gsub("[() ]", "", mm2$X) 
mm2$X = gsub("[[]", "", mm2$X) 
mm2$X = gsub("[]]", "", mm2$X)
mm2 = left_join(., mm2, by = )
mm2$AOU =strsplit(mm2$X,"_")

vecAOU = 1:1078
for(i in 1:length(mm2$AOU)){
  vecAOU[i] <- mm2$AOU[[i]][2]
}
mm2 = cbind(vecAOU, mm2)
mm2$vecAOU2 = as.numeric(as.character(mm2$vecAOU))
mm3 = left_join(mm2[,c("vecAOU2","X","mean" ,"mcse", "sd" ,"X2.5.","X25.","X50.","X75.","X97.5.","n_eff","Rhat")], nsw[,c("focalAOU", "Focal", "Family")], by = c("vecAOU2" = "focalAOU"))
write.csv(mm3, "data/tableS6.csv", row.names = FALSE)

mmint = read.csv("data/bayesian_sum_mod_output_full_11_14.csv", header = TRUE)
mm2 = read.csv("data/mm_slope_full.csv", header = TRUE)
modoutput2 = subset(mm2, mean > -2.52e+05)
modoutput2$X = gsub("b", "", modoutput2$X) 
modoutput2$X = gsub("(Intercept)", "", modoutput2$X) 
modoutput2$X = gsub("FocalAOU", "", modoutput2$X) 
modoutput2$X = gsub(":", "", modoutput2$X) 
modoutput2$X = gsub("", "", modoutput2$X) 


hist(modoutput2$mean, breaks = 20)

mm3 = modoutput2[7:181,] #removed global vals, only have ind spp
ggplot(data = mm3) + #geom_point() + geom_ribbon(aes(ymin = X2.5., ymax = X97.5.))
stat_density(aes(mm3$mean))

mm_fixed = modoutput2[1:6,]

occumatrix$cspred = inv.logit(occumatrix$comp_scaled * mm_fixed$mean[2] + mm_fixed$mean[1])
occumatrix$temppred = inv.logit(occumatrix$abTemp * mm_fixed$mean[3] + mm_fixed$mean[1])
occumatrix$elevpred = inv.logit(occumatrix$abElev * mm_fixed$mean[4] + mm_fixed$mean[1])
occumatrix$precippred = inv.logit(occumatrix$abPrecip * mm_fixed$mean[5] + mm_fixed$mean[1])
occumatrix$ndvipred = inv.logit(occumatrix$abNDVI * mm_fixed$mean[6] + mm_fixed$mean[1])
# temp_pred <- predict(stan_glmer(cbind(sp_success, sp_fail) ~ c_s + 
#abTemp + abElev + abPrecip + abNDVI + (1|FocalAOU), family = binomial(link = logit), data = occumatrix, iter = 10000, prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1)), newdata = data.frame(range(occumatrix$abTemp)), type = "response")
tempsub = occumatrix[occumatrix$abTemp < quantile(occumatrix$abTemp, 0.95), ]
tempsub = tempsub[,c("FocalAOU", "FocalOcc", "abTemp", "abElev", "abPrecip", "abNDVI", "temppred")]
temp = ggplot(data = tempsub, aes(x = abTemp, y = temppred)) + scale_x_continuous(limits = c(0,2))  +
  geom_point(data = tempsub, aes(x = abTemp, y = FocalOcc), shape=18, alpha = 0.05,position=position_jitter(width=0,height=.02)) + geom_line(color = "#2ca25f", lwd = 3) + theme_classic() + xlab("Temperature") + ylab("Focal Occupancy") + theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36, vjust = 2), axis.text.x=element_text(size=32), axis.text.y=element_text(size=32))
#geom_smooth(stat = "smooth", formula = FocalOcc ~ inv.logit(mm_fixed$mean)[3]*abTemp + 1, ymin = inv.logit(mm_fixed$X2.5)[3], ymax = inv.logit(mm_fixed$X97.5)[3], data = tempsub)  
ggsave("C:/Git/Biotic-Interactions/Figures/temp.pdf", height = 8, width = 12)

elevsub = occumatrix[occumatrix$abElev < quantile(occumatrix$abElev, 0.95), ]
elev = ggplot(data = elevsub, aes(x = abElev, y = elevpred)) + theme_classic()  +
  geom_point(data = elevsub, aes(x = abElev, y = FocalOcc), shape=18, alpha = 0.05,position=position_jitter(width=0,height=.02)) + geom_line(color = "#2ca25f", lwd = 3) + xlab("Elevation") + ylab("Focal Occupancy") + theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36, vjust = 2), axis.text.x=element_text(size=32), axis.text.y=element_text(size=32))
ggsave("C:/Git/Biotic-Interactions/Figures/elev.pdf", height = 8, width = 12)

precipsub = occumatrix[occumatrix$abPrecip < quantile(occumatrix$abPrecip, 0.95), ]
precip = ggplot(data = precipsub, aes(x = abPrecip, y = precippred)) + scale_x_continuous(limits = c(0,2))  + geom_point(data = elevsub, aes(x = abPrecip, y = FocalOcc), shape=18, alpha = 0.05,position=position_jitter(width=0,height=.02)) + geom_line(color = "#2ca25f", lwd = 3) + theme_classic()+ xlab("Precipitation") + ylab("Focal Occupancy") + theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36, vjust = 2), axis.text.x=element_text(size=32), axis.text.y=element_text(size=32))
ggsave("C:/Git/Biotic-Interactions/Figures/precip.pdf", height = 8, width = 12)

ndvisub = occumatrix[occumatrix$abNDVI < quantile(occumatrix$abNDVI, 0.95), ]
NDVI = ggplot(data = ndvisub, aes(x = abNDVI, y = ndvipred)) +
  geom_point(data = ndvisub, aes(x = abNDVI, y = FocalOcc), shape=18, alpha = 0.05,position=position_jitter(width=0,height=.02)) + geom_line(color = "#2ca25f", lwd = 3) +theme_classic()+ xlab("NDVI") + ylab("Focal Occupancy") + theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36, vjust = 2), axis.text.x=element_text(size=32), axis.text.y=element_text(size=32))

ggsave("C:/Git/Biotic-Interactions/Figures/ndvi.pdf", height = 8, width = 12)

cssub = occumatrix[abs(occumatrix$comp_scaled) < quantile(abs(occumatrix$comp_scaled), 0.95), ]
comp = ggplot(data = cssub, aes(x = comp_scaled, y = cspred)) + xlim(0,1) +
  geom_point(data = cssub, aes(x = comp_scaled, y = FocalOcc), shape=18, alpha = 0.05,position=position_jitter(width=0,height=.02)) + geom_line(color = "#dd1c77", lwd = 3) + theme_classic()+ xlab("Competitor Abundance") + ylab("Focal Occupancy") + theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36, vjust = 2), axis.text.x=element_text(size=32), axis.text.y=element_text(size=32))
ggsave("C:/Git/Biotic-Interactions/Figures/comp.pdf", height = 8, width = 12)
#   geom_segment(aes(x = 0, y =  inv.logit(mm_fixed$mean)[1], xend = inv.logit(mm_fixed$mean[2]), yend = 0), col = "dark green", lwd=2) 


z <- plot_grid(comp + theme(legend.position="none"),
               temp + theme(legend.position="none"),
               NDVI + theme(legend.position="none"),
               elev + theme(legend.position="none"),
               precip + theme(legend.position="none"),
               align = 'hv',
               # labels = c("A","B", "C", "D", "E"),
               nrow = 2)
ggsave("C:/Git/Biotic-Interactions/Figures/cowplotabiotic.pdf", height = 20, width = 30)

mm_fixed2 = mm_fixed[2:6,]
ggplot(mm_fixed2, aes(x = X, y = mean)) + geom_point(pch=15, size = 5, col = "dark blue") + theme_classic() + geom_hline(yintercept = 0, lty = 2, color = "red") + geom_errorbar(ymin = mm_fixed2$X2.5., ymax = mm_fixed2$X97.5., width=0.2, size=1, color="black") + scale_y_continuous(limits = c(-0.6, .1)) + xlab("Parameter Estimate") + ylab("Value") +
  scale_x_discrete("Parameter Estimate", labels = c("Elevation","NDVI","Precipitation", "Temperature", "Competitor Abundance"))+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) + 
  theme(axis.line=element_blank(),axis.text.x=element_text(size=25),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))
ggsave("C:/Git/Biotic-Interactions/Figures/abiotic_estimates.pdf", height = 8, width = 11)

#### new fig 1 ####
occ1b = occuenv %>% filter(FocalAOU == 6860|FocalAOU  == 7222|FocalAOU  == 5840) %>%
        filter(stateroute == 68015)

fig1a = ggplot(data = occuenv, aes(x = log10(FocalAbundance), y = FocalOcc)) +geom_point()+ geom_jitter(width = 0, height = 0.02)  +xlab("log10(Focal Abundance)")+ylab("Focal Occupancy") + geom_hline(yintercept = 0.5, lwd = 1, col = "dark red", lty = 2)+ geom_vline(xintercept = median(log10(occuenv$FocalAbundance)), lwd = 1, col = "dark red", lty = 2) +theme_classic() +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90), axis.text=element_text(size=12)) + theme(plot.margin = unit(c(.5,6,.5,.5),"lines")) + geom_point(data = occ1b, aes(x = log10(FocalAbundance), y = FocalOcc, color = as.factor(occ1b$Focal), size = 20)) +theme(legend.position = "none")+ scale_color_manual(breaks = c("Swamp Sparrow",  "Canada Warbler", "Winter Wren"), values=c("#e7298a","#cb181d", "#fd8d3c"), labels=c("Swamp Sparrow",  "Canada Warbler", "Winter Wren"))
#+geom_point(data = bbs_sub4, color = "red")
# +geom_text(label = occuenv$Species)
ggExtra::ggMarginal(fig1a , type = "histogram", fill = "dark gray")
ggsave("C:/Git/Biotic-Interactions/Figures/fig1.png", height = 8, width = 12)

bbs_sub3.5 = bbs_abun %>% filter(aou == 6860|aou == 7222|aou == 5840) %>%
  filter(stateroute == 68015)
bbs_sub4 = read.csv("data/bbs_route68015.csv", header = TRUE)
bbs_sub4$speciestotal[bbs_sub4$speciestotal == 0] <- NA
fig1b = ggplot(data = bbs_sub4, aes(x = year, y = speciestotal))+ geom_line(aes(color = as.factor(bbs_sub4$SpeciesName)), lwd = 1.5) + geom_point(aes(color = as.factor(bbs_sub4$SpeciesName), size = 12))+theme_classic()+ xlim(c(0, 30)) +xlab("Year")+ylab("Abundance at Route") +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90),legend.title=element_blank(), axis.text=element_text(size=16), legend.text = element_text(size = 12)) + theme(plot.margin = unit(c(.5,6,.5,.5),"lines"))  + scale_x_continuous(breaks = c(2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015))+ scale_y_continuous(breaks = c(0, 5, 10, 20, 29))+ scale_color_manual(breaks = c("Swamp Sparrow",  "Canada Warbler", "Winter Wren"), values=c("#e7298a","#cb181d", "#fd8d3c"), labels=c("Swamp Sparrow",  "Canada Warbler", "Winter Wren")) 
ggsave("C:/Git/Biotic-Interactions/Figures/1b.pdf", height = 7, width = 12)

fig1 = plot_grid(fig1a + theme(legend.position="none"),
               fig1b + theme(legend.position="none"), 
               labels = c("A","B"),
               align = 'h', 
               rel_widths = c(1, 1.3))


bbs_sub4 = bbs_abun %>%
  filter(stateroute == 68015)  %>% left_join(tax_code, by = c("aou" = "AOU_OUT"))%>% filter(aou %in%  c(4980, 6280
, 6130, 6160, 5980))
bbs_sub4$speciestotal[bbs_sub4$speciestotal == 0] <- NA
ggplot(data = bbs_sub4, aes(x = year, y = speciestotal))+ geom_line(aes(color = as.factor(bbs_sub4$PRIMARY_COM_NAME)), lwd = 1.5) + theme_classic() +xlab("Year")+ylab("Abundance at Route") +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90),legend.title=element_blank(), axis.text=element_text(size=16), legend.text = element_text(size = 12)) + theme(plot.margin = unit(c(.5,6,.5,.5),"lines"))  + scale_x_continuous(breaks = c(2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015))+ scale_color_manual(breaks = c("Bank Swallow",  "Barn Swallow", "Indigo Bunting", "Red-winged Blackbird", "Yellow-throated Vireo"), values=c("#e41a1c","#377eb8", "#4daf4a", "#984ea3", '#ff7f00', "#fc8d62")) 
ggsave("C:/Git/Biotic-Interactions/Figures/forLB.pdf", height = 7, width = 12)


ggplot(envoutput1, aes(x = ENV, y = COMP)) + geom_point() + geom_abline(intercept = 0, slope = 1, col = "navy", lwd = 1.25)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=20, hjust = 1, vjust = 0.5), legend.position = c(0.2,0.9))

##### Variance Partitioning Plot #####
envloc$EW <- 0
envloc$EW[envloc$Long > -98.583333] <- 1 ## mid point of USA
# 1 = East
envloc$COMPSC = envloc$COMP/(envloc$COMP+envloc$ENV)

envloc1 = merge(envloc, occumatrix[,c("FocalAOU", "Family")], by.x = "FocalAOU", by.y="FocalAOU", all.y = FALSE)
# cutting down to unique subset
envloc1 = unique(envloc1)
envloc1 = merge(envloc1, tax_code[,c("AOU_OUT", "ALPHA.CODE")], by.x = "FocalAOU", by.y = "AOU_OUT", all.x=TRUE)
envloc1$ALPHA.CODE = as.character(envloc1$ALPHA.CODE)
envloc1$ALPHA.CODE[envloc1$FocalAOU == 2920] <- 'MOUQ' #Mountain Quail
envloc1$ALPHA.CODE[envloc1$FocalAOU == 6720] <- 'PAWA' #Palm Warbler
envloc1$ALPHA.CODE = as.factor(envloc1$ALPHA.CODE)

nrank = envloc1 %>% 
  dplyr::mutate(rank = row_number(-COMP))# change here for comp
nrank$NONE = 0
envflip = tidyr::gather(nrank, "Type", "value", 2:5)
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)])
envflip = dplyr::arrange(envflip,rank)

# envflip = merge(envflip, envloc[,c("FocalAOU", "EW")], by = "FocalAOU")

envrank = envflip %>% 
  dplyr::group_by(Type == 'COMP') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "COMP") # change here for comp

### CREATE LABEL OF FAMilY ########
envrank$Fam_abbrev = envrank$Family
envrank$Fam_abbrev = gsub('Emberizidae','E', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Turdidae','Tu', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Fringillidae','F', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Tyrannidae','Ty', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Mimidae','M', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Vireonidae','V', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Aegithalidae','A', envrank$Fam_abbrev)                        
envrank$Fam_abbrev = gsub('Corvidae','Co', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Timaliidae','Ti', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Troglodytidae','T', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Cuculidae','Cu', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Icteridae','I', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Picidae','Pi', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Motacillidae','M', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Columbidae','Cl', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Trochilidae','Tr', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Cardinalidae','Ca', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Paridae','Pa', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Sylviidae','Sy', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Parulidae','P', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Sittidae','Si', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Regulidae','R', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Hirundinidae','H', envrank$Fam_abbrev)
envrank$Fam_abbrev = gsub('Certhiidae','Ce', envrank$Fam_abbrev)


envrank$Fam_abbrevf = as.factor(as.character(envrank$Fam_abbrev))
envrank$Fam_abbrevf = gsub('E','#000000', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Tu','#a6bddb', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('F','#67a9cf', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Tr','#048691', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Ty','#9ecae1', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Ti','#02818a', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('V','#016c59', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('A','#014636', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Co','#081d58', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Cu','#253494', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('I','#225ea8', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Pi','#1d91c0', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('M','#41b6c4', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('O','#7f7fff', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Cl','#0000ff', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Ca','#016c59', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Pa','#02818a', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('Sy','#014636', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('P','#3690c0', envrank$Fam_abbrevf)
envrank$Fam_abbrevf = gsub('T','#0080ff', envrank$Fam_abbrevf) 
famlabel= envrank$Fam_abbrev
####### OTHER LABEL ######
envrank$mig_abbrev = envrank$migclass
envrank$mig_abbrev = gsub("neotrop", 'L', envrank$mig_abbrev)
envrank$mig_abbrev = gsub("resid", 'R', envrank$mig_abbrev)
envrank$mig_abbrev = gsub("short", 'S', envrank$mig_abbrev)
envrank$mig_abbrevf = as.factor(as.character(envrank$mig_abbrev))

envrank$mig_abbrevf = gsub('L','#bae4b3', envrank$mig_abbrevf)
envrank$mig_abbrevf = gsub('R','#31a354', envrank$mig_abbrevf)
envrank$mig_abbrevf = gsub('S','#006d2c', envrank$mig_abbrevf)
miglabel= envrank$mig_abbrev

envrank$trophlabel = envrank$Trophic.Group
envrank$trophlabel = gsub("frugivore", 'F', envrank$trophlabel)
envrank$trophlabel = gsub("granivore", 'G', envrank$trophlabel)
envrank$trophlabel = gsub("herbivore", 'H', envrank$trophlabel)
envrank$trophlabel = gsub("insct/om", 'X', envrank$trophlabel)
envrank$trophlabel = gsub("insectivore", 'I', envrank$trophlabel)
envrank$trophlabel = gsub("nectarivore", 'N', envrank$trophlabel)
envrank$trophlabel = gsub("omnivore", 'O', envrank$trophlabel)
envrank$trophlabelf = as.factor(as.character(envrank$trophlabel))

envrank$trophlabelf = gsub('F','#fbb4b9', envrank$trophlabelf)
envrank$trophlabelf = gsub('G','#f768a1', envrank$trophlabelf)
envrank$trophlabelf = gsub('H','#c51b8a', envrank$trophlabelf)
envrank$trophlabelf = gsub('X','#7a0177', envrank$trophlabelf)
envrank$trophlabelf = gsub('I','#dd1c77', envrank$trophlabelf)
envrank$trophlabelf = gsub('N','#ce1256', envrank$trophlabelf)
envrank$trophlabelf = gsub('O','#67001f', envrank$trophlabelf)

envrank$EW[envrank$EW == 1] <- "E"
envrank$EW[envrank$EW == 0] <- "W" 
###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE","SHARED", "ENV","COMP"),ordered = TRUE)
envflip$value = abs(envflip$value)
# Plot with ENV ranked in decreasing order - had to flip everything to plot right
t = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE","SHARED", "ENV","COMP")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#2ca25f","#dd1c77"), labels=c("","Shared Variance","Environment", "Competition")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30, angle=90),legend.title=element_blank(), legend.text=element_text(size=28, hjust = 1, vjust = 0.5), legend.position = c(0.5,.8)) + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))

tt = t + annotate("text", x = 1:175, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=0.5,hjust = 0.8, color = "black") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks = c(0,0.2,0.4,0.6, 0.8))

# scales::pretty_breaks()(0:1)

#+ annotate("text", x = 1:104, y = -.08, label = envrank$mig_abbrev, size=6,vjust=0.5, color = envrank$mig_abbrevf, fontface =2) + annotate("text", x = 1:104, y = -.1, label = envrank$trophlabel, size=6,vjust=0.5, color = envrank$trophlabelf, fontface =2) + annotate("text", x = 1:104, y = -.12, label = envrank$EW, angle=90,size=6,vjust=0.5, color = "black", fontface =2)+ annotate("text", x = 1:104, y = -.06, label = envrank$Fam_abbrev, size=6,vjust=0.5, color = "black", fontface =2) 

plot(tt)

ggsave("Figures/barplotc.pdf", height = 35, width = 48)

#### top 10
maincomp = read.csv("data/shapefile_areas_w_comp.csv", header = TRUE)
maincomp2 = subset(maincomp, mainCompetitor == 1)
envflip_labs = subset(envflip_sub, Type == "COMP")
envflip_sub2 = left_join(envflip_labs, maincomp2, by = c("FocalAOU" = "focalAOU"))
envflip_sub2.5 = left_join(envflip_sub2, tax_code[, c("AOU_OUT", "PRIMARY_COM_NAME")], by = c("FocalAOU" = "AOU_OUT"))
envflip_sub3 = unique(left_join(envflip_sub2.5, tax_code, by = c("compAOU" = "AOU_OUT")))


envflip_sub = envflip[1:60,]
c = ggplot(data=envflip_sub, aes(factor(rank), y=abs(value), fill=factor(Type, levels = c("NONE","SHARED", "ENV","COMP")))) + geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#2ca25f","#dd1c77"), labels=c("","Shared Variance","Environment", "Competition")) +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90),legend.title=element_blank(), legend.text=element_text(size=22, hjust = 1, vjust = 0.5), legend.position = c(.8,.8)) + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))+ theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 24))   + annotate("text", x = 1:15, y = 0.3, label = envflip_sub3$PRIMARY_COM_NAME.y, angle=90,size=8,vjust=0.5,hjust = 1.4, color = "white", fontface = "bold") + ylim(0, 0.8) + annotate("text", x = 1:15, y = -.01, label = envrank$ALPHA.CODE[1:15], angle=90,size=6,vjust=0.5,hjust = 1, color = "black") #+ scale_y_continuous(breaks = c(0,0.52,0.5,0.75, 1))
ggsave("Figures/barplotc_sub.pdf", height = 18, width = 16)

geom_histogram(envoutput$ENV + envoutput$SHARED)
ggplot(envoutput, aes(x = ENV+SHARED)) + geom_histogram(binwidth = 0.05, fill = "#2ca25f") + xlab("Environment and Shared Variance Explained") + ylab("Frequency")
hist(envoutput$COMP + envoutput$SHARED)

#### ENV ####
nrank = envloc1 %>% 
  dplyr::mutate(rank = row_number(-ENV))# change here for comp
nrank$NONE = 0
envflip = tidyr::gather(nrank, "Type", "value", 2:5)
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)])
envflip = dplyr::arrange(envflip,rank)

# envflip = merge(envflip, envloc[,c("FocalAOU", "EW")], by = "FocalAOU")

envrank = envflip %>% 
  dplyr::group_by(Type == 'ENV') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "ENV") # change here for comp

# CREATE LABEL DF FAMilY 
# OTHER LABEL #
envrank$mig_abbrev = envrank$migclass
envrank$mig_abbrev = gsub("neotrop", 'L', envrank$mig_abbrev)
envrank$mig_abbrev = gsub("resid", 'R', envrank$mig_abbrev)
envrank$mig_abbrev = gsub("short", 'S', envrank$mig_abbrev)
envrank$mig_abbrevf = as.factor(as.character(envrank$mig_abbrev))

envrank$mig_abbrevf = gsub('L','#bae4b3', envrank$mig_abbrevf)
envrank$mig_abbrevf = gsub('R','#31a354', envrank$mig_abbrevf)
envrank$mig_abbrevf = gsub('S','#006d2c', envrank$mig_abbrevf)
miglabel= envrank$mig_abbrev

envrank$trophlabel = envrank$Trophic.Group
envrank$trophlabel = gsub("frugivore", 'F', envrank$trophlabel)
envrank$trophlabel = gsub("granivore", 'G', envrank$trophlabel)
envrank$trophlabel = gsub("herbivore", 'H', envrank$trophlabel)
envrank$trophlabel = gsub("insct/om", 'X', envrank$trophlabel)
envrank$trophlabel = gsub("insectivore", 'I', envrank$trophlabel)
envrank$trophlabel = gsub("nectarivore", 'N', envrank$trophlabel)
envrank$trophlabel = gsub("omnivore", 'O', envrank$trophlabel)
envrank$trophlabelf = as.factor(as.character(envrank$trophlabel))

envrank$trophlabelf = gsub('F','#fbb4b9', envrank$trophlabelf)
envrank$trophlabelf = gsub('G','#f768a1', envrank$trophlabelf)
envrank$trophlabelf = gsub('H','#c51b8a', envrank$trophlabelf)
envrank$trophlabelf = gsub('X','#7a0177', envrank$trophlabelf)
envrank$trophlabelf = gsub('I','#dd1c77', envrank$trophlabelf)
envrank$trophlabelf = gsub('N','#ce1256', envrank$trophlabelf)
envrank$trophlabelf = gsub('O','#67001f', envrank$trophlabelf)

envrank$EW[envrank$EW == 1] <- "E"
envrank$EW[envrank$EW == 0] <- "W" 
###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE", "SHARED","COMP","ENV"),ordered = TRUE)
envflip$value = abs(envflip$value)
# Plot with ENV ranked in decreasing order - had to flip everything to plot right
e = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE","SHARED","COMP", "ENV")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#dd1c77","#2ca25f"), labels=c("","Shared Variance", "Competition","Environment")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30, angle=90),legend.title=element_blank(), legend.text=element_text(size=50, hjust = 1, vjust = 0.5), legend.position = c(0.5,0.9)) # + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))

ee = e + annotate("text", x = 1:175, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=1,hjust = 0.8, color = "black") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks = c(0,0.2,0.4,0.6, 0.8))

ee

ggsave("C:/Git/Biotic-Interactions/Figures/barplot.pdf", height = 25, width = 48)


z <- plot_grid(tt+ theme(legend.position="top"),
               ee + theme(legend.position="none"),
               nrow = 2,
               align = 'hv',
               labels = c("A","B"),
               label_size = 36,
               hjust = -6)
ggsave("C:/Git/Biotic-Interactions/Figures/barplotboth.pdf", height = 25, width = 36)

env_sub = envflip[1:60,]
env_labs = unique(subset(env_sub, Type == "ENV"))
env_sub2 = left_join(env_labs, maincomp2, by = c("FocalAOU" = "focalAOU"))
env_sub2.5 = left_join(env_sub2, tax_code[, c("AOU_OUT", "PRIMARY_COM_NAME")], by = c("FocalAOU" = "AOU_OUT"))
env_sub3 = unique(left_join(env_sub2.5, tax_code, by = c("compAOU" = "AOU_OUT")))

w = ggplot(data=env_sub, aes(factor(rank), y=abs(value), fill=factor(Type, levels = c("NONE","SHARED","COMP", "ENV")))) + geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#dd1c77","#2ca25f"), labels=c("","Shared Variance", "Competition","Environment")) +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90),legend.title=element_blank(), legend.text=element_text(size=22, hjust = 1, vjust = 0.5), legend.position = c(.8,.8)) + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))+ theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 24)) + annotate("text", x = 1:15, y = -.01, label = envrank$ALPHA.CODE[1:15], angle=90,size=6,vjust=0.5,hjust = 1, color = "black") +ylim(0,0.8)# + scale_y_continuous(breaks = c(0,0.2,0.4,0.6, 0.8))
ggsave("Figures/barplote_sub.pdf", height = 18, width = 16)


plot_grid(c+ theme(legend.position="none"),
          w + theme(legend.position="none"),
          align = 'hv')
ggsave("C:/Git/Biotic-Interactions/Figures/subbarplotboth.pdf", height = 10, width = 40)
##################### TRAITS Model ####################################
logit = function(x) log(x/(1-x))

env_lm = subset(envflip, Type == 'ENV')

# env_traits = lm(logit(value) ~ Trophic.Group + migclass + EW, data = env_lm)
# summary(env_traits) 

comp_lm = subset(envflip, Type == 'COMP')

comp_traits = lm(COMPSC ~ Trophic.Group + migclass + Lat + Long, data = comp_lm, weights = n)
summary(comp_traits) 

env_sum = subset(envflip, Type != 'NONE')
total = env_sum %>% 
  dplyr::group_by(FocalAOU) %>%
  summarise(sum(value))

env_sum$edgeval = (env_sum$value * (1 - 2*edge_adjust)) + edge_adjust
total_traits = lm(logit(edgeval) ~ Trophic.Group + migclass + EW, data = env_sum)
summary(total_traits)

env_traits = lm(logit(value) ~ EW, data = env_lm)
anova(env_traits)


# creating env traits model to compare to comp and weighted traits mods
env_cont = merge(env_lm, shapefile_areas, by.x = "FocalAOU",by.y = "focalAOU")
env_cont2 = merge(env_cont, unique(occuenv[,c("FocalAOU", "Mean.Temp","Mean.Precip","Mean.Elev","Mean.NDVI")]), by.x = "FocalAOU", by.y = "FocalAOU")
econt = lm(COMPSC ~ FocalArea  + area_overlap + Mean.Temp + Mean.Precip + Mean.Elev + Mean.NDVI , data = env_cont2, weights = n)
env_est = summary(econt)$coef[,"Estimate"]
colname = c("Intercept","FocalArea", "area_overlap","Mean.Temp","Mean.Precip","Mean.Elev", "Mean.NDVI")
env = data.frame(colname, env_est)
env$env_lower =  as.vector(summary(econt)$coefficients[,"Estimate"]) - as.vector(summary(econt)$coef[,"Std. Error"])
env$env_upper = as.vector(summary(econt)$coefficients[,"Estimate"]) + as.vector(summary(econt)$coef[,"Std. Error"])

env_trait_rank = env %>% 
  dplyr::mutate(rank = row_number(-env_est)) 
env_trait_rank2 <- env_trait_rank[order(env_trait_rank$rank),]

# making comp mod - not working.
colname = c("Intercept","FocalArea","area_overlap", "Mean.Temp","Mean.Precip","Mean.Elev","Mean.NDVI")

env_trait_rank = env %>% 
  dplyr::mutate(rank = row_number(-env_est)) 
env_trait_rank2 <- env_trait_rank[order(env_trait_rank$rank),]
env_trait_rank2$colname = factor(env_trait_rank2$colname,
                              levels = c("Mean.NDVI","Mean.Temp","Mean.Elev","FocalArea","area_overlap", "Mean.Precip","Intercept"),ordered = TRUE)

ggplot(env_trait_rank2, aes(colname, env_est)) + geom_point(pch=15, size = 5, col = "dark blue") + 
  geom_errorbar(data=env_trait_rank2, mapping=aes(ymin=env_lower, ymax=env_upper), width=0.2, size=1, color="black") + ylab(bquote("Competitor R"^"2")) + xlab("Env") + theme_classic() + theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))  + 
  theme(axis.line=element_blank(),axis.text.x=element_text(size=25),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + 
  guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))



#column names to manipulate in plot
colname = c("Intercept","Sum Overlap","Temp","Precip","Elev","NDVI","Resident","Short", "Insct/Om","Insectivore","Nectarivore","Omnivore")
# this is the trait mod scaled by comp/env. there are > 183 rows bc of the competitors (FocalArea, area_overalp)
# trait_mod_scale = lm(COMPSC ~ log10(sum_overlap) + Mean.Temp + Mean.Precip + Mean.Elev + Mean.NDVI + migclass + Trophic.Group, data = comp_cont4)
comp_cont4 = comp_lm
comp_cont4 = filter(comp_cont4, Trophic.Group != "nectarivore" & Trophic.Group != "herbivore")
colname = c("Granivore", "Insectivore/\nOmnivore","Insectivore","Omnivore")
trait_mod_scale = lm(COMPSC ~ Trophic.Group, data = comp_cont4, weights = n)
scaled_est1 = summary(trait_mod_scale)$coef[,"Estimate"]
scaled_est2 = c(scaled_est1[1], scaled_est1[2:4] + scaled_est1[1])
scaled_est = data.frame(colname, scaled_est2)
scaled_est$scaled_lower =  as.vector(scaled_est$scaled_est) - as.vector(summary(trait_mod_scale)$coef[,"Std. Error"])

scaled_est$scaled_upper = as.vector(scaled_est$scaled_est) + as.vector(summary(trait_mod_scale)$coef[,"Std. Error"])

scaled_rank = scaled_est %>% 
  dplyr::mutate(rank = row_number(-scaled_est2)) 
scaled_rank2 <- scaled_rank[order(scaled_rank$rank),]
scaled_rank2$colname = factor(scaled_rank2$colname,
       levels = c("Insectivore","Insectivore/\nOmnivore","Omnivore","Granivore","Herbivore"),ordered = TRUE)

troph = ggplot(scaled_rank2, aes(colname, scaled_est2)) + geom_point(pch=16, size = 10, col = "dark blue") + 
  geom_errorbar(data=scaled_rank2, mapping=aes(ymin=scaled_lower, ymax=scaled_upper), width=0.2, size=1, color="black") + ylab(bquote("R"["c"])) + xlab("Trophic Group") + theme_classic() + theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) + ylim(-.05,0.5) + 
  theme(axis.line=element_blank(),axis.text.x=element_text(size=25),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + 
  guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))
ggsave("C:/Git/Biotic-Interactions/Figures/traitestimateplot.pdf", height = 8, width = 12)

# ggplot(comp_cont4, aes(as.factor(Trophic.Group), COMPSC)) + geom_violin(linetype = "blank", aes(fill = factor(comp_cont4$Trophic.Group)))


#### mig mod ####
#column names to manipulate in plot
colname = c("Neotropical", "Resident" ,  "Short-distance")
trait_mod_scale = lm(COMPSC ~ migclass, data = comp_cont4, weights = n)
scaled_est1 = summary(trait_mod_scale)$coef[,"Estimate"]
scaled_est2 = c(scaled_est1[1], scaled_est1[2:3] + scaled_est1[1])
scaled_est = data.frame(colname, scaled_est2)
scaled_est$scaled_lower =  as.vector(scaled_est$scaled_est) - as.vector(summary(trait_mod_scale)$coef[,"Std. Error"])
scaled_est$scaled_upper = as.vector(scaled_est$scaled_est) + as.vector(summary(trait_mod_scale)$coef[,"Std. Error"])


scaled_rank = scaled_est %>% 
  dplyr::mutate(rank = row_number(-scaled_est2)) 
scaled_rank2 <- scaled_rank[order(scaled_rank$rank),]
scaled_rank2$colname = factor(scaled_rank2$colname,
                              levels = c("Neotropical",  "Short-distance",  "Resident"),ordered = TRUE)

tukeys = aov(lm(COMPSC ~ migclass, data = comp_cont4, weights = n))
TukeyHSD(tukeys)
summary(tukeys)

mig = ggplot(scaled_rank2, aes(colname, scaled_est2)) + geom_point(pch=16, size = 10, col = "dark blue") + 
  geom_errorbar(data=scaled_rank2, mapping=aes(ymin=scaled_lower, ymax=scaled_upper), width=0.2, size=1, color="black") + ylab(bquote("R"["c"])) + xlab("Migratory Group") + theme_classic()+ ylim(0,0.5)+ theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) + 
  theme(axis.line=element_blank(),axis.text.x=element_text(size=25),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + 
  guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))
ggsave("C:/Git/Biotic-Interactions/Figures/traitestimate_mig.pdf", height = 8, width = 12)

suppl = merge(env_lm, nsw[,c("CompAOU", "focalAOU", "Competitor", "Focal")], by.x = "FocalAOU", by.y = "focalAOU")
# write.csv(suppl, "data/suppl_table.csv", row.names = FALSE)
# anova of traits
cor.test(envoutput$ENV, envoutputa$ENV)

plot_grid(troph + theme(legend.position="none"),
          mig + theme(legend.position="none"),
          align = 'hv')
ggsave("C:/Git/Biotic-Interactions/Figures/est_mods.pdf", height = 8, width = 18)

###### Figure 4 #####
# R2 plot - lm in ggplot
# X = occupancy, Y = abundance
R2plot = merge(envoutput, envoutputa, by = "FocalAOU")

tomerge = c()
for (s in subfocalspecies) {
  spsub = subset(R2plot,FocalAOU == s)
  print(s)
  total.x = sum(spsub$COMP.x + spsub$ENV.x + spsub$SHARED.x)
  total.y = sum(spsub$COMP.y + spsub$ENV.y + spsub$SHARED.y)
  tomerge = rbind(tomerge, c(s, total.x, total.y))
}
tomerge = data.frame(tomerge)
names(tomerge) = c("FocalAOU","Total.x", "Total.y")

# R2 plot
R2plot2 = merge(R2plot, tomerge, by = "FocalAOU")
R2plot2$violin_env = R2plot2$ENV.x + R2plot2$SHARED.x
R2plot2$violin_comp = R2plot2$COMP.x + R2plot2$SHARED.x
R2plot2$violin_total = R2plot2$ENV.x + R2plot2$COMP.x + R2plot2$SHARED.x

# need to change the slopes
cols = c("Competition" ="#dd1c77","Environment" = "#2ca25f","Total" = "dark gray")
r1 = ggplot(R2plot2, aes(x = COMP.x, y = COMP.y, col = "Competition")) +theme_classic()+ theme(axis.title.x=element_text(size=26),axis.title.y=element_text(size=26, angle=90)) + xlab(bquote("Occupancy R"^"2")) + ylab(bquote("Abundance R"^"2"))+
  geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash", lwd =2.5) +
      geom_point(data = R2plot2, aes(x = ENV.x, y = ENV.y, col = "Environment"), shape = 16, cex =4, stroke = 1)+geom_smooth(data = R2plot2, aes(x = ENV.x, y = ENV.y), method='lm', se=FALSE, col="#2ca25f",linetype="dotdash", lwd = 2.5) +
      geom_point(data = R2plot2, aes(Total.x,Total.y, col = "Total"), shape = 3, cex =5, stroke = 1)+geom_smooth(data = R2plot2, aes(x =Total.x, y = Total.y), method='lm', se=FALSE, col="dark gray",linetype="dotdash", lwd =2.5) +ylim(c(0, 0.8))+ xlim(c(0, 0.8))+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=24), legend.position = c(0.2,0.9))
ggsave("C:/Git/Biotic-Interactions/Figures/occvabun_lines.pdf", height = 8, width = 12)

R2plot2$occdiff = R2plot2$COMP.x - R2plot2$ENV.x
R2plot2$abundiff = R2plot2$COMP.y - R2plot2$ENV.y
R2plot2$totaldiff = R2plot2$abundiff - R2plot2$occdiff

r2 = ggplot(R2plot2, aes(x = occdiff, y = abundiff)) +theme_classic()+ geom_abline(intercept = 0, slope = 0, col = "black", lwd = 1.25, lty = "dashed")+ylim(c(-0.4, 0.6)) + geom_vline(xintercept = 0, col = "black", lwd = 1.25, lty = "dashed")+ geom_abline(intercept = 0, slope = 1, col = "navy", lwd = 1.25)+ theme(axis.title.x=element_text(size=26),axis.title.y=element_text(size=26)) + xlab("")+ ylab("") + geom_point(col = "black", shape=16, size = 3)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20)) 
#+ annotate("text", x = -.3, y = 0.5, label = "Abundance predicts \nmore competition") + annotate("text", x = 0.4, y = -0.3, label = "Occupancy predicts \nmore environment")
smooth_vals = predict(loess(COMP.x~COMP.y,R2plot2), R2plot2$COMP.y)
summary(lm(ENV.x~ENV.y,data = R2plot2))


p2 = plot_grid(r1,
               r2 + theme(legend.position="none"), 
               labels = c("A","B"),
               label_size = 35,
               align = 'hv')
ggsave("Figures/Figure4A_B.pdf", height = 10, width = 20)


ggplot(envoutput, aes(x = ENV, y = COMP)) +theme_classic()+ theme(axis.title.x=element_text(size=26),axis.title.y=element_text(size=26, angle=90)) + xlab(bquote("Environment R"^"2")) + ylab(bquote("Competitor R"^"2")) + geom_point(cex =4, shape=16)+geom_smooth(method='lm', se=FALSE, col="black",linetype="dotdash")+geom_abline(intercept = 0, slope = 1, col = "navy", lwd = 1.25)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=20, hjust = 1, vjust = 0.5), legend.position = c(0.2,0.9))



#### Figure 1 violin plots ####
# need to chance comp_scaled to all_comp scaled in envouputput loop for Fig 1B
R2plot2$COMPSC = R2plot2$COMP.x/(R2plot2$COMP.x+R2plot2$ENV.x)
R2violin.5 = left_join(R2plot2[,c("FocalAOU", "violin_env","violin_comp","violin_total")], envloc[,c("FocalAOU", "COMPSC")], by = c("FocalAOU" = "FocalAOU"))
R2violin = gather(R2violin.5, "type", "Rval", 2:5)

R2violin$type = factor(R2violin$type,
                              levels = c("violin_comp","violin_env", "violin_total","COMPSC"),ordered = TRUE)

ggplot(R2violin, aes(as.factor(type), Rval)) + geom_violin(linetype = "blank", aes(fill = factor(R2violin$type))) + xlab("") + ylab(bquote("Variance Explained"))+scale_fill_manual(values=c("#dd1c77","#2ca25f", "grey", "#636363"), labels=c("Competition","Environment", "Total Variance", "Scaled \nCompetition")) + theme_classic()+theme(axis.title.x=element_text(size=30, angle = 180),axis.title.y=element_text(size=30, vjust = 4))+scale_y_continuous(limits = c(0, 1)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25, angle = 90),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))  + stat_summary(aes(group=factor(R2violin$type)), fun.y=median, geom="point",fill="black", shape=21, size=3, position = position_dodge(width = .9)) 
# , sec.axis = sec_axis(~ . *1/1, name = "Variance Ratio")
ggsave("Figures/violin_mains.pdf", height = 8, width = 12)

# r2 plot for main vs all competitors
# envall = read.csv("data/envoutput_all.csv", header = TRUE) NEEDS TO BE CHANGED
mainvall = merge(envoutput, envall, by = "FocalAOU")
mainvall$total.x = mainvall$ENV.x + mainvall$COMP.x + mainvall$SHARED.x
mainvall$total.y = mainvall$ENV.y + mainvall$COMP.y + mainvall$SHARED.y

ggplot(mainvall, aes(x = COMP.y, y = COMP.x)) +geom_text(aes(label = mainvall$FocalAOU))+theme_bw()+ theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16, angle=90)) + xlab("main R2") + ylab("all R2") +geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash") + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25)
ggsave("C:/Git/Biotic-Interactions/Figures/mainvallcomp.pdf")


ggplot(mainvall, aes(x = COMP.y, y = COMP.x)) + geom_point(col = mainvall, shape=16)+theme_bw()+ theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16, angle=90)) + xlab("main R2") + ylab("all R2") +geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash") + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25)


foo = merge(mainvall, focal_competitor_table, by = "FocalAOU")
foo$difference = foo$COMP.y - foo$COMP.x ### y = all competitors, x = occupancy





##### Figure 6 non-competitor comparison ######
noncompdf = occuenv[,c("FocalAOU", "stateroute", "FocalOcc", "FocalAbundance", "Family", "FocalOcc_scale", "occ_logit")]
subfocspecies = unique(noncompdf$FocalAOU)
noncomps = nsw[,c("CompAOU", "Family")]
noncomps = unique(noncomps)


noncomps_output = c()
for (sp in subfocspecies){
  FocalAOU = sp
  temp = subset(noncompdf, FocalAOU == sp) 
  tempfam = unique(as.character(temp$Family))
  if(nrow(temp) > 0){
    ncomps = dplyr::filter(noncomps, Family != tempfam) %>%
      dplyr::select(CompAOU) %>% 
      unlist()
    comps = unique(noncomps$CompAOU)
    comps = subset(comps, !comps %in% sp)
    for(co in comps){
      mergespp = subset(bbs, aou == co) %>% 
        group_by(stateroute, aou) %>%
        summarize(meancomp = mean(speciestotal)) %>%
        right_join(temp, by = "stateroute")
      
      # Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) 
      mergespp$comp_scaled = mergespp$meancomp/(mergespp$FocalAbundance + mergespp$meancomp) 
      mergespp$comp_scaled[is.na(mergespp$comp_scaled)] = 0
      
      if(length(unique(mergespp$comp_scaled[!is.na(mergespp$comp_scaled)])) > 2){ 
        lms <- lm(mergespp$occ_logit ~  mergespp$comp_scaled) 
        lms_est = summary(lms)$coef[2,"Estimate"]
        lms_p = summary(lms)$coef[2,"Pr(>|t|)"]
        lms_r = summary(lms)$r.squared
        
        noncomps_output = rbind(noncomps_output, c(FocalAOU, co,lms_est, lms_p, lms_r))
      }
    }
  }
}         

noncomps_output = data.frame(noncomps_output)
names(noncomps_output) = c("FocalAOU", "CompetitorAOU", "Estimate","P", "R2")
# have to remove pairing of American REdstart/least flycatcher, non-familial pairing based on lit
noncomps_output = noncomps_output[!(noncomps_output$FocalAOU == 6870 & noncomps_output$CompetitorAOU == 4670),]

# write.csv(noncomps_output, "data/noncomps_output.csv", row.names = FALSE)
noncomps_output = read.csv("data/noncomps_output.csv", header = TRUE)
noncomps_output_bocc = left_join(noncomps_output, beta_occ[,c("FocalAOU", "Competition_R2", "Competition_Est", "Competition_P")], by = "FocalAOU")
nonps = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(R2 >= Competition_R2)
names(nonps) = c("FocalAOU", "main_g_non")

none = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(Estimate >= Competition_Est)
names(none) = c("FocalAOU", "main_g_non_e")

numcomps = na.omit(noncomps_output) %>% 
  count(FocalAOU)
names(numcomps) = c("FocalAOU", "Comp_count")

vec = c()
for(i in unique(noncomps_output_bocc$FocalAOU)){
  temp = subset(noncomps_output_bocc, FocalAOU == i)
  vec = append(vec, median(temp$R2))
} 
 
vec_meds = data.frame(FocalAOU = unique(noncomps_output_bocc$FocalAOU), Median = vec)  
noncomps_output_ttest = left_join(noncomps_output_bocc, vec_meds, by = "FocalAOU")
t.test(noncomps_output_ttest$Competition_R2, noncomps_output_ttest$Median, paired = TRUE, alternative= "two.sided")


noncompsdist  = merge(nonps, numcomps, by = ("FocalAOU"))
noncompsdist  = merge(none, numcomps, by = ("FocalAOU"))
noncompsdist$nullp = (noncompsdist$main_g_non + 1)/(noncompsdist$Comp_count + 1)
noncompsdist$nulle = (noncompsdist$main_g_non + 1)/(noncompsdist$Comp_count + 1)
nullpsub = filter(noncompsdist, nullp < 0.05) %>% 
  left_join(., envoutput1, by = "FocalAOU")

noncompsdist_trait = merge(noncompsdist, envoutput2[,c("FocalAOU", "migclass", "Trophic.Group")], by = "FocalAOU")

hist(noncompsdist$nullp,xlab = "", main = "Distribution of P-values of non-competitors")
abline(v=mean(noncompsdist$nullp), col = "blue", lwd = 2)

hist(noncomps_output$R2, main = "Distribution of R-squared of non-competitors", xlab = expression('R'^2))
abline(v=mean(beta_occ$Competition_R2), col = "blue", lwd = 2)
hist(na.omit(noncomps_output$Estimate), main = "Distribution of Estimates of non-competitors", xlab = 'Estimate', xlim = c(-40, 40))
abline(v=mean(na.omit(beta_occ$Competition_Est)), col = "blue", lwd = 2)


noncomps_output_bocc$Null = "Null"
noncomps_output_bocc$Comp = "Comp"
noncompsdist$Null = "Null"

R = ggplot(noncomps_output_bocc) +
  # stat_density(aes(Competition_R2, fill=factor(Comp, levels = c("Comp"))), alpha = 0.9) +
  stat_density(aes(R2, fill=factor(Null, levels = c("Null"))), alpha = 0.9) + 
  xlab(expression("Competitor R"^"2")) + ylab("Density") + theme_classic() +
  scale_fill_manual(values=c("purple4")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))
#ggsave("C:/Git/Biotic-Interactions/Figures/null_density_plot_R2.pdf", height = 7, width = 12)

E = ggplot(noncomps_output_bocc) +
  # stat_density(aes(Competition_Est, fill=factor(Comp, levels = c("Comp"))), alpha = 0.9) +
  stat_density(aes(Estimate, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  xlab("Competitor Estimate") + ylab("Density") + theme_classic() + 
  scale_fill_manual(values=c("purple4")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))
#ggsave("C:/Git/Biotic-Interactions/Figures/null_density_plot_Est.pdf", height = 7, width = 12)

P = ggplot(noncomps_output_bocc) +
  stat_density(aes(P, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  xlab("Competitor P-value") + ylab("Density") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("purple4"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))
#ggsave("C:/Git/Biotic-Interactions/Figures/null_density_plot_p.pdf", height = 7, width = 12)

plot_grid(P+ theme(legend.position="none"),
          R + theme(legend.position="none"),
          E + theme(legend.position="none"),
          align = 'h',
          labels = c("A","B", "C"),
          nrow = 1)
ggsave("C:/Git/Biotic-Interactions/Figures/densityplot_null.pdf", height = 7, width = 12)


#### example non-comp dist and main R2 ######
single_dist = subset(noncomps_output_bocc, FocalAOU == 3190)
n = ggplot(single_dist) +
  geom_histogram(bins = 15, aes(R2, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  geom_vline(xintercept = single_dist$Competition_R2, col = "black", lwd = 1.5, lty = 2) +
  xlab(expression("Variance Explained")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("#c994c7"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))

o = ggplot(single_dist) +
  geom_histogram(bins = 15, aes(Estimate, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  geom_vline(xintercept = single_dist$Competition_Est, col = "black", lwd = 1.5, lty = 2) +
  xlab(expression("Competitor Estimate")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("#c994c7"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))

p = ggplot(noncomps_output_bocc) +
  geom_histogram(bins = 20, aes(P, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  xlab(expression("Variance Explained")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("purple4"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))

noncomps_eplot = filter(noncomps_output_bocc, Estimate < 50)
q = ggplot(noncomps_eplot) +
  geom_histogram(bins = 25, aes(Estimate, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  xlab(expression("Competitor Estimate")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("purple4"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24))

plot_grid(n + theme(legend.position="none"),
          o + theme(legend.position="none"),
          p + theme(legend.position="none"),
          q + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C", "D"),
          label_size = 20,
          nrow = 2) 
ggsave("Figures/Figure6_null.pdf", height = 8, width = 12)

#### 1:1 plots ####
noncomps_points = noncomps_output_bocc %>% 
  select(FocalAOU, R2, Competition_R2) %>%
  group_by(FocalAOU) %>%
  summarise(nonr2 = mean(R2)) 

noncomps_1.1 = left_join(noncomps_points, noncomps_output_bocc[, c("FocalAOU", "Competition_R2")], by = "FocalAOU")

ggplot(aes(Competition_R2, nonr2), data = noncomps_1.1) + geom_point() + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25) + ylab("Mean Non-competitor R2") + xlab("Main Competitor R2")


noncomps_points_med = noncomps_output_bocc %>% 
  select(FocalAOU, R2, Competition_R2) %>%
  group_by(FocalAOU) %>%
  summarise(nonr2 = median(R2)) 

noncomps_med = left_join(noncomps_points_med, noncomps_output_bocc[, c("FocalAOU", "Competition_R2")], by = "FocalAOU")

ggplot(aes(Competition_R2, nonr2), data = noncomps_med) + geom_point() + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25) + ylab("Median Non-competitor R2") + xlab("Main Competitor R2")

#### non comp plots ####
# noncomps_output = merge(noncomps_output, nsw[,c("focalAOU", "Family")], by.x = "FocalAOU", by.y = "focalAOU")
noncomps_output = noncomps_output %>% arrange(FocalAOU)
pdf('Figures/noncomp_est.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
for (sp in unique(noncomps_output$FocalAOU)){
  temp = subset(noncomps_output, FocalAOU == sp) 
  comp_est = beta_occ$Competition_Est[beta_occ$FocalAOU == sp]
  hist(temp$Estimate, xlab = 'Estimate', 
       main = sp, xlim = c(min(temp$Estimate, comp_est), max(temp$Estimate, comp_est)))
  abline(v = comp_est, col = "blue", lwd = 2)
  
}
dev.off()

pdf('Figures/noncomp_r2.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
for (sp in unique(noncomps_output$FocalAOU)){
  temp = subset(noncomps_output, FocalAOU == sp) 
  comp_r2 = beta_occ$Competition_R2[beta_occ$FocalAOU == sp]
  hist(temp$R2, xlab = expression('R'^2), main = sp,
       xlim = c(0, max(temp$R2, comp_r2)))
  abline(v = comp_r2, col = "blue", lwd = 2)
}
dev.off()

noncomps_output$type = "Species"
ggplot(noncomps_output, aes(type, R2)) + geom_violin(linetype = "blank", aes(fill = factor(noncomps_output$type))) + xlab("Total Variance") + ylab("R2")+ theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_blank()) 

ggsave("C:/Git/Biotic-Interactions/Figures/violin_noncomps.png")

noncompsdist$type = "Species"
ggplot(noncompsdist, aes(type, nullp)) + geom_violin(linetype = "blank", aes(fill = factor(noncompsdist$type))) + xlab("Total Variance") + ylab("P-val")+ theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_blank()) 




#### ---- Plotting GLMs ---- ####
# Making pdf of ranges for each focal spp
pdf('precip_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occuenv[occuenv$Species == sp,]
  competition <- lm(psub$occ_logit ~  psub$comp_scaled) 
  # z scores separated out for env effects (as opposed to multivariate variable)
  env_z = lm(psub$occ_logit ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = psub)
  # z scores separated out for env effects
  both_z = lm(psub$occ_logit ~  comp_scaled + abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = psub)
  
  tes = ggplot(data = psub, aes(x =  psub$comp_scaled, y = psub$occ_logit)) +stat_smooth(data=env_z, lwd = 1.5,se = FALSE) +xlab(psub$Species)+theme_bw()
  plot(tes)
}
dev.off()
# Making pdf of ranges for each focal spp
pdf('Figures/allspp_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  temp = occumatrix[occumatrix$Species == sp,]
  lm = lm(temp$occ_logit ~  abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
  
  tes = ggplot(data = temp, aes(x = zTemp, y = FocalOcc)) +stat_smooth(method = "lm", lwd = 1.5,se = FALSE) +xlab(temp$Species) +ylim(0,1) +theme_bw()
  plot(tes)
}
dev.off()
# Making pdf of ranges for each focal spp
pdf('rangmaps.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occumatrix[occumatrix$Species == sp,]
  glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ comp_scaled + 
                              abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI) + (1|stateroute:Species), family = binomial(link = logit), data = psub)
  
  tes = ggplot(data = psub, aes(x = comp_scaled, y = FocalOcc)) +stat_smooth(data=glm_occ_rand_site, lwd = 1.5,se = FALSE) +theme_bw()
  plot(tes)
}
dev.off()

pTemp = predict(glm_occ_rand_site, newdata=with(occumatrix,data.frame(zTemp=0,comp_scaled,zPrecip,zElev,zNDVI,stateroute,Species, FocalOcc)), allow.new.levels = TRUE) #predict values assuming zTemp=0

inverselogit <- function(p) {exp(p)/(1+exp(p))} 
newintercept <- function(p) {mean(exp(p)/(1+exp(p)))} 

# this relationship should be negative
temp = ggplot(data = occumatrix, aes(x = abs(zTemp), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zTemp)), yend = 0), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02)) + theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logittemp.png")

ndvi = ggplot(data = occumatrix, aes(x = abs(zNDVI), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zNDVI)), yend = 0), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitndvi.png")

elev = ggplot(data = occumatrix, aes(x = abs(zElev), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zElev)), yend = 0), col = "dark green", lwd=2)  + 
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitelev.png")

precip = ggplot(data = occumatrix, aes(x = abs(zPrecip), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zPrecip)), yend = 0), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitprecip.png")

z <- plot_grid(precip+ theme(legend.position="none"),
               ndvi + theme(legend.position="none"),
               align = 'h',
               labels = c("A","B"),
               nrow = 1)
p2 = plot_grid(elev + theme(legend.position="none"),
               temp + theme(legend.position="none"), 
               labels = c("C","D"),
               align = 'h', 
               rel_widths = c(1, 1.3))

plot_grid(z, p2, ncol = 1, rel_heights = c(1, 1))
ggsave("C:/Git/Biotic-Interactions/Figures/cowplotabiotic.pdf")

ggplot(data = occumatrix, aes(x = comp_scaled, y = FocalOcc)) + 
  stat_function(fun=inverselogit, color = "blue") + 
  geom_point(colour="black", shape=18, alpha = 0.02,position=position_jitter(width=0,height=.02))+ theme_classic()

tableS1 = read.csv("//bioark/HurlbertLab/Snell/2018 BI MS/Tables/Table S1.csv", header = TRUE)
tableS1 = tableS1 %>%
  group_by(focalAOU) %>%
  arrange(desc(mainCompetitor), .by_group = TRUE) %>%
  group_by(CompAOU)
tableS1 = left_join(tableS1, subsetocc[,c("AOU", "migclass", "Foraging", "Trophic.Group")], by = c("focalAOU" ="AOU"))
# write.csv(tableS1, "data/Table S1.csv", row.names = FALSE)
