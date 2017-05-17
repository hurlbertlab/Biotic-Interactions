library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)

# read in files created in data cleaning script
tax_code = read.csv("data/Tax_AOU_Alpha.csv", header = TRUE) # Hurlbert Lab
temp_occ = read.csv("data/bbs_sub1.csv", header=TRUE) # BBS occ script
centroid=read.csv("data/centroid.csv", header=TRUE) # GIS script
occuenv= read.csv("data/all_expected_pres.csv", header = TRUE) # Data cleaning script
Hurlbert_o = read.csv('data/Master_RO_Correlates_20110610.csv', header = T) # Hurlbert Lab
subsetocc = Hurlbert_o[Hurlbert_o$X10yr.Prop > .3 & Hurlbert_o$X10yr.Prop < .7,]

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

beta_lm = matrix(NA, nrow = 102, ncol = 10)
beta_abun = matrix(NA, nrow = 102,ncol = 10)

# create beta output data frame

beta_lm[sp,1] = sp
beta_lm[sp,2] = summary(competition)$coef[1,"Estimate"]
beta_lm[sp,3] = summary(competition)$coef[1,"Pr(>|t|)"]
beta_lm[sp,4] = summary(competition)$r.squared #using multiple rsquared
beta_lm[sp,5] = summary(env_z)$coef[2,"Estimate"]
beta_lm[sp,6] = summary(env_z)$coef[2,"Pr(>|t|)"]
beta_lm[sp,7] = summary(env_z)$r.squared 
beta_lm[sp,8] = summary(both_z)$coef[2,"Estimate"]
beta_lm[sp,9] = summary(both_z)$coef[2,"Pr(>|t|)"]
beta_lm[sp,10] = summary(both_z)$r.squared 

beta_abun[sp,1] = subfocalspecies[sp]
beta_abun[sp,2] = summary(competition_abun)$coef[2,"Estimate"]
beta_abun[sp,3] = summary(competition_abun)$coef[2,"Pr(>|t|)"]
beta_abun[sp,4] = summary(competition_abun)$r.squared #using multiple rsquared
beta_abun[sp,5] = summary(env_abun)$coef[2,"Estimate"]
beta_abun[sp,6] = summary(env_abun)$coef[2,"Pr(>|t|)"]
beta_abun[sp,7] = summary(env_abun)$r.squared 
beta_abun[sp,8] = summary(both_abun)$coef[2,"Estimate"]
beta_abun[sp,9] = summary(both_abun)$coef[2,"Pr(>|t|)"]
beta_abun[sp,10] = summary(both_abun)$r.squared

occuenv = na.omit(occuenv)
subfocalspecies = unique(occuenv$Species)

for (sp in 1:length(subfocalspecies)){
  print(sp)
  temp = subset(occuenv,occuenv$Species == subfocalspecies[sp])
  
  competition <- lm(temp$occ_logit ~  temp$comp_scaled) 
  # z scores separated out for env effects (as opposed to multivariate variable)
  env_z = lm(temp$occ_logit ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  # z scores separated out for env effects
  both_z = lm(temp$occ_logit ~  temp$comp_scaled + abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
  
  # abundance, not temp occ - same results?
  competition_abun <- lm(temp$FocalAbundance ~  temp$comp_scaled) 
  # z scores separated out for env effects - abundance
  env_abun = lm(temp$FocalAbundance ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  # z scores separated out for env effects - abundance
  both_abun = lm(temp$FocalAbundance ~  comp_scaled + abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  
  #variance_partitioning 
  ENV = summary(both_z)$r.squared - summary(competition)$r.squared
  print(ENV) #env only
  COMP = summary(both_z)$r.squared - summary(env_z)$r.squared
  print(COMP) #competition only
  SHARED = summary(competition)$r.squared - COMP
  print(SHARED) #shared variance
  NONE = 1 - summary(both_z)$r.squared
  print(NONE) #neither variance
  sp1 = unique(temp$Species)
  sum = sum(ENV, COMP, SHARED)
  envoutput = rbind(envoutput, c(sp1, ENV, COMP, SHARED, NONE))
  
  #variance_partitioning 
  ENVa = summary(both_abun)$r.squared - summary(competition_abun)$r.squared
  
  COMPa = summary(both_abun)$r.squared - summary(env_abun)$r.squared
  
  SHAREDa = summary(competition_abun)$r.squared - COMP
  
  NONEa = 1 - summary(both_abun)$r.squared
  
  sp1 = unique(temp$Species)
  envoutputa = rbind(envoutputa, c(sp1, ENVa, COMPa, SHAREDa, NONEa))
}         


envoutput = data.frame(envoutput)
envoutputa = data.frame(envoutputa)
names(envoutput) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")
names(envoutputa) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")
# relabel dark-eyed junco
envoutput$FocalAOU[envoutput$FocalAOU == 5660] <- 5677
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
subsetocc$AOU[subsetocc$AOU == 5660] <- 5677
subsetocc$AOU[subsetocc$AOU == 7220] <- 7222

envoutput1 = merge(envoutput, tax_code[,c('AOU_OUT', 'ALPHA.CODE')], by.x = 'FocalAOU', by.y = "AOU_OUT", all.x = TRUE) # losing Alder Flycatcher

envoutput2 = merge(envoutput, subsetocc[,c("AOU", "migclass", "Trophic.Group")], by.x='FocalAOU', by.y='AOU', all.x = TRUE)

envloc = merge(envoutput2, centroid[, c("FocalAOU", "Long", "Lat")], by = 'FocalAOU', all.x = TRUE)

#write.csv(envoutput, "envoutput.csv", row.names = FALSE)
#write.csv(envoutputa, "envoutputa.csv", row.names = FALSE)
beta_lm = data.frame(beta_lm)
names(beta_lm) = c("FocalAOU", "Competition_Est", "Competition_P", "Competition_R2", "EnvZ_Est", "EnvZ_P", "EnvZ_R2", "BothZ_Est", "BothZ_P", "BothZ_R2")
beta_abun = data.frame(beta_abun)
names(beta_abun) = c("FocalAOU", "Competition_Est", "Competition_P", "Competition_R2", "EnvZ_Est", "EnvZ_P", "EnvZ_R2", "BothZ_Est", "BothZ_P", "BothZ_R2")

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
occumatrix$nyears = occumatrix$FocalOcc*15
occumatrix$sp_success = as.factor(occumatrix$nyears * occumatrix$FocalOcc)
occumatrix$sp_fail = as.factor(occumatrix$nyears * (1 - occumatrix$FocalOcc))

#### GLM of all matrices not just subset ####
glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ c_s + 

                            abTemp + abElev + abPrecip + abNDVI + (1|stateroute:Species), family = binomial(link = logit), data = occumatrix)
summary(glm_occ_rand_site) 

### FIX, Neg Binom
# glm_abun_rand_site = glmer.nb(FocalAbundance ~ c_s + 
                            #    abTemp + abElev + abPrecip + abEVI + (1|stateroute:Species), link=log, data = occumatrix)
# summary(glm_abundance_rand_site) 

#### PLOTTING MODELS ####
ggplot(data = occumatrix, aes(x = comp_scaled, y = FocalOcc)) +stat_smooth(data=glm_occ_rand_site, lwd = 1.5) +xlab("Scaled Competitor Abundance")+ylab("Focal Occupancy") +theme_bw() +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90), axis.text=element_text(size=12)) + theme(plot.margin = unit(c(.5,6,.5,.5),"lines")) 
ggsave("C:/Git/Biotic-Interactions/Figures/glmoutput.png")

####### GLM FIT PLOTS #################################################################################################
pTemp = predict(glm_occ_rand_site, newdata=with(occumatrix,data.frame(zTemp=0,comp_scaled,zPrecip,zElev,zEVI,stateroute,Species, FocalOcc)), allow.new.levels = TRUE) #predict values assuming zTemp=0

inverselogit <- function(p) {exp(p)/(1+exp(p))} 
newintercept <- function(p) {mean(exp(p)/(1+exp(p)))} 

# this relationship should be negative
ggplot(data = occumatrix, aes(x = abs(zTemp), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 0.892997, xend = abs(max(occumatrix$zTemp)), yend = 0.892997 +(-0.044095*max(abs(occumatrix$zTemp)))), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02)) + theme_classic()
#geom_abline(intercept=0.97569, slope= -0.05176, lwd=1, col="blue")
ggsave("C:/Git/Biotic-Interactions/Figures/logittemp.png")

ggplot(data = occumatrix, aes(x = abs(zEVI), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 0.892997, xend = abs(max(occumatrix$zEVI)), yend = 0.892997 +(-0.094849*max(abs(occumatrix$zEVI)))), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
# geom_abline(intercept=0.97569, slope= -0.05117, lwd=1, col="blue")
ggsave("C:/Git/Biotic-Interactions/Figures/logitevi.png")

ggplot(data = occumatrix, aes(x = abs(zElev), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 0.892997, xend = abs(max(occumatrix$zElev)), yend = 0.892997 +(-0.010785*max(abs(occumatrix$zElev)))), col = "dark green", lwd=2)  + 
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitelev.png")

ggplot(data = occumatrix, aes(x = abs(zPrecip), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 0.892997, xend = abs(max(occumatrix$zPrecip)), yend = 0.892997 +(0.001575*max(abs(occumatrix$zPrecip)))), col = "dark green", lwd=2) +
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitprecip.png")

ggplot(data = occumatrix, aes(x = comp_scaled, y = FocalOcc)) + 
  stat_function(fun=inverselogit, color = "blue") + 
  geom_point(colour="black", shape=18, alpha = 0.02,position=position_jitter(width=0,height=.02))+ theme_classic()

#### ---- Plotting GLMs ---- ####
# Making pdf of ranges for each focal spp
pdf('precip_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occumatrix[occumatrix$Species == sp,]
  glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ comp_scaled + 
                              abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zEVI) + (1|stateroute:Species), family = binomial(link = logit), data = psub)
  
  tes = ggplot(data = psub, aes(x = zPrecip, y = FocalOcc)) +stat_smooth(data=glm_occ_rand_site, lwd = 1.5,se = FALSE) +xlab(psub$Species)+theme_bw()
  plot(tes)
}
dev.off()
# Making pdf of ranges for each focal spp
pdf('Temp_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occumatrix[occumatrix$Species == sp,]
  glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ comp_scaled + 
                              abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zEVI) + (1|stateroute:Species), family = binomial(link = logit), data = psub)
  
  tes = ggplot(data = psub, aes(x = zTemp, y = FocalOcc)) +stat_smooth(data=glm_occ_rand_site, lwd = 1.5,se = FALSE) +xlab(psub$Species) +theme_bw()
  plot(tes)
}
dev.off()
# Making pdf of ranges for each focal spp
pdf('GLM_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occumatrix[occumatrix$Species == sp,]
  glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ comp_scaled + 
                              abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zEVI) + (1|stateroute:Species), family = binomial(link = logit), data = psub)
  
  tes = ggplot(data = psub, aes(x = comp_scaled, y = FocalOcc)) +stat_smooth(data=glm_occ_rand_site, lwd = 1.5,se = FALSE) +theme_bw()
  plot(tes)
}
dev.off()

##### Variance Partitioning Plot #####
envloc$EW <- 0

envloc$EW[envloc$Long > -98.583333] <- 1 ## mid point of USA
# 1 = East

envloc1 = merge(envloc, occumatrix[,c("Species", "Family")], by.x = "FocalAOU", by.y="Species", all.y = FALSE)
# cutting down to unique subset
envloc1 = unique(envloc1)
envloc1 = merge(envloc1, tax_code[,c("AOU_OUT", "ALPHA.CODE")], by.x = "FocalAOU", by.y = "AOU_OUT", all.x=TRUE)
envloc1$ALPHA.CODE = as.character(envloc1$ALPHA.CODE)
envloc1$ALPHA.CODE[envloc1$FocalAOU == 2920] <- 'MOUQ' #Mountain Quail
envloc1$ALPHA.CODE[envloc1$FocalAOU == 6720] <- 'PAWA' #Palm Warbler
envloc1$ALPHA.CODE = as.factor(envloc1$ALPHA.CODE)

nrank = envloc1 %>% 
  mutate(rank = row_number(-ENV))# change here for comp
envflip = tidyr::gather(nrank, "Type", "value", 2:5)
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)])
envflip = plyr::arrange(envflip,(envflip$rank),envflip$FocalAOU)

envflip = merge(envflip, envloc[,c("FocalAOU", "EW")], by = "FocalAOU")

envrank = envflip %>% 
  dplyr::group_by(Type == 'ENV') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "ENV")

### CREATE LABEL DF FAMilY ########
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

envrank$EW.x[envrank$EW.x == 1] <- "E"
envrank$EW.x[envrank$EW.x == 0] <- "W" 
###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("ENV","COMP","SHARED","NONE"),ordered = TRUE)
# Plot with ENV ranked in decreasing order
t = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("ENV","COMP","SHARED","NONE")))) + 
  geom_bar(stat = "identity")  + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("#2ca25f","#dd1c77","#43a2ca","white"), labels=c("Environment", "Competition","Shared Variance", "")) +theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20, angle=90),legend.title=element_text(size=12), legend.text=element_text(size=20), legend.position="top", legend.justification=c(0, 1), legend.key.width=unit(1, "lines")) + guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))

tt = t + annotate("text", x = 1:102, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=0.5, color = "black") + annotate("text", x = 1:102, y = -.08, label = envrank$mig_abbrev, size=6,vjust=0.5, color = envrank$mig_abbrevf, fontface =2) + annotate("text", x = 1:102, y = -.1, label = envrank$trophlabel, size=6,vjust=0.5, color = envrank$trophlabelf, fontface =2) + annotate("text", x = 1:102, y = -.12, label = envrank$EW.x, angle=90,size=6,vjust=0.5, color = "black", fontface =2)+ annotate("text", x = 1:102, y = -.06, label = envrank$Fam_abbrev, size=6,vjust=0.5, color = "black", fontface =2) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks=scales::pretty_breaks(9))

plot(tt)

ggsave("C:/Git/Biotic-Interactions/Figures/barplot.pdf", height = 30, width = 38)

##################### TRAITS Model ####################################
logit = function(x) log(x/(1-x))

env_lm = subset(envflip, Type == 'ENV')

env_traits = lm(logit(value) ~ Trophic.Group + migclass + EW.x, data = env_lm)
summary(env_traits) 

comp_lm = subset(envflip, Type == 'COMP')

comp_traits = lm(logit(value) ~ Trophic.Group + migclass + EW.x, data = comp_lm[comp_lm$value > 0,])
summary(comp_traits) 

env_sum = subset(envflip, Type != 'NONE')
total = env_sum %>% 
  dplyr::group_by(FocalAOU) %>%
  summarise(sum(value))

env_sum$edgeval = (env_sum$value * (1 - 2*edge_adjust)) + edge_adjust
total_traits = lm(logit(edgeval) ~ Trophic.Group + migclass + EW.x, data = env_sum)
summary(total_traits)

# R2 plot - lm in ggplot
#envoutputa = read.csv("envoutputa.csv", header = TRUE)
#names(envoutputa) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")
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

ggplot(R2plot2, aes(x = COMP.x, y = COMP.y)) +theme_bw()+ theme(axis.title.x=element_text(size=35),axis.title.y=element_text(size=35, angle=90)) + xlab("Occupancy R2") + ylab("Abundance R2") + geom_point(col = "#dd1c77", cex =4, shape=24) + geom_point(data = R2plot2, aes(x = ENV.x, y = ENV.y), shape = 16, col = "#2ca25f", cex =4, stroke = 1) + geom_point(data = R2plot2, aes(Total.x,Total.y), shape = 3, col = "dark gray", cex =5, stroke = 1) +geom_abline(intercept = 0, slope = 1, col = "red", lwd = 1.25)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=2))
ggsave("C:/Git/Biotic-Interactions/Figures/occvabun.png")

# need to change the slopes
ggplot(R2plot2, aes(x = COMP.x, y = COMP.y)) +theme_bw()+ theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16, angle=90)) + xlab("Occupancy R2") + ylab("Abundance R2") + geom_point(col = "#dd1c77", cex =4, shape=24) + geom_point(data = R2plot2, aes(x = ENV.x, y = ENV.y), shape = 16, col = "#2ca25f", cex =4, stroke = 1) + geom_point(data = R2plot2, aes(Total.x,Total.y), shape = 3, col = "#43a2ca", cex =5, stroke = 1) +geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash") + geom_abline(intercept= 0.03536,slope=0.43437, col = "#2ca25f", lwd = 1,linetype="dotdash")+ geom_smooth(method="lm", se= F, size = 1, aes(linetype = "dotdash", group = ENV.y))+geom_abline(intercept=0.1424,slope=0.4193, col = "#43a2ca", lwd = 1,linetype="dotdash")



# R2 plot - glm violin plots
ggplot(R2plot2, aes(x = FocalAOU, y = Total.x)) + geom_violin(lwd = 2, fill = "grey", color = "grey") + xlab("Total Variance") + ylab("R2")+ theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))+scale_y_continuous(limits = c(0, 1)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25))

ggplot(R2plot2, aes(x = FocalAOU, y = COMP.x)) + geom_violin(lwd = 2, fill = "#dd1c77", color = "#dd1c77") + xlab("Competition") + ylab("R2")+ theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30)) +scale_y_continuous(limits = c(0, 1))+ theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25))

ggplot(R2plot2, aes(x = FocalAOU, y = ENV.x)) + geom_violin(lwd = 2, fill = "#2ca25f", color = "#2ca25f") + xlab("Environment") + ylab("R2")+ theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),legend.title=element_text(size=12), legend.text=element_text(size=12)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25)) +scale_y_continuous(limits = c(0, 1))

#Coyle fig 1: Z:\Coyle\Projects\BBS Core\Final Analysis

#### basement ### poster
forplots = merge(Hurlbert_o[, c("AOU","Trophic.Group","Foraging","migclass")], AOUsub4[, c("AOU","Family", "CommonName")], by = "AOU")
# midpoint long of US is -98.5795, so 1 indicates east of that line, 0 = west

# Dark-eyed Junko and Winter Wren need AOUs changed
tax_code$AOU_OUT[tax_code$AOU_OUT == 5677] = 5660
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] = 7222
forplots$AOU[forplots$AOU == 7220] = 7222

envloc = merge(envoutput, centroid[, c("FocalAOU", "Long", "Lat")], by = 'FocalAOU', all = TRUE)
envloc$EW <- 0
envloc$EW[envloc$Long > -98.583333] <- 1 ## from https://tools.wmflabs.org/geohack/geohack.php?pagename=Geographic_center_of_the_contiguous_United_States&params=39_50_N_98_35_W_region:US-KS_type:landmark&title=Geographic+Center+of+the+Contiguous+United+States
# 1 = East
#### ---- Plotting LMs ---- ####
# Making pdf of ranges for each focal spp
pdf('Figures/Lin_Reg.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))
# Plotting basic lms to understand relationships
for(sp in subfocalspecies){ 
  print(sp)
  psub = occuenv[occuenv$Species == sp,]
  #psub = filter(psub, occ_logit < 9) # eliminating 100% occ?
  title = unique(psub$FocalSciName)
  #ggplot(psub, aes(x = psub$occ_logit, y = psub$MainCompSum)) + geom_point(data=psub, pch = 16)+geom_smooth(method = "lm", col = "red")+ theme_classic()+ xlab("Focal Occupancy")+ylab("Competitor Abundance")+ggtitle(title)
  
  #+ ggtitle(title[1])
  competition <- lm(psub$occ_logit ~  psub$comp_scaled) 
  # z scores separated out for env effects (as opposed to multivariate variable)
  env_z = lm(occ_logit ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zEVI), data = psub)
  # z scores separated out for env effects
  both_z = lm(psub$occ_logit ~  psub$comp_scaled + abs(psub$zTemp)+abs(psub$zElev)+abs(psub$zPrecip)+abs(psub$zEVI), data = psub)
  
  plot(psub$comp_scaled, psub$occ_logit, pch = 20, xlab = "Main Competitor Abundance", ylab = "Focal Occupancy (logit link)", main = psub$FocalSciName[1], sub = "Competition", abline(competition, col = "red"))
  plot(psub$comp_scaled,psub$occ_logit,  pch = 20, xlab = "Main Competitor Abundance", ylab = "Focal Occupancy (logit link)", main = psub$FocalSciName[1], sub = "Environment", abline(env_z, col = "red"))
  plot(psub$comp_scaled, psub$occ_logit,  pch = 20, xlab = "Main Competitor Abundance", ylab = "Focal Occupancy (logit link)", main = psub$FocalSciName[1], sub = "Both", abline(both_z, col = "red"))
}
dev.off()

# Plotting basic lm hists to understand relationships btwn occ and abun
hist(beta_lm$Competition_R2, 10, main = "R Squared Distribution for Competition", xlab = "Competition R Squared")
hist(beta_lm$EnvZ_R2, 10, main = "R Squared Distribution for Env", xlab = "Env R Squared")
hist(beta_lm$BothZ_R2, 10, main = "R Squared Distribution for Both", xlab = "Both R Squared")

hist(beta_lm$Competition_Est, 10, main = "Slope Distribution for Competition", xlab = "Competition Slope")
abline(v = mean(beta_lm$Competition_Est), col = "red", lwd = 3)
hist(beta_lm$EnvZ_Est, 10, main = "Slope Distribution for Environment", xlab = "Environment Slope")
abline(v = mean(beta_lm$EnvZ_Est), col = "red", lwd = 3)
hist(beta_lm$BothZ_Est, 10, main = "Slope Distribution for Both", xlab = "Both Slope")
abline(v = mean(beta_lm$BothZ_Est), col = "red", lwd = 3)

hist(beta_abun$Competition_R2, 10, main = "R Squared Distribution for Competition", xlab = "Competition R Squared")
hist(beta_abun$EnvZ_R2, 10, main = "R Squared Distribution for Env", xlab = "Env R Squared")
hist(beta_abun$BothZ_R2, 10, main = "R Squared Distribution for Both", xlab = "Both R Squared")

hist(beta_abun$Competition_Est, 10, main = "Slope Distribution for Competition", xlab = "Competition Slope")
abline(v = mean(beta_abun$Competition_Est), col = "red", lwd = 3)
hist(beta_abun$EnvZ_Est, 10, main = "Slope Distribution for Environment", xlab = "Environment Slope")
abline(v = mean(beta_abun$EnvZ_Est), col = "red", lwd = 3)
hist(beta_abun$BothZ_Est, 10, main = "Slope Distribution for Both", xlab = "Both Slope")
abline(v = mean(beta_abun$BothZ_Est), col = "red", lwd = 3)

beta_lm$sumR2 = beta_lm$BothZ_R2+beta_lm$Competition_R2+beta_lm$EnvZ_R2  # tbd if we'll keep this