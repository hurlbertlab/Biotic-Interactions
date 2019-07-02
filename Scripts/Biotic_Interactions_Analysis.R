library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggExtra)
library(rstanarm)
library(boot)

# We can flip back and forth between comp_scaled and all_comp_scaled with the find and replace buttons. There should be 26 occurrences

# read in files created in data cleaning script
temp_occ <- read.csv("data/bbs_sub1.csv", header=TRUE) # BBS occ script
centroid<-read.csv("data/centroid.csv", header=TRUE) # GIS script
occuenv <- read.csv("data/all_expected_pres.csv", header = TRUE) # Data cleaning script
subsetocc <- read.csv('data/subsetocc.csv', header = T) # Hurlbert Lab
tax_code <- read.csv("data/Tax_AOU_Alpha.csv", header = TRUE) # Hurlbert Lab
bbs_abun <- read.csv("data/bbs_abun.csv", header =TRUE)
nsw <- read.csv("data/new_spec_weights.csv", header = TRUE)
shapefile_areas <- read.csv("data/shapefile_areas.csv", header =TRUE)
sppGT50rtes <- read.csv("data/sppGT50rtes.csv", header = TRUE)

#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
subsetocc$AOU[subsetocc$AOU == 4810] = 4812
temp_occ$Aou[temp_occ$Aou == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4123] = 4120

# have to omit focal/competitors with no overlap
shapefile_areas = na.omit(shapefile_areas)
# sum all comps for each focal, divide by focal range
shapefile_overlap = shapefile_areas %>%
  group_by(focalAOU) %>%
  summarise(all_overlap = sum(area_overlap))
shapefile_overlap = left_join(shapefile_overlap, shapefile_areas[,c("Focal", "focalAOU", "FocalArea")], by = "focalAOU")
shapefile_overlap$prop_overlap = shapefile_overlap$all_overlap/shapefile_overlap$FocalArea

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
  
  competition <- lm(temp$occ_logit ~  temp$all_comp_scaled)  # changes between main and all comps
  # z scores separated out for env effects 
  env_z = lm(temp$occ_logit ~ abs(zTemp) + abs(zElev) + abs(zPrecip) + abs(zNDVI), data = temp)
  # z scores separated out for env effects
  both_z = lm(temp$occ_logit ~  temp$all_comp_scaled + abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
  
  # abundance, not temp occ - same results?
  competition_abun <- lm(temp$FocalAbundance ~  temp$all_comp_scaled) 
  # z scores separated out for env effects - abundance
  env_abun = lm(temp$FocalAbundance ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  # z scores separated out for env effects - abundance
  both_abun = lm(temp$FocalAbundance ~  all_comp_scaled + abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
  
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
  SHAREDa = summary(competition_abun)$r.squared - COMPa
  NONEa = 1 - summary(both_abun)$r.squared
  
  sp1 = unique(temp$FocalAOU)
  
  
  envoutputa = rbind(envoutputa, c(sp1, ENVa, COMPa, SHAREDa, NONEa))
  
  if(length(unique(temp$all_comp_scaled[!is.na(temp$all_comp_scaled)])) > 2){
  # saving model output into separate data frames
  occ_comp_est = summary(competition)$coef[2,"Estimate"]
  occ_comp_p = summary(competition)$coef[2,"Pr(>|t|)"]
  occ_comp_r = summary(competition)$r.squared
  occ_env_est = mean(summary(env_z)$coef[,"Estimate"])
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

envoutput2 = merge(envoutput, subsetocc[,c("AOU", "CommonName", "migclass", "Trophic.Group")], by.x='FocalAOU', by.y='AOU', all.x = TRUE)

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

beta_plot <- left_join(beta_occ, occuenv, by = "FocalAOU")
ggplot(beta_plot, aes(x = NDVI_est, y = occ)) + geom_point()

#### ---- GLM fitting  ---- ####
# add on success and failure columns by creating # of sites where birds were found
# and # of sites birds were not found from original bbs data
# occumatrix = merge(temp_occ, occuenv, by.x=c("Aou", "stateroute"),by.y=c("Species", "stateroute"))
occumatrix=occuenv
occumatrix$c_s = scale(occumatrix$all_comp_scaled, scale = T, center = T)
occumatrix$abTemp=abs(occumatrix$zTemp)
occumatrix$abElev=abs(occumatrix$zElev)
occumatrix$abPrecip=abs(occumatrix$zPrecip)
occumatrix$abNDVI=abs(occumatrix$zNDVI)

# using equation species sum*Focal occ to get success and failure for binomial anlaysis
occumatrix$sp_success = as.integer(15 * occumatrix$FocalOcc)
occumatrix$sp_fail = as.integer(15 * (1 - occumatrix$FocalOcc))

#### Bayesian of all matrices not just subset ####
library(brms)
library(rstudioapi)
# mmslope <- brm(sp_success | trials(sp_success+sp_fail) ~ c_s +  abTemp + abElev + abPrecip + abNDVI + (c_s + abTemp + abElev + abPrecip + abNDVI|FocalAOU), family = binomial(link = logit), data = occumatrix , cores = 2, chains=4, iter=500,warmup=200,control = list(max_treedepth = 15),set_prior("lkj(1)", class = "cor"))

# save(mmslope, filename ="mmslope.rda")

mod2 <- readRDS("Z:/Snell/2018 BI MS/mmslope_all_5000.rds")
TableS13 <- readRDS("Z:/Snell/2018 BI MS/mmslope_main_5000.rds")
samples1 <- posterior_summary(mod2)
sample_13 = posterior_summary(TableS13)
samples2 <- posterior_samples(mod2)
# write.csv(sample_13, "data/TableS13.csv", row.names = TRUE)

Table8 = mod2[["fit"]]
mod2sum = summary(mod2)


mmslope5000 = read.csv("data/TableS12.csv", header = TRUE)
mmslope5000$X = gsub("r_", "", mmslope5000$X) 
mmslope5000$X = gsub("(Intercept)", "", mmslope5000$X) 
mmslope5000$X = gsub("FocalAOU", "", mmslope5000$X) 
mmslope5000$X = gsub("[[]", "", mmslope5000$X) 
mmslope5000$X = gsub("[]]", "", mmslope5000$X)
mmslope5000 = mmslope5000 %>% separate(X, c("AOU", "Var"), ",", remove= TRUE)
mmslope5000$AOU = as.factor(mmslope5000$AOU)
mmslope5000 = mmslope5000[,c("AOU", "Var", "Estimate", "Q2.5", "Q97.5")]

# mmslope5000$AOU = as.integer(mmslope5000$AOU)
# mm3 = left_join(mmslope5000, nsw[,c("focalAOU", "Focal", "Family")], by = c("AOU" = "focalAOU"))
mm3 = unique(mm3) %>% group_by(Family)
# write.csv(mm3, "data/tableS12.csv", row.names = FALSE)

occumatrix$cspred = inv.logit(occumatrix$all_comp_scaled * mm_fixed$mean[2] + mm_fixed$mean[1])
occumatrix$temppred = inv.logit(occumatrix$abTemp * mm_fixed$mean[3] + mm_fixed$mean[1])
occumatrix$elevpred = inv.logit(occumatrix$abElev * mm_fixed$mean[4] + mm_fixed$mean[1])
occumatrix$precippred = inv.logit(occumatrix$abPrecip * mm_fixed$mean[5] + mm_fixed$mean[1])
occumatrix$ndvipred = inv.logit(occumatrix$abNDVI * mm_fixed$mean[6] + mm_fixed$mean[1])

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

envrank = envflip %>% 
  dplyr::group_by(Type == 'COMP') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "COMP") # change here for comp

###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE","SHARED", "ENV","COMP"),ordered = TRUE)
envflip$value = abs(envflip$value)
# Plot with ENV ranked in decreasing order - had to flip everything to plot right
t = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE","SHARED", "ENV","COMP")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#2ca25f","#dd1c77"), labels=c("","Shared Variance","Environment", "Competition")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30),legend.title=element_blank(), legend.text=element_text(size=28, hjust = 1, vjust = 0.5), legend.position = c(0.5,.9)) + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))

tt = t + annotate("text", x = 1:175, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=0.5,hjust = 0.8, color = "black") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks = c(0,0.2,0.4,0.6, 0.8))

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
c = ggplot(data=envflip_sub, aes(factor(rank), y=abs(value), fill=factor(Type, levels = c("NONE","SHARED", "ENV","COMP")))) + geom_bar(stat = "identity") + theme_classic() + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#2ca25f","#dd1c77"), labels=c("","Shared","Environment", "Competition")) +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=30, vjust = 3),legend.title=element_blank(), legend.text=element_text(size=34, hjust = 1, vjust = 0.5), legend.position = c(.3,.92)) +theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.key.width=unit(4, "point"), legend.key.height =unit(3, "lines")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size = 30,color= "black")) + scale_y_continuous(breaks = c(0, 0.2,  0.4,  0.6,  0.8),limits = c(0, 0.8),sec.axis = sec_axis(~ . *1)) + annotate("text", x = 1:15, y = -.01, label = envrank$ALPHA.CODE[1:15], angle=90,size=6,vjust=0.5,hjust = 1, color = "black")  # + annotate("text", x = 1:15, y = 0.3, label = envflip_sub3$PRIMARY_COM_NAME.y, angle=90,size=8,vjust=0.5,hjust = 1.4, color = "white", fontface = "bold") + theme(legend.background = element_rect(color = "black"))
ggsave("Figures/barplotc_sub.pdf", width = 16, height = 11)

#### ENV ####
nrank = envloc1 %>% 
  dplyr::mutate(rank = row_number(-ENV))# change here for comp
nrank$NONE = 0
envflip = tidyr::gather(nrank, "Type", "value", 2:5)
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)])
envflip = dplyr::arrange(envflip,rank)

envrank = envflip %>% 
  dplyr::group_by(Type == 'ENV') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "ENV") # change here for comp

###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE", "SHARED","COMP","ENV"),ordered = TRUE)
envflip$value = abs(envflip$value)
# Plot with ENV ranked in decreasing order - had to flip everything to plot right
e = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE","SHARED","COMP", "ENV")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#dd1c77","#2ca25f"), labels=c("","Shared Variance", "Competition","Environment")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30, angle=90),legend.title=element_blank(), legend.text=element_text(size=60, hjust = 1, vjust = 1), legend.position = c(0.7,0.9))+ guides(fill = guide_legend(keywidth = 2, keyheight = 4, legend.key.size = unit(5,"line")),title="")

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

w = ggplot(data=env_sub, aes(factor(rank), y=abs(value), fill=factor(Type, levels = c("NONE","SHARED","COMP", "ENV")))) + geom_bar(stat = "identity") + theme_classic() + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#dd1c77","#2ca25f"), labels=c("","Shared", "Competition","Environment")) +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=30, vjust =3),legend.title=element_blank(), legend.text=element_text(size=34, hjust = 1, vjust = 0.5), legend.position = c(.8,.9)) +theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.key.width=unit(4, "point"), legend.key.height =unit(3, "lines")) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 30,color= "black")) +scale_y_continuous(breaks = c(0,   0.2,  0.4, 0.6),limits = c(0, 0.6),sec.axis = sec_axis(~ . *1/1)) + annotate("text", x = 1:15, y = -.01, label = envrank$ALPHA.CODE[1:15], angle=90,size=6,vjust=0.5,hjust = 1, color = "black")
  
 
ggsave("Figures/barplote_sub.pdf", height = 11, width = 16)

##################### TRAITS Model ####################################
logit = function(x) log(x/(1-x))

env_lm = subset(envflip, Type == 'ENV')
comp_lm = subset(envflip, Type == 'COMP')

env_sum = subset(envflip, Type != 'NONE')
total = env_sum %>% 
  dplyr::group_by(FocalAOU) %>%
  summarise(sum(value))

# creating env traits model to compare to comp and weighted traits mods
env_cont = merge(env_lm, shapefile_overlap, by.x = "FocalAOU",by.y = "focalAOU")
env_cont2 = merge(env_cont, unique(occuenv[,c("FocalAOU", "Mean.Temp","Mean.Precip","Mean.Elev","Mean.NDVI")]), by.x = "FocalAOU", by.y = "FocalAOU")

# rescaling all occupancy values  - odds ratio
# need to get rid of ones in order to not have infinity values 
edge_adjust = .005 
env_cont2$COMPSC_sc = env_cont2$COMPSC * (1 - 2*edge_adjust) + edge_adjust
# create logit transformation function, did on rescaled vals
env_cont2$COMPSC_logit =  log(env_cont2$COMPSC_sc/(1-env_cont2$COMPSC_sc)) 
combined_mod = filter(env_cont2, Trophic.Group != "nectarivore" & Trophic.Group != "herbivore")

#### 5/7 CHANGED to be one big model incl trophic and mig
econt = lm(COMPSC_logit ~ log10(FocalArea)  + prop_overlap + Mean.Temp + Mean.Precip + Mean.Elev + Mean.NDVI + Trophic.Group + migclass, data = combined_mod, weights = n)
env_est = summary(econt)$coef[,"Estimate"]
colname = c("Intercept","FocalArea", "area_overlap","Mean.Temp","Mean.Precip","Mean.Elev", "Mean.NDVI", "Trophic.Groupinsct/om","Trophic.Groupinsectivore", "Trophic.Groupomnivore", "migclassresid", "migclassshort")
env = data.frame(colname, env_est)
# add the constant here
env$env_lower =  as.vector(summary(econt)$coefficients[,"Estimate"]) - 1.96*as.vector(summary(econt)$coef[,"Std. Error"])
env$env_upper = as.vector(summary(econt)$coefficients[,"Estimate"]) + 1.96*as.vector(summary(econt)$coef[,"Std. Error"])

env_trait_rank = env %>% 
  dplyr::mutate(rank = row_number(-env_est)) 
env_trait_rank2 <- env_trait_rank[order(env_trait_rank$rank),]

env_trait_rank = env %>% 
  dplyr::mutate(rank = row_number(-env_est)) 
env_trait_rank2 <- env_trait_rank[order(env_trait_rank$rank),]
env_trait_rank2$colname = factor(env_trait_rank2$colname,
     levels = c("Intercept","Mean.NDVI","Trophic.Groupinsct/om","Trophic.Groupinsectivore","migclassshort","area_overlap","Mean.Precip","Mean.Elev","Mean.Temp","Trophic.Groupomnivore","migclassresid", "FocalArea"),ordered = TRUE)

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
       levels = c("Insectivore","Insectivore/\nOmnivore","Granivore","Omnivore","Herbivore"),ordered = TRUE)

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

tukeys = aov(lm(COMPSC_logit ~ log10(FocalArea)  + prop_overlap + Mean.Temp + Mean.Precip + Mean.Elev + Mean.NDVI + Trophic.Group + migclass, data = combined_mod, weights = n))
TukeyHSD(tukeys)
summary(tukeys)



suppl = merge(env_lm, nsw[,c("CompAOU", "focalAOU", "Competitor", "Focal")], by.x = "FocalAOU", by.y = "focalAOU")
# write.csv(suppl, "data/suppl_table.csv", row.names = FALSE)
# anova of traits
cor.test(envoutput$ENV, envoutputa$ENV)





table_s2 = left_join(envloc1, unique(nsw[,c("focalAOU", "FocalMass")]), by = c("FocalAOU" = "focalAOU")) %>%
  left_join(., unique(shapefile_overlap), by = c("FocalAOU" = "focalAOU")) %>%
  left_join(., unique(occuenv[,c("FocalAOU", "Mean.Temp","Mean.Precip","Mean.Elev","Mean.NDVI")]), by = "FocalAOU") %>% left_join(., sppGT50rtes, by = "FocalAOU") 


table_s2_cont = occuenv %>%
  group_by(FocalAOU, Focal) %>%
  summarise(Median_Occ = median(FocalOcc),
            Variance_Occ = var(FocalOcc),
            Min_Occ = min(FocalOcc),
            Max_Occ = max(FocalOcc))
# write.csv(table_s2_cont, "data/table_s2_cont.csv", row.names = FALSE)

# write.csv(table_s2, "data/table_s2.csv", row.names = FALSE)
TableS2 <- read.csv("Z:/Snell/2019 BI MS/Tables/Table S2 Traits.csv", header = TRUE) %>%
  left_join(., sppGT50rtes, by = c("Focal.Common.Name" = "Focal")) 
# write.csv(TableS2, "data/TableS2.csv", row.names = FALSE) 

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
R2plot2$Environment = R2plot2$ENV.x + R2plot2$SHARED.x
R2plot2$Competition = R2plot2$COMP.x + R2plot2$SHARED.x
R2plot2$Total = R2plot2$ENV.x + R2plot2$COMP.x + R2plot2$SHARED.x

R2plot2$env_a = R2plot2$ENV.y + R2plot2$SHARED.y
R2plot2$comp_a = R2plot2$COMP.y + R2plot2$SHARED.y

# need to change the slopes
cols = c("Competition" ="#dd1c77","Environment" = "#2ca25f","Total" = "dark gray")
r1 = ggplot(R2plot2, aes(x = comp_o, y = comp_a, col = "Competition")) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy R"^"2")) + ylab(bquote("Abundance R"^"2"))+
  geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="longdash", lwd =2.5) +
      geom_point(data = R2plot2, aes(x = env_o, y = env_a, col = "Environment"), shape = 16, cex =4, stroke = 1)+geom_smooth(data = R2plot2, aes(x = env_o, y = env_a), method='lm', se=FALSE, col="#2ca25f",linetype="longdash", lwd = 2.5) + 
      geom_point(data = R2plot2, aes(Total.x,Total.y, col = "Total"), shape = 3, cex =5, stroke = 1)+geom_smooth(data = R2plot2, aes(x =Total.x, y = Total.y), method='lm', se=FALSE, col="dark gray",linetype="longdash", lwd =2.5) + xlim(c(0, 0.8))+ theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.position = c(0.8,0.2), legend.key.width=unit(2, "lines"), legend.key.height =unit(3, "lines")) + scale_y_continuous(limits = c(0, 0.8), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1)) #+geom_label(data = R2plot2, aes(Total.x,Total.y, label = FocalAOU))
ggsave("C:/Git/Biotic-Interactions/Figures/occvabun_lines.pdf", height = 8, width = 12)

leg4 <- gather(R2plot2, "Occ", "val", c(Environment,Competition,Total))
ggplot(leg4, aes(x = Occ, y = val, shape = Occ, color = Occ)) + geom_point() + theme_classic() + scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray")) + scale_shape_manual("", values = c(24, 16, 3))+theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.position = c(0.5,0.8), legend.key.width=unit(2, "lines"), legend.key.height =unit(3, "lines")) + guides(shape = guide_legend(override.aes = list(size = 7)))

R2plot2$occdiff = R2plot2$COMP.x - R2plot2$ENV.x
R2plot2$abundiff = R2plot2$COMP.y - R2plot2$ENV.y
R2plot2$totaldiff = R2plot2$abundiff - R2plot2$occdiff

#### Figure 1 violin plots ####
# Figure 1B is in suppl script
# To get supplement, change x to y and use the correct COMPSC!
R2plot2$COMPSC = R2plot2$COMP.x/(R2plot2$COMP.x+R2plot2$ENV.x)
# R2violin.5 = left_join(R2plot2[,c("FocalAOU", "violin_env","violin_comp","violin_total")], envloc[,c("FocalAOU", "COMPSC")], by = c("FocalAOU" = "FocalAOU"))

# create Table S3
Table_S3 = left_join(R2plot2, tax_code, by =c ("FocalAOU" = "AOU_OUT"))
# write.csv(Table_S3, "data/Table_s3.csv", row.names = FALSE) 

R2plot2$Total.x = R2plot2$ENV.x + R2plot2$COMP.x + R2plot2$SHARED.x
R2violin.5 = left_join(R2plot2[,c("FocalAOU", "ENV.x","COMP.x","Total.x")], envloc[,c("FocalAOU", "COMPSC")], by = c("FocalAOU" = "FocalAOU"))

R2violin = gather(R2violin.5, "type", "Rval", 2:5)

R2violin$type = factor(R2violin$type,
                              levels = c("ENV.x","COMP.x","Total.x","COMPSC" ),ordered = TRUE)

# Figure 1B in Supplemental Material
ggplot(R2violin, aes(as.factor(type), Rval)) + geom_violin(linetype = "blank", scale ="count", aes(fill = factor(R2violin$type))) + xlab("") + ylab(bquote("Variance Explained"))+scale_fill_manual(values=c("#2ca25f","#dd1c77", "grey", "#636363"), labels=c("Environment","Competition", "Total Variance", "Scaled \nCompetition")) + theme_classic()+theme(axis.title.x=element_text(size=30, angle = 180),axis.title.y=element_text(size=30, vjust = 4))+scale_y_continuous(limits = c(0, 1)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size=25, color = "black"),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) +  guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))  + stat_summary(aes(group=factor(R2violin$type)), fun.y=median, geom="point",fill="black", shape=16, size=4,position = position_dodge(width = .9)) 
ggsave("Figures/violin_all.pdf", height = 8, width = 12)

##### Figure 6 non-competitor comparison ######
noncompdf = occuenv[,c("FocalAOU", "stateroute", "FocalOcc", "FocalAbundance", "Family", "FocalOcc_scale", "occ_logit")]
subfocspecies = unique(noncompdf$FocalAOU)
noncomps = nsw[,c("CompAOU", "Family")]
noncomps = unique(noncomps)


noncomps_output = c()
# for (sp in subfocspecies){
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
      mergespp$all_comp_scaled = mergespp$meancomp/(mergespp$FocalAbundance + mergespp$meancomp) 
      mergespp$all_comp_scaled[is.na(mergespp$all_comp_scaled)] = 0
      
      if(length(unique(mergespp$all_comp_scaled[!is.na(mergespp$all_comp_scaled)])) > 2){ 
        lms <- lm(mergespp$occ_logit ~  mergespp$all_comp_scaled) 
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

# write.csv(noncomps_output, "data/noncomps_output_all.csv", row.names = FALSE)

noncomps_output = read.csv("data/noncomps_output.csv", header = TRUE)
# filtering to species where p <0.05
beta_occ_main = read.csv("Z:/Snell/2019 BI MS/Tables/Table S8 beta occupancy widespread.csv", header = TRUE)
comp_abun_data <- read.csv("Z:/Snell/2019 BI MS/data/comp_abun_data.csv", header = TRUE)
noncomps_sub = left_join(beta_occ_main, tax_code, by = c("Focal.Common.Name" = "PRIMARY_COM_NAME"))
noncomps_sub2 = filter(noncomps_sub, Competitor.R2 >= 0.1) 
noncomps_sub3 = unique(noncomps_sub2)
noncomps_plot = subset(noncomps_output, FocalAOU %in% noncomps_sub3$AOU_OUT)

noncomps_output_bocc = left_join(noncomps_plot, noncomps_sub[,c("AOU_OUT", "Competitor.R2", "Competitor.Estimate", "Competitor.P.value")], by = c("FocalAOU" = "AOU_OUT"))

comp_abun_data$FocalOcc_scale = (comp_abun_data$occ * (1 - 2*edge_adjust)) + edge_adjust
# create logit transformation function, did on rescaled vals
comp_abun_data$occ_logit =  log(comp_abun_data$FocalOcc_scale/(1-comp_abun_data$FocalOcc_scale)) 
comp_abun_data$comp_abun[is.na(comp_abun_data$comp_abun)] <- 0
comp_abun_data$occ_logit[is.na(comp_abun_data$occ_logit)] <- 0

beta_occ_abun = data.frame(FocalAOU = c(), FocalSciName = c(), CompAOU = c(), CompSciName = c(), Estimate = c(), P = c(), R2 = c()) 
for(i in unique(envoutput$FocalAOU)){
  print(i)
  temp = subset(comp_abun_data, focalAOU == i) 
  for(j in unique(temp$CompAOU)){
    ctemp = subset(temp, CompAOU == j)
    com <- unique(ctemp$CompSciName)
    foc <- unique(ctemp$FocalSciName)
    length(na.omit(ctemp$abundance.y))
    if(sum(ctemp$comp_abun) > 2){
      competition <- lm(ctemp$occ_logit ~  ctemp$comp_abun)  # changes between main and all comps
      occ_comp_est = summary(competition)$coef[2,"Estimate"]
      occ_comp_p = summary(competition)$coef[2,"Pr(>|t|)"]
      occ_comp_r = summary(competition)$r.squared
     beta_occ_abun = rbind(beta_occ_abun, data.frame(FocalAOU = i, FocalSciName = foc, CompAOU = j, CompSciName = com, Estimate = occ_comp_est, P = occ_comp_p, R2 = occ_comp_r))
    } 
  }
}
# write.csv(beta_occ_abun, "data/beta_occ_abun.csv", row.names = FALSE)



nonps = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(R2 >= Competition.R2)
names(nonps) = c("FocalAOU", "main_g_non")

none = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(Estimate <= Competition.Estimate)
names(none) = c("FocalAOU", "main_g_non_e")

numcomps = na.omit(noncomps_output_bocc) %>% 
  count(FocalAOU)
names(numcomps) = c("FocalAOU", "Comp_count")

vec = c()
for(i in unique(noncomps_output_bocc$FocalAOU)){
  temp = subset(noncomps_output_bocc, FocalAOU == i)
  vec = append(vec, median(temp$R2))
} 

vec_meds = data.frame(FocalAOU = unique(noncomps_output_bocc$FocalAOU), Median = vec)  
noncomps_output_ttest = left_join(noncomps_output_bocc, vec_meds, by = "FocalAOU")
t.test(noncomps_output_ttest$Competition.R2, noncomps_output_ttest$Median, paired = TRUE, alternative= "two.sided")

noncompsdistp  = merge(nonps, numcomps, by = ("FocalAOU"))
noncompsdiste  = merge(none, numcomps, by = ("FocalAOU"))
noncompsdistp$nullp = (noncompsdistp$main_g_non)/(noncompsdistp$Comp_count + 1)
noncompsdiste$Widespread_e = (noncompsdiste$main_g_non_e)/(noncompsdiste$Comp_count + 1)



##### for main comp analysis #####
nonps = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(R2 >= Competitor.R2)
names(nonps) = c("FocalAOU", "main_g_non")

none = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(Estimate <= Competitor.Estimate)
names(none) = c("FocalAOU", "main_g_non_e")

numcomps = na.omit(noncomps_output_bocc) %>% 
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

noncompsdistp  = merge(nonps, numcomps, by = ("FocalAOU"))
noncompsdiste  = merge(none, numcomps, by = ("FocalAOU"))
noncompsdistp$Widespread_p = (noncompsdistp$main_g_non)/(noncompsdistp$Comp_count + 1)
noncompsdiste$Widespread_e = (noncompsdiste$main_g_non_e)/(noncompsdiste$Comp_count + 1)

noncomps_output_bocc$Null = "Null"
noncomps_output_bocc$Comp = "Comp"
noncompsdistp$Null = "Null"
noncompsdiste$Null = "Null"

noncompsdiste_posthoc <- read.csv("data/noncompsdiste_posthoc.csv", header = TRUE)
noncompsdistp_posthoc <- read.csv("data/noncompsdistp_posthoc.csv", header = TRUE)


noncompdist_pcombined <- full_join(noncompsdistp[,c("FocalAOU", "Widespread_p")], noncompsdistp_posthoc[,c("FocalAOU", "nullp")], by = "FocalAOU") %>%
  gather(p_val, value, Widespread_p:nullp)

noncompdist_ecombined <- full_join(noncompsdiste[,c("FocalAOU", "Widespread_e")], noncompsdiste_posthoc[,c("FocalAOU", "nulle")], by = "FocalAOU") %>%
  gather(e_val, value, Widespread_e:nulle)

#### Figure 6 example non-comp dist and main R2 ######
single_dist = subset(noncomps_output_bocc, FocalAOU == 4020)
n = ggplot(single_dist) +
  geom_histogram(bins = 15, aes(R2, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  geom_vline(xintercept = single_dist$Competitor.R2, col = "black", lwd = 1.5, lty = 2) +
  xlab(expression("Variance Explained")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("#c994c7"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"))

o = ggplot(single_dist) +
  geom_histogram(bins = 15, aes(Estimate, fill=factor(Null, levels = c("Null"))), alpha = 0.9) +
  geom_vline(xintercept = single_dist$Competitor.Estimate, col = "black", lwd = 1.5, lty = 2) +
  xlab(expression("Competitor Estimate")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("#c994c7"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"))

p = ggplot(noncompdist_pcombined, aes(x=value,fill=p_val)) +
  geom_histogram(bins = 15, position = "identity", alpha = 0.9) +
  xlab(expression('Proportion of non-competitors with higher R'^2)) + ylab("Frequency") + theme_classic()  +
  scale_fill_manual(values = c("#330066","#756bb1"), labels = c("Maximum Non-competitors", "Widespread Non-competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 24)) + theme(axis.title.x=element_text(size=24, vjust = -1.5),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"), legend.position = c(.5, .7))


q = ggplot(noncompdist_ecombined, aes(x=value,fill=e_val)) +
  geom_histogram(bins = 15, position = "identity", alpha = 0.9) + 
  xlab("Proportion of non-competitors \n with more negative slope") + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(values = c("#330066","#756bb1"), labels = c("Predictive", "Widespread")) + theme(legend.title=element_blank(), legend.text=element_text(size = 20)) + theme(axis.title.x=element_text(size=24, vjust = -5),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"))


legend <- get_legend(q)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(n + theme(legend.position="none"),
          o + theme(legend.position="none"),
          p + theme(legend.position= "none"),
          q + theme(legend.position=c(0.5, 0.7)),
          align = 'hv',
          labels = c("a","b", "c", "d"),
          label_size = 20,
          nrow = 2) 

ggsave("Figures/Figure6_null_all.pdf", height = 8, width = 12)

# NOTE: All figures were created using this script, but final edits were made in powerpoint. See ppt file in repository for final figures. 