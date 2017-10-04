library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggExtra)

# read in files created in data cleaning script
temp_occ = read.csv("data/bbs_sub1.csv", header=TRUE) # BBS occ script
centroid=read.csv("data/centroid.csv", header=TRUE) # GIS script
occuenv= read.csv("data/all_expected_pres.csv", header = TRUE) # Data cleaning script
subsetocc = read.csv('data/subsetocc.csv', header = T) # Hurlbert Lab
tax_code = read.csv("data/Tax_AOU_Alpha.csv", header = TRUE) # Hurlbert Lab
#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
# rbind Pacific wren to data frame
pacific = data.frame("Pacific Wren", 7221,  "PAWR")
colnames(pacific) = c("PRIMARY_COM_NAME", "AOU_OUT" ,"ALPHA.CODE")
tax_code = rbind(tax_code, pacific)

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

subfocalspecies = unique(occuenv$Species)

for (sp in 1:length(subfocalspecies)){
  temp = subset(occuenv,occuenv$Species == subfocalspecies[sp])
  
  competition <- lm(temp$occ_logit ~  temp$comp_scaled)  # changes between main and all comps
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
  ENV = summary(both_z)$r.squared - summary(competition)$r.squared #env only
  COMP = summary(both_z)$r.squared - summary(env_z)$r.squared #competition only
  SHARED = summary(competition)$r.squared - COMP #shared variance
  NONE = 1 - summary(both_z)$r.squared # neither variance
  sp1 = unique(temp$Species)
  sum = sum(ENV, COMP, SHARED)
  envoutput = rbind(envoutput, c(sp1, ENV, COMP, SHARED, NONE)) #, sum
  
  #variance_partitioning 
  ENVa = summary(both_abun)$r.squared - summary(competition_abun)$r.squared
  COMPa = summary(both_abun)$r.squared - summary(env_abun)$r.squared
  SHAREDa = summary(competition_abun)$r.squared - COMP
  NONEa = 1 - summary(both_abun)$r.squared
  
  sp1 = unique(temp$Species)
  envoutputa = rbind(envoutputa, c(sp1, ENVa, COMPa, SHAREDa, NONEa))
  # saving model output into separate data frames
  occ_comp_est = summary(competition)$coef[1,"Estimate"]
  occ_comp_p = summary(competition)$coef[1,"Pr(>|t|)"]
  occ_comp_r = summary(competition)$r.squared
  occ_env_est = summary(env_z)$coef[2,"Estimate"]
  occ_env_p = summary(env_z)$coef[2,"Pr(>|t|)"]
  occ_env_r = summary(env_z)$r.squared 
  occ_b_est = summary(both_z)$coef[2,"Estimate"]
  occ_b_p = summary(both_z)$coef[2,"Pr(>|t|)"]
  occ_b_r = summary(both_z)$r.squared 

  abun_comp_est = summary(competition_abun)$coef[1,"Estimate"]
  abun_comp_p = summary(competition_abun)$coef[1,"Pr(>|t|)"]
  abun_comp_r = summary(competition_abun)$r.squared #using multiple rsquared
  abun_env_est = summary(env_abun)$coef[2,"Estimate"]
  abun_env_p = summary(env_abun)$coef[2,"Pr(>|t|)"]
  abun_env_r = summary(env_abun)$r.squared 
  abun_both_est = summary(both_abun)$coef[2,"Estimate"]
  abun_both_p = summary(both_abun)$coef[2,"Pr(>|t|)"]
  abun_both_r = summary(both_abun)$r.squared
  
  beta_occ = rbind(beta_occ, c(sp1, occ_comp_est, occ_comp_p, occ_comp_r, occ_env_est, occ_env_p, occ_env_r,occ_b_est, occ_b_p , occ_b_r))
  beta_abun = rbind(beta_abun, c(sp1, abun_comp_est, abun_comp_p, abun_comp_r, abun_env_est, abun_env_p, abun_env_r,abun_both_est, abun_both_p , abun_both_r))
  
}         


envoutput = data.frame(envoutput)
envoutputa = data.frame(envoutputa)

beta_occ = data.frame(beta_occ)
beta_abun = data.frame(beta_abun)

names(envoutput) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")
names(envoutputa) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")

envoutput1 = merge(envoutput, tax_code[,c('AOU_OUT', 'ALPHA.CODE')], by.x = 'FocalAOU', by.y = "AOU_OUT", all.x = TRUE) 

envoutput2 = merge(envoutput, subsetocc[,c("AOU", "migclass", "Trophic.Group")], by.x='FocalAOU', by.y='AOU', all.x = TRUE)

envloc = merge(envoutput2, centroid[, c("FocalAOU", "Long", "Lat")], by = 'FocalAOU', all.x = TRUE)
### supp table goes here, just need to add main competitor!


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
occumatrix$sp_success = 15 * occumatrix$FocalOcc
occumatrix$sp_fail = 15 * (1 - occumatrix$FocalOcc)

#### GLM of all matrices not just subset ####
glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ c_s + 

                            abTemp + abElev + abPrecip + abNDVI + (1|stateroute:Species), family = binomial(link = logit), data = occumatrix)
summary(glm_occ_rand_site)                                    

#### new fig 1 ####
fig1 = ggplot(data = occuenv, aes(x = log10(FocalAbundance), y = FocalOcc)) +geom_point() + geom_jitter(width = 0, height = 0.02) +xlab("log10(Focal Abundance)")+ylab("Focal Occupancy") + geom_hline(yintercept = median(occuenv$FocalOcc), lwd = 1, col = "red")+ geom_vline(xintercept = median(log10(occuenv$FocalAbundance)), lwd = 1, col = "red") +theme_classic() +theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24, angle=90), axis.text=element_text(size=12)) + theme(plot.margin = unit(c(.5,6,.5,.5),"lines")) 
ggExtra::ggMarginal(fig1 , type = "histogram", fill = "dark gray")
ggsave("C:/Git/Biotic-Interactions/Figures/fig1.pdf")

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
  dplyr::mutate(rank = row_number(-COMP))# change here for comp
envflip = tidyr::gather(nrank, "Type", "value", 2:5)
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)])
envflip = dplyr::arrange(envflip,rank)

# envflip = merge(envflip, envloc[,c("FocalAOU", "EW")], by = "FocalAOU")

envrank = envflip %>% 
  dplyr::group_by(Type == 'COMP') %>% # change here for comp
  dplyr::mutate(rank = row_number(-value)) # need to get just the envs to rank, then plot
envrank <- envrank[order(envrank$rank),]

envrank <- subset(envrank,Type == "COMP") # change here for comp

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

envrank$EW[envrank$EW == 1] <- "E"
envrank$EW[envrank$EW == 0] <- "W" 
###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE", "SHARED","ENV","COMP"),ordered = TRUE)
envflip$value = abs(envflip$value)
# Plot with ENV ranked in decreasing order - had to flip everything to plot right
t = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE", "SHARED","ENV","COMP")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(angle=90,size=10,vjust=0.5)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","yellow3","#2ca25f","#dd1c77"), labels=c("","Shared Variance", "Environment", "Competition")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30, angle=90),legend.title=element_blank(), legend.text=element_text(size=60), legend.position = "top",legend.key.width=unit(1, "lines")) + guides(fill=guide_legend(fill = guide_legend(keywidth = 2, keyheight = 2),title=""))

tt = t + annotate("text", x = 1:104, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=0.5, color = "black") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks=scales::pretty_breaks()(0:1))

#+ annotate("text", x = 1:104, y = -.08, label = envrank$mig_abbrev, size=6,vjust=0.5, color = envrank$mig_abbrevf, fontface =2) + annotate("text", x = 1:104, y = -.1, label = envrank$trophlabel, size=6,vjust=0.5, color = envrank$trophlabelf, fontface =2) + annotate("text", x = 1:104, y = -.12, label = envrank$EW, angle=90,size=6,vjust=0.5, color = "black", fontface =2)+ annotate("text", x = 1:104, y = -.06, label = envrank$Fam_abbrev, size=6,vjust=0.5, color = "black", fontface =2) 

plot(tt)

ggsave("C:/Git/Biotic-Interactions/Figures/barplotc.pdf", height = 25, width = 36)

##################### TRAITS Model ####################################
logit = function(x) log(x/(1-x))

env_lm = subset(envflip, Type == 'ENV')

env_traits = lm(logit(value) ~ Trophic.Group + migclass + EW, data = env_lm)
summary(env_traits) 

comp_lm = subset(envflip, Type == 'COMP')

comp_traits = lm(logit(value) ~ Trophic.Group + migclass + EW, data = comp_lm[comp_lm$value > 0,])
summary(comp_traits) 

env_sum = subset(envflip, Type != 'NONE')
total = env_sum %>% 
  dplyr::group_by(FocalAOU) %>%
  summarise(sum(value))

env_sum$edgeval = (env_sum$value * (1 - 2*edge_adjust)) + edge_adjust
total_traits = lm(logit(edgeval) ~ Trophic.Group + migclass + EW, data = env_sum)
summary(total_traits)

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
r1 = ggplot(R2plot2, aes(x = COMP.x, y = COMP.y)) +theme_classic()+ theme(axis.title.x=element_text(size=26),axis.title.y=element_text(size=26, angle=90)) + xlab("Occupancy R2") + ylab("Abundance R2") + geom_point(col = "#dd1c77", cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash") +
      geom_point(data = R2plot2, aes(x = ENV.x, y = ENV.y), shape = 16, col = "#2ca25f", cex =4, stroke = 1)+geom_smooth(data = R2plot2, aes(x = ENV.x, y = ENV.y), method='lm', se=FALSE, col="#2ca25f",linetype="dotdash") +
      geom_point(data = R2plot2, aes(Total.x,Total.y), shape = 3, col = "dark gray", cex =5, stroke = 1)+geom_smooth(data = R2plot2, aes(x =Total.x, y = Total.y), method='lm', se=FALSE, col="dark gray",linetype="dotdash") +
      geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))
ggsave("C:/Git/Biotic-Interactions/Figures/occvabun_lines.png")

R2plot2$occdiff = R2plot2$COMP.x - R2plot2$ENV.x
R2plot2$abundiff = R2plot2$COMP.y - R2plot2$ENV.y

r2 = ggplot(R2plot2, aes(x = occdiff, y = abundiff)) +theme_classic()+ geom_abline(intercept = 0, slope = 0, col = "black", lwd = 1.25, lty = "dashed") + geom_vline(xintercept = 0, col = "black", lwd = 1.25, lty = "dashed")+ theme(axis.title.x=element_text(size=26),axis.title.y=element_blank()) + xlab("Occupancy R2") + geom_point(col = "black", shape=16, size = 3)+ theme(axis.text.x=element_text(size = 20),axis.ticks=element_blank(), axis.text.y=element_text(size=20))

p2 = plot_grid(r1 + theme(legend.position="none"),
               r2 + theme(legend.position="none"), 
               labels = c("A","B"),
               align = 'h')
ggsave("C:/Git/Biotic-Interactions/Figures/Figure4A_B.pdf")


# r2 plot for main vs all competitors
envmain = read.csv("data/envoutput.csv", header = TRUE)
mainvall = merge(envoutput, envmain, by = "FocalAOU")

ggplot(mainvall, aes(x = COMP.y, y = COMP.x)) +theme_bw()+ theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16, angle=90)) + xlab("main R2") + ylab("all R2") + geom_point(col = "#dd1c77", cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="dotdash") + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.25)


# R2 plot - glm violin plots
R2violin = gather(R2plot2, "type", "Rval", 12:14)

ggplot(R2violin, aes(as.factor(type), Rval)) + geom_violin(linetype = "blank", aes(fill = factor(R2violin$type))) + xlab("Total Variance") + ylab("R2")+scale_fill_manual(values=c("#2ca25f","#dd1c77", "grey"), labels=c("Environment", "Competition", "Total Variance")) + theme_bw()+theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))+scale_y_continuous(limits = c(0, 0.8)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(), axis.text.y=element_text(size=25),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) + guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))

ggsave("C:/Git/Biotic-Interactions/Figures/violin_mains.png")

#Coyle fig 1: Z:\Coyle\Projects\BBS Core\Final Analysis

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
  psub = occumatrix[occumatrix$Species == sp,]
  glm_occ_rand_site = glmer(cbind(sp_success, sp_fail) ~ comp_scaled + 
                              abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI) + (1|stateroute:Species), family = binomial(link = logit), data = psub)
  
  tes = ggplot(data = psub, aes(x = zTemp, y = FocalOcc)) +stat_smooth(method = "lm", lwd = 1.5,se = FALSE) +xlab(psub$Species) +ylim(0,1) +theme_bw()
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
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zElev)), yend = 1 +(-0.0144*max(abs(occumatrix$zElev)))), col = "dark green", lwd=2)  + 
  geom_point(colour="black", shape=18, alpha = 0.1,position=position_jitter(width=0,height=.02))+ theme_classic()
ggsave("C:/Git/Biotic-Interactions/Figures/logitelev.png")

precip = ggplot(data = occumatrix, aes(x = abs(zPrecip), y = FocalOcc)) + 
  geom_segment(aes(x = 0, y = 1, xend = abs(max(occumatrix$zPrecip)), yend = 1 +(0.005566*max(abs(occumatrix$zPrecip)))), col = "dark green", lwd=2) +
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


