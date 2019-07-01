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

post_hoc <- read.csv("data/post_hoc_comps.csv", header = TRUE)
post_hoc_output <- read.csv("data/envoutput_post_hoc_main.csv", header = TRUE)
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

l_focalspecies = inner_join(shapefile_areas, occuenv, c("focalAOU"="FocalAOU"))
subfocalspecies = unique(occuenv$FocalAOU) #[-99] excluding sage sparrow bc no good routes


beta_occ_all = read.csv("Z:/Snell/2019 BI MS/Tables/Table S5 beta occupancy all.csv", header = TRUE)
comp_abun_data <- read.csv("Z:/Snell/2019 BI MS/comp_abun_data.csv", header = TRUE) %>%
  left_join(occuenv[,c("stateroute","FocalAOU","zTemp", "zElev","zPrecip", "zNDVI")], by = c("focalAOU" = "FocalAOU", "stateroute" = "stateroute")) 

comp_abun_data$FocalOcc_scale = (comp_abun_data$occ * (1 - 2*edge_adjust)) + edge_adjust
# create logit transformation function, did on rescaled vals
comp_abun_data$occ_logit =  log(comp_abun_data$FocalOcc_scale/(1-comp_abun_data$FocalOcc_scale)) 
comp_abun_data$comp_abun[is.na(comp_abun_data$comp_abun)] <- 0
comp_abun_data$occ_logit[is.na(comp_abun_data$occ_logit)] <- 0
# Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) 
comp_abun_data$comp_scaled = comp_abun_data$comp_abun/(comp_abun_data$focal_abun + comp_abun_data$comp_abun) 
comp_abun_data$comp_scaled[is.na(comp_abun_data$comp_scaled)] = 0



subfocalspecies = unique(occuenv$FocalAOU)

allcomps_output = c()
for (sp in subfocalspecies){
  print("newspp")
  temp = filter(comp_abun_data, focalAOU == sp) 
  if(nrow(temp) > 0){
    ncomps = temp %>%
      dplyr::distinct(CompAOU) %>% 
      unlist()
    comps = unique(ncomps)
    for(co in comps){
      print(co)
      mergespp = filter(temp, CompAOU == co)
      
      if(length(unique(mergespp$comp_scaled[!is.na(mergespp$comp_scaled)])) > 2){ 
        lms <- lm(mergespp$occ_logit ~  mergespp$comp_scaled) 
        lms_est = summary(lms)$coef[2,"Estimate"]
        lms_p = summary(lms)$coef[2,"Pr(>|t|)"]
        lms_r = summary(lms)$r.squared
        
        allcomps_output = rbind(allcomps_output, c(sp, co,lms_est, lms_p, lms_r))
      }
    }
  }
}         

allcomps_output = data.frame(allcomps_output)
names(allcomps_output) = c("FocalAOU", "CompetitorAOU", "Estimate","P", "R2")
# write.csv(allcomps_output, "data/allcomps_output.csv", row.names = FALSE)

highestR2 <- left_join(allcomps_output, maincomp1[,c("focalAOU", "compAOU", "mainCompetitor")], by = c("FocalAOU" = "focalAOU", "CompetitorAOU" = "compAOU")) %>%
  group_by(FocalAOU) %>%
  arrange(desc(R2)) %>% 
  # Pick the top 1 value
  slice(1) %>%
  # Remember to ungroup in case you want to do further work without grouping.
  ungroup()

highestR2$post_hoc_main <- 1

neg_est_comps <- left_join(allcomps_output, maincomp1[,c("focalAOU", "compAOU", "mainCompetitor")], by = c("FocalAOU" = "focalAOU", "CompetitorAOU" = "compAOU")) %>%
  group_by(FocalAOU) %>%
  arrange(desc(R2)) %>% 
  filter(Estimate < 0) %>%
  slice(1) %>%
  # Remember to ungroup in case you want to do further work without grouping.
  ungroup()
neg_est_comps$maxCompetitor <- 1
neg_est_join <- left_join(neg_est_comps, tax_code, by = c("FocalAOU" = "AOU_OUT")) %>%
  left_join(tax_code, by = c("CompetitorAOU" = "AOU_OUT")) %>%
  left_join(nsw[,c("focalAOU", "CompAOU", "Family", "FocalSciName", "CompSciName")], by = c("FocalAOU" = "focalAOU", "CompetitorAOU" ="CompAOU")) %>%
  select(FocalAOU, CompetitorAOU, maxCompetitor, PRIMARY_COM_NAME.x, PRIMARY_COM_NAME.y)
Table_S1 <- read.csv("Z:/Snell/2019 BI MS/Tables/Table S1 Focal Competitors.csv", header = TRUE)
TableS1 <- left_join(Table_S1, neg_est_join, by = c("Focal.Common.Name" = "PRIMARY_COM_NAME.x", "Competitor.Common.Name" = "PRIMARY_COM_NAME.y", "Focal.AOU" = "FocalAOU", "Competitor.AOU" = "CompetitorAOU"))
# write.csv(TableS1, "data/TableS1.csv", row.names = FALSE)

envoutput = c()
envoutputa = c()
beta_occ_abun = data.frame(FocalAOU = c(), Focal = c(), CompAOU = c(), Comp = c(), Estimate = c(), P = c(), R2 = c(), envr2 = c(), bothp = c(), bothr2 = c()) 
for(i in unique(highestR2$FocalAOU)){
  print(i)
  comp = subset(neg_est_comps, FocalAOU == i)$CompetitorAOU
  temp = subset(comp_abun_data, focalAOU == i & CompAOU == comp) 
    if(sum(temp$comp_abun) > 0){
      competition <- lm(temp$occ_logit ~  temp$comp_scaled)  # changes between main and all comps
      env_z = lm(temp$occ_logit ~ abs(zTemp) + abs(zElev) + abs(zPrecip) + abs(zNDVI), data = temp)
      # z scores separated out for env effects
      both_z = lm(temp$occ_logit ~  temp$comp_scaled + abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
      
      # abundance, not temp occ - same results?
      comp_abun <- lm(temp$focal_abun ~  temp$comp_scaled) 
      # z scores separated out for env effects - abundance
      env_abun = lm(temp$focal_abun  ~ abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
      # z scores separated out for env effects - abundance
      both_abun = lm(temp$focal_abun  ~  comp_scaled + abs(zTemp)+abs(zElev)+abs(zPrecip)+abs(zNDVI), data = temp)
      
      sp1 = i
      # variance_partitioning 
      ENV = summary(both_z)$r.squared - summary(competition)$r.squared #env only
      COMP = summary(both_z)$r.squared - summary(env_z)$r.squared #competition only
      SHARED = summary(competition)$r.squared - COMP #shared variance
      NONE = 1 - summary(both_z)$r.squared # neither variance
      sum = sum(ENV, COMP, SHARED)
      n = length(temp$stateroute)
      envoutput = rbind(envoutput, c(sp1, ENV, COMP, SHARED, NONE, n)) #, sum
      
      # variance_partitioning 
      ENVa = summary(both_abun)$r.squared - summary(comp_abun)$r.squared
      COMPa = summary(both_abun)$r.squared - summary(env_abun)$r.squared
      SHAREDa = summary(comp_abun)$r.squared - COMPa
      NONEa = 1 - summary(both_abun)$r.squared
    
      envoutputa = rbind(envoutputa, c(sp1, ENVa, COMPa, SHAREDa, NONEa))
      
      
      occ_comp_est = summary(competition)$coef[2,"Estimate"]
      occ_comp_p = summary(competition)$coef[2,"Pr(>|t|)"]
      occ_comp_r = summary(competition)$r.squared
      abun_comp_est = summary(comp_abun)$coef[2,"Estimate"]
      abun_comp_p = summary(comp_abun)$coef[2,"Pr(>|t|)"]
      abun_comp_r = summary(comp_abun)$r.squared
      occ_env_r = summary(env_z)$r.squared 
      occ_b_p = summary(both_z)$coef[2,"Pr(>|t|)"]
      occ_b_r = summary(both_z)$r.squared 
      
      beta_occ_abun = rbind(beta_occ_abun, data.frame(FocalAOU = i, Focal = temp$Focal, CompAOU = comp, Comp = temp$Competitor, Estimate = occ_comp_est, P = occ_comp_p, R2 = occ_comp_r, envr2 = occ_env_r, bothp = occ_b_p, bothr2 = occ_b_r))
    } 
 }

beta_occ_abun <- beta_occ_abun %>%
  distinct()
# write.csv(beta_occ_abun, "data/beta_occ_abun_posthoc.csv", row.names = FALSE)


envoutput = data.frame(envoutput)
envoutputa = data.frame(envoutputa)

names(envoutput) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE", "n")
names(envoutputa) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")

envoutput$COMPSC = envoutput$COMP/(envoutput$COMP+envoutput$ENV)

suppl_table <- left_join(envoutput, beta_occ_abun[,c("FocalAOU", "Focal", "Comp")], by = "FocalAOU") %>%
  distinct()
# write.csv(suppl_table, "data/post_hoc_comps.csv", row.names = FALSE)

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
R2plot2$env_o = R2plot2$ENV.x + R2plot2$SHARED.x
R2plot2$comp_o = R2plot2$COMP.x + R2plot2$SHARED.x
R2plot2$total_o = R2plot2$ENV.x + R2plot2$COMP.x + R2plot2$SHARED.x

R2plot2$env_a = R2plot2$ENV.y + R2plot2$SHARED.y
R2plot2$comp_a = R2plot2$COMP.y + R2plot2$SHARED.y
R2plot2$total_a = R2plot2$ENV.y + R2plot2$COMP.y + R2plot2$SHARED.y


R2plot2$occdiff = R2plot2$COMP.x - R2plot2$ENV.x
R2plot2$abundiff = R2plot2$COMP.y - R2plot2$ENV.y
R2plot2$totaldiff = R2plot2$abundiff - R2plot2$occdiff

#### Figure 1 violin plots ####
# need to chance comp_scaled to all_comp scaled in envouputput loop for Fig 1B
R2plot2$COMPSC = R2plot2$COMP.x/(R2plot2$COMP.x+R2plot2$ENV.x)

# create Table S2
Table_S2 = left_join(R2plot2, maincomp1.5, by =c ("FocalAOU" = "focalAOU"))
Table_S2.1 = left_join(Table_S2, tax_code, by =c ("compAOU" = "AOU_OUT"))
# write.csv(Table_S2.1, "data/Table_s2.csv", row.names = FALSE) 

# R2violin.5 = left_join(R2plot2[,c("FocalAOU", "violin_env","violin_comp","violin_total")], envloc[,c("FocalAOU", "COMPSC")], by = c("FocalAOU" = "FocalAOU"))

R2violin.5 = R2plot2[,c("FocalAOU", "ENV.x","COMP.x","Total.x", "COMPSC")]
R2violin = gather(R2violin.5, "type", "Rval", 2:5)

R2violin$type = factor(R2violin$type,
                       levels = c("ENV.x","COMP.x","Total.x","COMPSC"),ordered = TRUE)

ggplot(R2violin, aes(as.factor(type), Rval)) + geom_violin(linetype = "blank",scale = "count", aes(fill = factor(R2violin$type))) + xlab("") + ylab(bquote("Variance Explained"))+scale_fill_manual(values=c("#2ca25f","#dd1c77", "grey", "#636363"), labels=c("Environment","Competition", "Total Variance", "Scaled \nCompetition")) + theme_classic()+theme(axis.title.x=element_text(size=30, angle = 180),axis.title.y=element_text(size=30, vjust = 3))+scale_y_continuous(limits = c(0, 1)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size=25, color = "black"),legend.title=element_blank(), legend.text=element_text(size=27), legend.position = "top",legend.key.width=unit(1, "lines")) +  guides(fill=guide_legend(fill = guide_legend(keywidth = 3, keyheight = 1),title=""))  + stat_summary(aes(group=factor(R2violin$type)), fun.y=median, geom="point",fill="black", shape=16, size=3, position = position_dodge(width = .9)) 
# , sec.axis = sec_axis(~ . *1/1, name = "Variance Ratio")
ggsave("Figures/violin_post_hoc.pdf", height = 8, width = 12)


# need to change the slopes
cols = c("Competition" ="#dd1c77","Environment" = "#2ca25f","Total" = "dark gray")
r1 = ggplot(R2plot2, aes(x = comp_o, y = comp_a, col = "Competition")) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy R"^"2")) + ylab(bquote("Abundance R"^"2"))+
  geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="#dd1c77",linetype="longdash", lwd =2.5) +
  geom_point(data = R2plot2, aes(x = env_o, y = env_a, col = "Environment"), shape = 16, cex =4, stroke = 1)+geom_smooth(data = R2plot2, aes(x = env_o, y = env_a), method='lm', se=FALSE, col="#2ca25f",linetype="longdash", lwd = 2.5) + 
  geom_point(data = R2plot2, aes(Total.x,Total.y, col = "Total"), shape = 3, cex =5, stroke = 1)+geom_smooth(data = R2plot2, aes(x =Total.x, y = Total.y), method='lm', se=FALSE, col="dark gray",linetype="longdash", lwd =2.5) + xlim(c(0, 0.8))+ theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.position = c(0.8,0.2), legend.key.width=unit(2, "lines"), legend.key.height =unit(3, "lines")) + scale_y_continuous(limits = c(0, 0.8), breaks = c(0,0.2, 0.4, 0.6, 0.8, 1)) #+geom_label(data = R2plot2, aes(Total.x,Total.y, label = FocalAOU))
ggsave("C:/Git/Biotic-Interactions/Figures/occvabun_lines_posthoc.pdf", height = 8, width = 12)




# write.csv(noncomps_output, "data/noncomps_output.csv", row.names = FALSE)
noncomps_output = read.csv("data/noncomps_output.csv", header = TRUE)
names(noncomps_output) = c("FocalAOU", "CompetitorAOU", "non_Estimate", "non_P", "non_R2")
noncomps_output_bocc = left_join(noncomps_output, allcomps_output, by = c("FocalAOU", "CompetitorAOU"))
nonps = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(non_R2 >= R2)
names(nonps) = c("FocalAOU", "main_g_non")

none = na.omit(noncomps_output_bocc) %>% 
  group_by(FocalAOU) %>%
  tally(non_Estimate >= Estimate)
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

noncompsdistp_posthoc  = merge(nonps, numcomps, by = ("FocalAOU"))
noncompsdiste_posthoc  = merge(none, numcomps, by = ("FocalAOU"))
noncompsdistp_posthoc$nullp = (noncompsdistp_posthoc$main_g_non)/(noncompsdistp_posthoc$Comp_count + 1)
noncompsdiste_posthoc$nulle = (noncompsdiste_posthoc$main_g_non_e)/(noncompsdiste_posthoc$Comp_count + 1)

noncomps_output_bocc$Null = "Null"
noncomps_output_bocc$Comp = "Comp"
noncompsdistp_posthoc$Null = "Null"
noncompsdiste_posthoc$Null = "Null"

#### Figure 6 example non-comp dist and main R2 ######

p = ggplot(noncompsdistp_posthoc) +
  geom_histogram(bins = 10, aes(nullp), alpha = 0.9, fill="#330066") +
  xlab(expression('Proportion of non-competitors with higher R'^2)) + ylab("Frequency") + theme_classic()  +
  scale_fill_manual(breaks = c("Null"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"))


q = ggplot(noncompsdiste_posthoc) +
  geom_histogram(bins = 10, aes(nulle, fill=factor(Null, levels = c("Null"))), alpha = 0.9) + 
  xlab(expression("Proportion of non-competitors with more negative slope")) + ylab("Frequency") + theme_classic() + 
  scale_fill_manual(breaks = c("Null"), values=c("#330066"), labels=c("Non-Competitors")) + theme(legend.title=element_blank(), legend.text=element_text(size = 12)) + theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24), axis.text.x=element_text(size=24, color = "black"), axis.text.y=element_text(size=24, color = "black")) + theme(plot.margin=unit(c(1,1,1,1),"cm"))

theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(p + theme(legend.position="none"),
          q + theme(legend.position="none"),
          align = 'hv',
          labels = c("a","b", "c", "d"),
          label_size = 20,
          nrow = 2) 

ggsave("Figures/Figure6_null_posthoc.pdf", height = 8, width = 10)

write.csv(noncompsdiste_posthoc, "data/noncompsdiste_posthoc.csv", row.names = FALSE)
write.csv(noncompsdistp_posthoc, "data/noncompsdistp_posthoc.csv", row.names = FALSE)
