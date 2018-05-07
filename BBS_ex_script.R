### This script derives avian temporal occupancy data using ecoretriever to subset raw BBS data
### Also bbs abundance and lat/long data
### sample plot of ranked variance partitioning values

# setwd("C:/Git/Biotic-Interactions/data/Data_for_class")
library(rdataretriever)
library(dplyr)
library(tidyr)
library(ggplot2)

##### ecoretriever to download bbs data and derive occupancy values #####
# bbs_eco = rdataretriever::fetch("breed-bird-survey") # original data source
# read in r data retriever raw bbs data
bbs_raw = read.csv("bbs_eco.csv", header = TRUE)
latlongs = read.csv("latlongs.csv", header = TRUE)
# read in files created in data cleaning script
centroid=read.csv("centroid.csv", header=TRUE) # GIS script
occuenv= read.csv("all_expected_pres.csv", header = TRUE) # Data cleaning script
subsetocc = read.csv('subsetocc.csv', header = T) # Hurlbert Lab
tax_code = read.csv("Tax_AOU_Alpha.csv", header = TRUE) # Hurlbert Lab

# combine state number and route number to get a "stateroute" unique identifier for site
bbs_raw$stateroute = bbs_raw$statenum*1000 + bbs_raw$route

# Get subset of stateroutes that have been surveyed every year from 2001-2015
good_rtes = bbs_raw %>% 
  dplyr::filter(year > 2000, year < 2016) %>% # filtering window of years
  dplyr::select(year, stateroute) %>% # selecting specified columns from dataset
  unique() %>%    # selecting unique rows of dataset
  dplyr::count(stateroute) %>% # counting number of years of data for each stateroute (1-15)
  filter(n == 15) # filtering only routes that have data for all 15 years

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = bbs_raw %>% 
  filter(year > 2000, year < 2016, stateroute %in% good_rtes$stateroute) %>% # filtering to specific years AND stateroutes from good_rtes data frame
  dplyr::select(year, stateroute, aou) %>% # selecting specified columns of dataset (aou = numeric species id for birds)
  unique() %>%    # selecting unique rows of dataset
  dplyr::count(aou, stateroute) # count number of years of data for each stateroute and aou

bbs_sub1$occ = bbs_sub1$n/15 # new column of proportion value for each row (1/15-15/15)

bbs_sub1 = data.frame(bbs_sub1) #converting from tibble to data frame

# Save bbs abundance data
bbs_abun = bbs_raw %>% 
  filter(year > 2000, year < 2016, stateroute %in% good_rtes$stateroute) %>% # filtering to specific years AND stateroutes
  dplyr::select(year, stateroute, aou, speciestotal) # selecting relevant columns

# Save latlong data for bbs routes
latlong_rtes = latlongs %>% 
  dplyr::select(statenum, route, latitude, longitude) %>% # selecting specified columns of dataset
  unique() %>%     # selecting unique rows of dataset
  latlong_rtes$stateroute = latlong_rtes$statenum*1000 + latlong_rtes$route 

# writing final output
# write.csv(latlong_rtes, "data/latlong_rtes.csv", row.names=FALSE)

#### Analysis ####
#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
subsetocc$AOU[subsetocc$AOU == 4810] = 4812
bbs_sub1$aou[bbs_sub1$aou == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4810] = 4812
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

# for loop subsetting env data to expected occurrence for focal species
envoutput = c()
# create beta output data frame
beta_occ = c()

subfocalspecies = unique(occuenv$Species)

for (sp in 1:length(subfocalspecies)){
  print(sp)
  temp = subset(occuenv,Species == subfocalspecies[sp])
  
  competition <- lm(temp$occ_logit ~  temp$comp_scaled)  # changes between main and all comps
  # z scores separated out for env effects (as opposed to multivariate variable)
  env_z = lm(temp$occ_logit ~ abs(zTemp) + abs(zElev) + abs(zPrecip) + abs(zNDVI), data = temp)
  # z scores separated out for env effects
  both_z = lm(temp$occ_logit ~  temp$comp_scaled + abs(temp$zTemp)+abs(temp$zElev)+abs(temp$zPrecip)+abs(temp$zNDVI), data = temp)
  
  # abundance of competitor
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
    
    beta_occ = rbind(beta_occ, c(sp1, occ_comp_est, occ_comp_p, occ_comp_r,occ_env_r, occ_b_p, occ_b_r))
  } 
}         

envoutput = data.frame(envoutput)
names(envoutput) = c("FocalAOU", "ENV", "COMP", "SHARED", "NONE")

#### dplyr merge sequnce ####
envloc = left_join(envoutput, tax_code[,c('AOU_OUT', 'ALPHA.CODE')], by = c('FocalAOU' = "AOU_OUT")) %>%
  left_join(., subsetocc[,c("AOU", "migclass", "Trophic.Group")], by = c('FocalAOU'= "AOU")) %>%
  left_join(., centroid[, c("FocalAOU", "Long", "Lat")], by = 'FocalAOU')
envloc1$EW <- 0
envloc1$EW[envloc1$Long > -98.583333] <- 1 ## mid point of USA, 1 = East
envloc1 = left_join(envloc, occuenv[,c("Species", "Family")], by = c("FocalAOU"="Species"))

# selecting unique subset of data frame
envloc1 = unique(envloc1)
envloc1 = merge(envloc1, tax_code[,c("AOU_OUT", "ALPHA.CODE")], by.x = "FocalAOU", by.y = "AOU_OUT", all.x=TRUE)
envloc1$ALPHA.CODE = as.character(envloc1$ALPHA.CODE.y)
envloc1$ALPHA.CODE[envloc1$FocalAOU == 2920] <- 'MOUQ' #Mountain Quail
envloc1$ALPHA.CODE[envloc1$FocalAOU == 6720] <- 'PAWA' #Palm Warbler
envloc1$ALPHA.CODE = as.factor(envloc1$ALPHA.CODE)

##### Dplyr code ######
nrank = envloc1 %>% 
  dplyr::mutate(rank = row_number(-ENV)) # ordering dataframe by decreasing value of environmental variance
envflip = tidyr::gather(nrank, "Type", "value", 2:5) # creates new dataset in long format
envflip$rank <- factor(envflip$rank, levels = envflip$rank[order(envflip$rank)]) # making rank a factor (for plotting)
envflip = dplyr::arrange(envflip,rank) # arranges new dataframe by rank


envrank = envflip %>% 
  dplyr::group_by(Type == 'ENV') %>% # groups all the entries with type "ENV" 
  dplyr::mutate(rank = row_number(-value)) # ranking data by decreasing value of envrionmental variance
envrank <- envrank[order(envrank$rank),] # arranges new dataframe by rank

envrank <- subset(envrank,Type == "ENV") # change here for comp

###### PLOTTING #####
envflip$Type = factor(envflip$Type,
                      levels = c("NONE", "SHARED","COMP","ENV"),ordered = TRUE)
envflip$value = abs(envflip$value)

# Plot with ENV ranked in decreasing order - had to flip everything to plot right
Fig2b = ggplot(data=envflip, aes(factor(rank), y=value, fill=factor(Type, levels = c("NONE","SHARED","COMP", "ENV")))) + 
  geom_bar(stat = "identity") + theme_classic() +
  theme(axis.text.x=element_text(size=10,vjust=0.5),axis.text.y=element_text(angle=90,size=10)) + xlab("Focal Species") + ylab("Percent Variance Explained") +
  scale_fill_manual(values=c("white","lightskyblue","#dd1c77","#2ca25f"), labels=c("","Shared Variance", "Competition","Environment")) +theme(axis.title.x=element_text(size=40),axis.title.y=element_text(size=30, angle=90),legend.title=element_blank(), legend.text=element_text(size=50, hjust = 1, vjust = 0.5), legend.position = c(0.5,0.9)) # + guides(fill=guide_legend(fill = guide_legend(keywidth = 1, keyheight = 1),title=""))
Fig2b_final = Fig2b + annotate("text", x = 1:104, y = -.03, label = envrank$ALPHA.CODE, angle=90,size=6,vjust=1,hjust = 0.8, color = "black") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_text(size = 40)) + scale_y_continuous(breaks = c(0,0.2,0.4,0.6, 0.8))
Fig2b_final

# save output at specified width
# ggsave("barplot.pdf", height = 25, width = 36)




