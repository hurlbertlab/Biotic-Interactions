#setwd("C:/git/Biotic-Interactions")
library(dplyr)
library(tidyr)

# read in all dataset inputs
# dataset from Hurlbert and White 2007
Hurlbert_o = read.csv('data/Master_RO_Correlates_20110610.csv', header = T) # range occupancy dataset 

Hurlbert_sub = Hurlbert_o %>% dplyr::select(AOU, CommonName, Route.Expected, X10yr.Prop, migclass, Foraging, Trophic.Group, log_Brange_Area)
# subset species whose range occupancies were greater than 0.3 over a 10 year period
subsetocc = Hurlbert_sub[Hurlbert_sub$X10yr.Prop > .3,]
# Winter Wren had AOU code change (7220 to 7222), changing in occ code to reflect that
subsetocc$AOU[subsetocc$AOU == 7220] <- 7222
# updating dark-eyed junco AOU
subsetocc$AOU[subsetocc$AOU == 5660] <- 5677
# rbind Pacific wren to data frame
pacific = data.frame(7221, "Pacific Wren", 280, 0.589,"short", "ground glean", "insectivore", 6.711282)
colnames(pacific) = c("AOU", "CommonName", "Route.Expected", "X10yr.Prop","migclass","Foraging","Trophic.Group", "log_Brange_Area")
subset_occ = rbind(subsetocc, pacific)

write.csv(subset_occ, "data/subsetocc.csv", row.names = FALSE)


comps_output = c()
for (sp in family$AOU){
  FocalAOU = sp
  temp = subset(family, AOU == sp) 
  tempfam = unique(as.character(temp$FAMILY))
  ncomps = dplyr::filter(allspp, FAMILY == tempfam) %>%
      dplyr::select(AOU_OUT, FAMILY, PRIMARY_COM_NAME)
  ncomps = na.omit(ncomps)
  ncomps = cbind(ncomps, FocalAOU)
  comps_output = rbind(comps_output, ncomps)
  }

comps_output_focal = merge(comps_output, AOU[,c("AOU_OUT", "PRIMARY_COM_NAME")], by.x = "FocalAOU", by.y = "AOU_OUT")
comps_output_focal = unique(comps_output_focal)
write.csv(comps_output_focal, "data/comps_output_focal.csv", row.names = FALSE)
