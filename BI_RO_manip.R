#setwd("C:/git/Biotic-Interactions")
library(dplyr)
library(tidyr)

# read in all dataset inputs
# dataset from Hurlbert and White 2007
Hurlbert_o = read.csv('data/Master_RO_Correlates_20110610.csv', header = T) # range occupancy dataset 

Hurlbert_sub = Hurlbert_o %>% dplyr::select(AOU, CommonName, Route.Expected, X10yr.Prop, migclass, Foraging, Trophic.Group)
# subset species whose range occupancies were between 0.3 and 0.7 over a 10 year period
subsetocc = Hurlbert_sub[Hurlbert_sub$X10yr.Prop > .3 & Hurlbert_sub$X10yr.Prop < .7,]
# Winter Wren had AOU code change (7220 to 7222), changing in occ code to reflect that
subsetocc$AOU[subsetocc$AOU == 7220] <- 7222
# updating dark-eyed junco AOU
subsetocc$AOU[subsetocc$AOU == 5660] <- 5677
# rbind Pacific wren to data frame
pacific = data.frame(7221, "Pacific Wren", 280, 0.589,"short", "ground glean", "insectivore")
colnames(pacific) = c("AOU", "CommonName", "Route.Expected", "X10yr.Prop","migclass","Foraging","Trophic.Group")
subset_occ = rbind(subsetocc, pacific)

write.csv(subset_occ, "data/subsetocc.csv", row.names = FALSE)
