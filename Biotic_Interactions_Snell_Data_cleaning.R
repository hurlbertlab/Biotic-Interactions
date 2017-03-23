# Biotic Interactions script
# In this script we compare the differences in occupancy and abundance between
<<<<<<< HEAD
# focal species and their competitors using occupancy, abundance, and environmental data.
=======
# 63 focal species and their competitors using occupancy, abundance, and environmental data.
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
# Env data was formatted in Snell_code.R from BIOL 465 project. Occupancy data from BBS ecoretriever.

#setwd("C:/git/Biotic-Interactions")
#### ---- Inital Formatting ---- ####
library(dplyr)
library(tidyr)

# read in range occupancy dataset 
Hurlbert_o = read.csv('Master_RO_Correlates_20110610.csv', header = T)
# subset species whose range occupancies were between 0.3 and 0.7 over a 10 year period
subsetocc = Hurlbert_o[Hurlbert_o$X10yr.Prop > .3 & Hurlbert_o$X10yr.Prop < .7,]

# read in BBS abundance data - from Hurlbert Lab 
# bbs = read.csv('dataset_1.csv', header = T)
bbs = read.csv('bbs_abun.csv', header = T)

# subset bbs abundance columns
bbs = bbs[, (names(bbs) %in% c("stateroute", "Aou", "Year","SpeciesTotal",  'routeID', 'Lati', 'Longi'))]

# read in expected presence data based on BBS for 372 landbird species --- 2001-2015
expect_pres = read.csv('expect_pres.csv', header = T)

# read in BBS temporal occupancy data (just for 372 landbird species) --- 2001-2015
temp_occ = read.csv("bbs_sub1.csv", header=TRUE) 
temp_occ = subset(temp_occ, Aou %in% subsetocc$AOU)


############# ---- Set up pairwise comparison table ---- #############
# read in table with pairwise comparison of each focal species to several potential competitors - created by hand
focal_competitor_table = read.csv("focal spp.csv", header = TRUE)
focal_competitor_table = dplyr::select(focal_competitor_table, AOU, CommonName, Competitor)
names(focal_competitor_table)= c("FocalAOU", "Focal", "Competitor")

# create data frame of unique focal species
focal_unique = data.frame(unique(focal_competitor_table[, c("Focal", "FocalAOU")]))
names(focal_unique) = c("Focal_Common","FocalAOU")

# read in taxonomy data
AOU = read.csv("Bird_Taxonomy.csv", header = TRUE)
AOU=AOU[, (names(AOU) %in% c("SCI_NAME", "AOU_OUT", "PRIMARY_COM_NAME","FAMILY"))]
names(AOU) = c("SciName","CommonName",  "AOU", "Family")

# remove duplicates/subspecies from AOU data frame
AOUsub = AOU[-grep("sp\\.", AOU$CommonName),] 
AOUsub2 = AOUsub[-grep("\\)", AOUsub$CommonName),]
AOUsub3 = AOUsub2[-grep(" \\(", AOUsub2$CommonName),]
AOUsub4 = unique(AOUsub3)
sp_list = na.omit(AOUsub4)

############# ----  nomenclature corrections for shapefiles,correct name is "match" column ---- ######
sp_list$match = as.character(sp_list$SciName)
# renaming to get latest scientific names for mismatch spp
sp_list$match[sp_list$match =="Oreothlypis peregrina"] = "Vermivora peregrina"

sp_list$match[sp_list$match =="Vermivora pinus"] = "Vermivora cyanoptera"

sp_list$match[sp_list$match =="Stellula calliope"] = "Selasphorus calliope"

sp_list$match = gsub('Setophaga ', 'Dendroica ', sp_list$match)

sp_list$match[sp_list$match =="Dendroica ruticilla"] = "Setophaga ruticilla"

sp_list$match[sp_list$match =="Picoides nuttallii"] = "Dryobates nuttallii"

sp_list$match[sp_list$match =="Cardellina canadensis"] = "Wilsonia canadensis"

sp_list$match[sp_list$match =="Geothlypis philadelphia"] = "Oporornis philadelphia"

sp_list$match[sp_list$match =="Oreothlypis ruficapilla"] = "Vermivora ruficapilla"

sp_list$match[sp_list$match =="Oreothlypis celata"] = "Vermivora celata"

sp_list$match[sp_list$match =="Cardellina pusilla"] = "Wilsonia pusilla"

sp_list$match[sp_list$match =="Oreothlypis virginiae"] = "Vermivora virginiae"

sp_list$match[sp_list$match =="Poecile hudsonica"] = "Parus hudsonicus"

sp_list$match[sp_list$match =="Pica hudsonia"] = "Pica pica"

sp_list$match = gsub('Poecile ', 'Parus ', sp_list$match)

sp_list$match[sp_list$match =="Dendroica citrina"] = "Wilsonia citrina"

sp_list$match[sp_list$match =="Geothlypis formosus"] = "Oporornis formosus"

sp_list$match[sp_list$match =="Oreothlypis luciae"] = "Vermivora luciae"

sp_list$match[sp_list$match =="Geothlypis tolmiei"] = "Oporornis tolmiei"

sp_list$match[sp_list$match =="Troglodytes hiemalis"] = "Troglodytes troglodytes"

# Winter Wren had AOU code change (7220 to 7222), changing in occ code to reflect that
temp_occ$Aou[temp_occ$Aou == 7220] <- 7222

###### ---- Create final focal-comp table ----######
#merge pairwise table with taxonomy info
comp_AOU = merge(focal_competitor_table, sp_list, by.x = "Competitor", by.y = "CommonName")
names(comp_AOU) = c("Competitor", "focalAOU", "Focal", "old", "CompAOU", "CompFamily","CompSciName")
comp_AOU = dplyr::select(comp_AOU, -c(old, CompFamily))

# merging in focal sci name to table
focal_AOU = merge(comp_AOU, sp_list[,c("match", "AOU", "Family")], by.x = "focalAOU", by.y = "AOU")
names(focal_AOU)[6] = "FocalSciName"

# import body size data from Dunning 2008
bsize = read.csv("DunningBodySize_old_2008.11.12.csv", header = TRUE)
bsize$AOU[bsize$AOU == 7220] <- 7222 # Winter Wren

# merge in competitor and focal body size
spec_w_bsize = merge(focal_AOU, bsize[,c("AOU", "Mass.g.")], by.x = "focalAOU", by.y = "AOU")
spec_w_weights = merge(spec_w_bsize, bsize[,c("AOU", "Mass.g.")], by.x = "CompAOU", by.y = "AOU")

names(spec_w_weights)[8] = "FocalMass"
names(spec_w_weights)[9] = "CompMass"

# want to compare body size - if competitor is double or more in size to focal, then delete
new_spec_weights = subset(spec_w_weights, spec_w_weights$FocalMass / spec_w_weights$CompMass >= 0.5 &
                            spec_w_weights$FocalMass / spec_w_weights$CompMass <= 2)

# adding in underscore for file name matching
new_spec_weights$focalcat = gsub(" ", "_", new_spec_weights$FocalSciName)
new_spec_weights$compcat = gsub(" ", "_", new_spec_weights$CompSciName)

#write this data frame for GIS script
write.csv(new_spec_weights, "new_spec_weights.csv", row.names=FALSE) 
############# ---- Generate total species occupancies ---- #############
# Take range overlap area to assign "main competitor" for each focal species
# "area.df" with cols: FocalAOU, CompAOU, focalArea, compArea, intArea, intProp
# read in area shapefile if not running GIS code 
shapefile_areas = read.csv("shapefile_areas.csv", header = TRUE)
shapefile_areas = dplyr::select(shapefile_areas, -X)

# calculate proportion of overlap between focal range and overlap range
shapefile_areas$PropOverlap = shapefile_areas$area_overlap/shapefile_areas$FocalArea

# Which competitor has greatest area of overlap? -- main competitor
shapefile_areas$mainCompetitor = 0 # set up main competitor column, 0 = not the primary competitor
for (s in unique(shapefile_areas$focalAOU)) {
  maxOverlap = max(shapefile_areas$PropOverlap[shapefile_areas$focalAOU == s], na.rm = TRUE) #largest area of proportion overlap
  shapefile_areas$mainCompetitor[shapefile_areas$focalAOU == s & shapefile_areas$PropOverlap == maxOverlap] = 1 # 1 assigns main competitor
}

# pull out stateroutes that have been continuously sampled 2001-2015
routes = unique(temp_occ$stateroute)

#### ---- Gathering Occupancy and Abundance Data for Biotic Comparisons ---- ####
focal_and_comp_species = unique(c(new_spec_weights$focalAOU, new_spec_weights$CompAOU))

# need to change winter wren AOU to 7222 from 7220 in bbs_pool
bbs$Aou[bbs$Aou == 7220] <- 7222

# pooling BBS mean abundance by AOU/stateroute and by year  
bbs_pool = bbs %>% 
  group_by(stateroute, Aou) %>% 
  dplyr::summarize(abundance = mean(SpeciesTotal)) %>%
  filter(Aou %in% focal_and_comp_species) 
names(bbs_pool)[names(bbs_pool)=="Aou"] <- "AOU"
bbs_pool = data.frame(bbs_pool)

# merge bbs pooled abundances with expected presences for all species
bbs_ep = full_join(expect_pres[,c("stateroute", "spAOU")], bbs_pool, by = c('stateroute' = 'stateroute', 'spAOU' = 'AOU')) %>%
  filter(stateroute %in% routes)
bbs_ep$abundance[is.na(bbs_ep$abundance)] = 0

# Calculating summed competitor abundance for each focal species
prefull_data = left_join(bbs_ep, focal_AOU, by = c('spAOU' = 'focalAOU')) %>% 
  left_join(temp_occ[temp_occ$Aou %in% focal_and_comp_species, ], 
            by = c('spAOU' = 'Aou', 'stateroute' = 'stateroute')) %>%
  left_join(bbs_pool, by = c('CompAOU' = 'AOU', 'stateroute' = 'stateroute')) %>%
  group_by(spAOU, Focal, Family, stateroute, occ, abundance.x) %>%
  dplyr::summarize(allCompN = sum(abundance.y, na.rm = T))

# create focalcompoutput table that adds MainCompN column to indicate primary competitors to the 
# prefull_data with focal/comp/stroute/abundance/occ/summed abundance
focalcompoutput = bbs_ep %>%    
  dplyr::select(stateroute, spAOU) %>%
  left_join(subset(shapefile_areas, mainCompetitor == 1, 
                   select = c('focalAOU', 'compAOU', 'mainCompetitor')), 
            by = c('spAOU' = 'focalAOU')) %>%
  left_join(bbs_pool, by = c('stateroute' = 'stateroute', 'compAOU' = 'AOU')) %>%
  left_join(prefull_data, by = c('spAOU' = 'spAOU', 'stateroute' = 'stateroute')) %>%
  dplyr::select(stateroute, Focal, spAOU, Family, abundance.x, occ, abundance, allCompN)
names(focalcompoutput) = c("stateroute","Focal", "FocalAOU", "Family", "FocalAbundance", "FocalOcc","MainCompN", "AllCompN")

focalcompoutput$FocalOcc[is.na(focalcompoutput$FocalOcc)] = 0
focalcompoutput$MainCompN[is.na(focalcompoutput$MainCompN)] = 0
focalcompoutput = focalcompoutput[focalcompoutput$FocalOcc & focalcompoutput$AllCompN > 0,] 

# Filter number to spp present at 20+ routes for better model results
# Subset to get the count of routes for each spp
sppGT20rtes = focalcompoutput %>%
  group_by(FocalAOU) %>%
  summarise(n = n_distinct(stateroute)) %>%
  filter(n>=20) %>%
  dplyr::select(FocalAOU)

# Merge with focalcompoutput data table, new # of focal spp is 171 with route filters applied
focalcompsub = filter(focalcompoutput, FocalAOU %in% sppGT20rtes$FocalAOU)

# Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) ### FOR MAIN
focalcompsub$comp_scaled = focalcompsub$MainCompN/(focalcompsub$FocalAbundance + focalcompsub$MainCompN)

# Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) ### FOR ALL
focalcompsub$all_comp_scaled = focalcompsub$AllCompN/(focalcompsub$FocalAbundance + focalcompsub$AllCompN)


#### ---- Processing Environmental Data - Re-done from Snell_abiotic_code.R ---- ####
<<<<<<< HEAD
# read in raw env data UPDATED from gimms script
all_env = read.csv('Z:/Snell/occuenv.csv', header = T)
=======
# read in raw env data UPDATED from MODIS script
all_env = read.csv('occuenv.csv', header = T)
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
# merge in ENV
all_expected_pres = merge(all_env, focalcompsub, by.x = c("stateroute", "Species"), by.y = c("stateroute", "FocalAOU"))

write.csv(all_expected_pres,"all_expected_pres.csv", row.names= F)

####### END DATA CLEANING, see analysis script ##########