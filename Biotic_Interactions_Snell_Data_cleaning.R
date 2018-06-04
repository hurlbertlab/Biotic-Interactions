# Biotic Interactions script
# In this script we compare the differences in occupancy and abundance between
# focal species and their competitors using occupancy, abundance, and environmental data.
# Env data was formatted in Snell_code.R from BIOL 465 project. Occupancy data from BBS ecoretriever.

#setwd("C:/git/Biotic-Interactions")
#### ---- Inital Formatting ---- ####
library(dplyr)
library(tidyr)

# read in all dataset inputs
bbs = read.csv('data/bbs_abun.csv', header = T) # BBS abundance data - from Hurlbert Lab 
expect_pres = read.csv('data/expect_pres.csv', header = T) # expected presence data based on BBS for 372 landbird species --- 2001-2015 from NatureServ
bbs_occ = read.csv("data/bbs_sub1.csv", header=TRUE) # BBS temporal occupancy data (just for landbird species) --- 2001-2015
bsize = read.csv("data/DunningBodySize_old_2008.11.12.csv", header = TRUE) # import body size data from Dunning 2008
focal_competitor_table = read.csv("data/focal spp.csv", header = TRUE)
AOU = read.csv("data/Bird_Taxonomy.csv", header = TRUE) # taxonomy data
shapefile_areas = read.csv("data/shapefile_areas.csv", header = TRUE) # area shapefile if not running GIS code 
subsetocc = read.csv("data/subsetocc.csv", header = TRUE)
# read in raw env data UPDATED from gimms script
all_env = read.csv('data/occuenv.csv', header = T)
# in 2016 Western Scrub jay split into CA and Island scrub jay, not reflected in most dataframes
expect_pres$spAOU[expect_pres$spAOU == 4810] = 4812
bbs_occ$Aou[bbs_occ$Aou == 4810] = 4812
bsize$AOU[bsize$AOU == 4810] = 4812
AOU$AOU_OUT[AOU$AOU_OUT == 4810] = 4812
shapefile_areas$compAOU[shapefile_areas$compAOU == 4810] = 4812
subsetocc$AOU[subsetocc$AOU == 4810] = 4812
all_env$Species[all_env$Species == 4810] = 4812

# subset temporal occupancy
temp_occ = subset(bbs_occ, Aou %in% subsetocc$AOU)
############# ---- Set up pairwise comparison table ---- #############
# read in table with pairwise comparison of each focal species to several potential competitors - created by hand
focal_competitor_table = dplyr::select(focal_competitor_table, AOU, CommonName, Competitor)
names(focal_competitor_table)= c("FocalAOU", "Focal", "Competitor")

# create data frame of unique focal species
focal_unique = data.frame(unique(focal_competitor_table[, c("Focal", "FocalAOU")]))
names(focal_unique) = c("Focal_Common","FocalAOU")

AOU=AOU[, (names(AOU) %in% c("SCI_NAME", "AOU_OUT", "PRIMARY_COM_NAME","FAMILY"))]
names(AOU) = c("SciName","CommonName",  "AOU", "Family")

# remove duplicates/subspecies from AOU data frame
AOUsub = AOU[-grep("sp\\.", AOU$CommonName),] 
AOUsub2 = AOUsub[-grep("\\)", AOUsub$CommonName),]
AOUsub3 = AOUsub2[-grep(" \\(", AOUsub2$CommonName),]
AOUsub4 = AOUsub3[-grep("/", AOUsub3$CommonName),] 
AOUsub5 = unique(AOUsub4)
sp_list = na.omit(AOUsub5)

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

sp_list$match[sp_list$match =="Calcarius mccownii"] = "Rhynchophanes mccownii"

sp_list$match[sp_list$match =="Picoides villosus"] = "Leuconotopicus villosus"

sp_list$match[sp_list$match =="Picoides pubescens"] = "Dryobates pubescens"

sp_list$match[sp_list$match =="Picoides dorsalis"] = "Picoides tridactylus"

sp_list$match[sp_list$match =="Picoides scalaris"] = "Dryobates scalaris"

sp_list$match[sp_list$match =="Picoides albolarvatus"] = "Leuconotopicus albolarvatus"

sp_list$match[sp_list$match =="Picoides borealis"] = "Leuconotopicus borealis"

sp_list$match[sp_list$match =="Carduelis hornemanni"] = "Carduelis flammea"

sp_list$match[sp_list$match =="Aimophila cassinii"] = "Peucaea cassinii"

sp_list$match[sp_list$match =="Aimophila aestivalis"] = "Peucaea aestivalis"

sp_list$match[sp_list$match =="Aimophila botterii"] = "Peucaea botterii"

sp_list$match[sp_list$match =="Aimophila carpalis"] = "Peucaea carpalis"

sp_list$match[sp_list$match =="Oreothlypis crissalis"] = "Vermivora crissalis"

sp_list$match[sp_list$match =="Ixoreus naevius"] = "Zoothera naevia"
###### ---- Create final focal-comp table ----######
#merge pairwise table with taxonomy info
comp_AOU = merge(focal_competitor_table, sp_list, by.x = "Competitor", by.y = "CommonName")
names(comp_AOU) = c("Competitor", "focalAOU", "Focal", "old", "CompAOU", "CompFamily","CompSciName")
comp_AOU = dplyr::select(comp_AOU, -c(old, CompFamily))
comp_AOU$CompAOU[comp_AOU$CompAOU == 5660] = 5677
# merging in focal sci name to table
focal_AOU = merge(comp_AOU, sp_list[,c("match", "AOU", "Family")], by.x = "focalAOU", by.y = "AOU")
names(focal_AOU)[6] = "FocalSciName"
focal_AOU$focalAOU[focal_AOU$focalAOU == 5660] = 5677
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
write.csv(new_spec_weights, "data/new_spec_weights.csv", row.names=FALSE) 
############# ---- Generate total species occupancies ---- #############
# Take range overlap area to assign "main competitor" for each focal species
# "area.df" with cols: FocalAOU, CompAOU, focalArea, compArea, intArea, intProp
shapefile_areas = na.omit(shapefile_areas)
shapefile_areas$focalAOU[shapefile_areas$focalAOU == 5660] = 5677
shapefile_areas$compAOU[shapefile_areas$compAOU == 5660] = 5677

# calculate proportion of overlap between focal range and overlap range
shapefile_areas$PropOverlap = shapefile_areas$area_overlap/shapefile_areas$FocalArea

# Which competitor has greatest area of overlap? -- main competitor
shapefile_areas$mainCompetitor = 0 # set up main competitor column, 0 = not the primary competitor
for (s in unique(shapefile_areas$focalAOU)) {
  maxOverlap = max(shapefile_areas$PropOverlap[shapefile_areas$focalAOU == s], na.rm = TRUE) #largest area of proportion overlap
  shapefile_areas$mainCompetitor[shapefile_areas$focalAOU == s & shapefile_areas$PropOverlap == maxOverlap] = 1 # 1 assigns main competitor
}

#### ---- Gathering Occupancy and Abundance Data for Biotic Comparisons ---- ####
# pull out stateroutes that have been continuously sampled 2001-2015
routes = unique(temp_occ$stateroute)

focal_and_comp_species = unique(c(new_spec_weights$focalAOU, new_spec_weights$CompAOU))

# pooling BBS mean abundance by AOU/stateroute and by year  
bbs_pool = bbs %>% 
  group_by(stateroute, aou) %>% 
  dplyr::summarize(abundance = mean(speciestotal)) %>%
  dplyr::filter(aou %in% focal_and_comp_species) 
names(bbs_pool)[names(bbs_pool)=="aou"] <- "AOU"
bbs_pool = data.frame(bbs_pool)

# merge bbs pooled abundances with expected presences for all species
bbs_ep = full_join(expect_pres[,c("stateroute", "spAOU")], bbs_pool, 
                   by = c('stateroute' = 'stateroute', 'spAOU' = 'AOU')) %>%
  filter(stateroute %in% routes)
bbs_ep$abundance[is.na(bbs_ep$abundance)] = 0
bbs_ep$spAOU[bbs_ep$spAOU == 5660] = 5677

# Calculating summed competitor abundance for each focal species
prefull_data = left_join(focal_AOU, bbs_ep, by = c('focalAOU'= 'spAOU')) %>% 
  left_join(temp_occ[temp_occ$Aou %in% focal_and_comp_species, ], 
            by = c('focalAOU' = 'Aou', 'stateroute' = 'stateroute')) %>%
  left_join(bbs_pool, by = c('CompAOU' = 'AOU', 'stateroute' = 'stateroute')) %>%
  group_by(focalAOU, Focal, Family, stateroute, occ, abundance.x) %>% 
  dplyr::summarize(allCompN = sum(abundance.y, na.rm = T)) %>% na.omit(.)

# subsetting prefull data to species in new_spec_weights
prefull_data2 = subset(prefull_data, focalAOU %in% new_spec_weights$focalAOU)
prefull_data2 = data.frame(prefull_data2)
# create focalcompoutput table that adds MainCompN column to indicate primary competitors to the 
# prefull_data with focal/comp/stroute/abundance/occ/summed abundance
focalcompoutput.5 = bbs_ep %>%    
  dplyr::select(stateroute, spAOU) %>%
  left_join(subset(shapefile_areas, mainCompetitor == 1, 
                   select = c('focalAOU', 'compAOU', 'mainCompetitor')), 
            by = c('spAOU' = 'focalAOU')) %>%
  left_join(bbs_pool, by = c('stateroute' = 'stateroute', 'compAOU' = 'AOU')) %>%
  left_join(prefull_data, by = c('spAOU' = 'focalAOU', 'stateroute' = 'stateroute')) %>%
  dplyr::select(stateroute, Focal, spAOU, Family, abundance.x, occ, abundance, allCompN)
names(focalcompoutput.5) = c("stateroute","Focal", "FocalAOU", "Family", "FocalAbundance", "FocalOcc","MainCompN", "AllCompN")

focalcompoutput = focalcompoutput.5[!is.na(focalcompoutput.5$FocalOcc)  & !is.na(focalcompoutput.5$AllCompN),] 
focalcompoutput$FocalOcc[is.na(focalcompoutput$FocalOcc)] = 0
focalcompoutput$MainCompN[is.na(focalcompoutput$MainCompN)] = 0


#### selecting competitors for noncomp analysis ####
uniq_comps = unique(shapefile_areas$compAOU)
uniq_foc = filter(shapefile_areas, !focalAOU %in% uniq_comps)[,2]
uniq_spp = data.frame(c(uniq_comps, uniq_foc))
noncomps = left_join(uniq_spp, AOU[c("AOU", "Family")], by = c("c.uniq_comps..uniq_foc." = "AOU"))
noncomps = noncomps[!duplicated(noncomps),]
noncomps$AOU = noncomps$c.uniq_comps..uniq_foc.
# write.csv(noncomps[,c(2,3)], "data/noncomps.csv", row.names = FALSE)



# Filter number to spp present at 40+ routes for better model results
# Subset to get the count of routes for each spp
sppGT50rtes = focalcompoutput %>%
  group_by(FocalAOU, Focal) %>%
  summarise(n = n_distinct(stateroute)) %>%
  filter(n>=50) %>% 
  dplyr::select(FocalAOU)

# Merge with focalcompoutput data table, new # of focal spp is 171 with route filters applied
focalcompsub = filter(focalcompoutput, FocalAOU %in% sppGT50rtes$FocalAOU)

# Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) ### FOR MAIN
focalcompsub$comp_scaled = focalcompsub$MainCompN/(focalcompsub$FocalAbundance + focalcompsub$MainCompN)

# Create scaled competitor column = main comp abundance/(focal abundance + main comp abundance) ### FOR ALL
focalcompsub$all_comp_scaled = focalcompsub$AllCompN/(focalcompsub$FocalAbundance + focalcompsub$AllCompN)


#### ---- Processing Environmental Data - Re-done from Snell_abiotic_code.R ---- ####
# merge in ENV

all_expected_pres = left_join(focalcompsub, all_env, by = c("stateroute" = "stateroute", "FocalAOU" = "Species"))
all_expected_pres = na.omit(all_expected_pres)
all_expected_pres = all_expected_pres[all_expected_pres$Focal != "Rock Pigeon",]

write.csv(all_expected_pres,"data/all_expected_pres.csv", row.names= F)

####### END DATA CLEANING, see analysis script ##########