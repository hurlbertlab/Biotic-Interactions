### This script derives avian temporal occupancy data using ecoretriever to subset raw BBS data
### Also bbs abundance and lat/long data

# setwd("C:/git/Biotic-Interactions")
library(rdataretriever)
##### ecoretriever to download bbs data and derive occupancy values #####
bbs_eco = rdataretriever::fetch("breed-bird-survey") # takes forever!!

Years = (bbs_eco$breed_bird_survey_counts)
# bbs_eco$breed_bird_survey_counts = as.numeric(bbs_eco$breed_bird_survey_counts$Year)
Years$stateroute = Years$statenum*1000 + Years$route

# Get subset of stateroutes that have been surveyed every year from 2001-2015
good_rtes = Years %>% 
  dplyr::filter(year > 2000, year < 2016) %>% 
  dplyr::select(year, stateroute) %>%
  unique() %>%    
  dplyr::count(stateroute) %>% 
  filter(n == 15) # have to stay at 15 to keep # of years consistent

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = Years %>% 
  filter(year > 2000, year < 2016, stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(year, stateroute, aou) %>%
  dplyr::count(aou, stateroute) 

bbs_sub1$occ = bbs_sub1$n/15 # new occupancy values calculated
write.csv(bbs_sub1, "2001_2015_bbs_occupancy.csv", row.names=FALSE)

# Save bbs abundance data
bbs_abun = Years %>% 
  filter(year > 2000, year < 2016, stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(year, stateroute, aou, speciestotal)
write.csv(bbs_abun, "bbs_abun.csv", row.names=FALSE)

# Save latlong data for bbs routes
latlong_rtes = bbs_eco$breed_bird_survey_routes %>% 
  dplyr::select(statenum, route, latitude, longitude) %>%
  unique() %>%    
  group_by(statenum, route) 
latlong_rtes$stateroute = latlong_rtes$statenum*1000 + latlong_rtes$route 
write.csv(latlong_rtes, "latlong_rtes.csv", row.names=FALSE)


# get occupancy by stop
read.csv("Z:/Gartland/BBS scaled/output.csv", header= TRUE)
bbs_abun = merge(bbs_eco, output, by = c("stateroute", "scale", "Aou"))


