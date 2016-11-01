### This script derives avian temporal occupancy data using ecoretriever to subset raw BBS data
### Also bbs abundance and lat/long data

# setwd("C:/git/Biotic-Interactions")

##### ecoretriever to download bbs data and derive occupancy values #####
bbs_eco = ecoretriever::fetch("BBS") # takes forever!!

Years = (bbs_eco$counts$Year)
bbs_eco$counts$Year = as.numeric(bbs_eco$counts$Year)
bbs_eco$counts$stateroute = bbs_eco$counts$statenum*1000 + bbs_eco$counts$Route

# Get subset of stateroutes that have been surveyed every year from 1996-2010
good_rtes = bbs_eco$counts %>% 
  filter(Year >= 2001, Year <= 2015) %>% 
  dplyr::select(Year, stateroute) %>%
  unique() %>%    
  group_by(Year) %>% 
  dplyr::count(stateroute) %>% 
  filter(n == 15) # have to stay at 15 to keep # of years consistent

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = bbs_eco$counts %>% 
  filter(Year > 2000, Year < 2016, stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(Year, stateroute, Aou) %>%
  unique() %>%
  dplyr::count(Aou, stateroute) 

bbs_sub1$occ = bbs_sub1$n/15 # new occupancy values calculated
write.csv(bbs_sub1, "bbs_sub1.csv", row.names=FALSE)

# Save bbs abundance data
bbs_abun = bbs_eco$counts %>% 
  filter(Year > 2000, Year < 2016, stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(Year, stateroute, Aou, SpeciesTotal)
write.csv(bbs_abun, "bbs_abun.csv", row.names=FALSE)

# Save latlong data for bbs routes
latlong_rtes = bbs_eco$routes %>% 
  dplyr::select(statenum, route, latitude, longitude) %>%
  unique() %>%    
  group_by(statenum, route) 
latlong_rtes$stateroute = latlong_rtes$statenum*1000 + latlong_rtes$route 
write.csv(latlong_rtes, "latlong_rtes.csv", row.names=FALSE)