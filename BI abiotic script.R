###This script imports environmental data for BBS routes

library(raster)
library(sp)
library(rgdal)
library(RCurl)
library(dplyr)

# Define lat/long window for raster data
box = c(-170,-54,24,74) # North America

# Fill in which box you want to use
my.extent = extent(box) 

# Define projection to be used throughout analysis
prj.string <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"

####################################################################################
#### Import route locations and draw sample circles around them
# setwd("C:/Git/Biotic-Interactions")
# derived from BBS_occ script
routes = read.csv("data/latlong_rtes.csv",header =TRUE)
# routes = routes[row.names(unique(routes[,c('latitude', 'longitude', 'stateroute')])),]
routes$latitude = abs(routes$latitude)
# Makes routes into a spatialPointsDataframe
coordinates(routes)=c('longitude','latitude')
projection(routes) = CRS("+proj=longlat +ellps=WGS84")

# Transforms routes to an equal-area projection - see previously defined prj.string
routes.laea = spTransform(routes, CRS(prj.string))

# A function that draws a circle of radius r around a point: p (x,y)
RADIUS = 40

make.cir = function(p,r){
  points=c()
  for(i in 1:360){
    theta = i*2*pi/360
    y = p[2] + r*cos(theta)
    x = p[1] + r*sin(theta)
    points = rbind(points,c(x,y))
  }
  points=rbind(points,points[1,])
  circle=Polygon(points,hole=F)
  circle
}

#Draw circles around all routes
circs = sapply(1:nrow(routes.laea), function(x){
  circ =  make.cir(routes.laea@coords[x,],RADIUS)
  circ = Polygons(list(circ),ID=routes.laea@data$stateroute[x])
}
)
circs.sp = SpatialPolygons(circs, proj4string=CRS(prj.string))

# Check that circle locations look right
plot(circs.sp)

#####################################################################################
start_yr <- 2001
end_yr <- 2015
min_num_yrs <- 15
richness_w_env <- get_richness_ts_env_data(start_yr, end_yr, min_num_yrs)


setwd('C:/Git/Biotic-Interactions/ENV DATA')
# Read in stack of layers from all 12 months 
tmean <- getData("worldclim", var = "tmean", res = 2.5)
files<-paste('tmean',1:12,'.bil',sep='')
tmeans<-stack(files)
mat = calc(tmeans, mean)

#### ----Precip ----#####
#read in precip data from world clim, stack data to get 1 MAP value
# Read in stack of layers from all 12 months
prec <- getData("worldclim", var = "prec", res = 2.5)
pfiles<-paste('prec',1:12,'.bil',sep='')
pmeans<-stack(pfiles)
map = calc(pmeans, mean)

setwd('C:/Git/Biotic-Interactions')
#### ----Elev ----#####
elev <- raster("Z:/GIS/DEM/sdat_10003_1_20170424_102000103.tif")
NorthAm = readOGR("Z:/GIS/geography", "continent")
NorthAm2 = spTransform(NorthAm, CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))

plot(elevNA2)
plot(NorthAm2)

elevNA2 = projectRaster(elev, crs = prj.string) #UNMASKED!
elevNA3 <- raster::mask(elevNA2, NorthAm2)
# Extract Data
mat.point = raster::extract(mat, routes)
mat.mean = raster::extract(mat, circs.sp, fun = mean, na.rm=T)
mat.var = raster::extract(mat, circs.sp, fun = var, na.rm=T)
#env.data = rbind(env.data, c(env.i,env.point, env.mean,env.var))
#}
env_mat = data.frame(stateroute = names(circs.sp), mat.point = mat.point, mat.mean = mat.mean, mat.var = mat.var)

env_mat = read.csv("data/env_mat.csv", header = TRUE)
# Extract Data
elev.point = raster::extract(elev, routes)
elev.mean = raster::extract(elev, circs.sp, fun = mean, na.rm=T)
elev.var = raster::extract(elev, circs.sp, fun = var, na.rm=T)

env_elev = data.frame(stateroute = names(circs.sp), elev.point = elev.point, elev.mean = elev.mean, elev.var = elev.var)
env_elev=read.csv("data/env_elev.csv", header = TRUE)

# Extract Data
map.point = raster::extract(map, routes)
map.mean = raster::extract(map, circs.sp, fun = mean, na.rm=T)
map.var = raster::extract(map, circs.sp, fun = var, na.rm=T)

env_map = data.frame(stateroute = names(circs.sp), map.point = map.point, map.mean = map.mean, map.var = map.var)
env_map = read.csv("data/env_map.csv", header = TRUE)

# NDVI 
gimms_ndvi = read.csv("ENV DATA/gimms_ndvi_bbs_data.csv", header = TRUE)
gimms_agg = gimms_ndvi %>% filter(month == c("may", "jun", "jul")) %>% 
  group_by(site_id)  %>%  summarise(ndvi.mean=mean(ndvi))
gimms_agg$stateroute = gimms_agg$site_id
ndvi = gimms_agg[,c("stateroute", "ndvi.mean")]


# merge together
all_env = Reduce(function(x, y) merge(x, y, by = "stateroute"), list(env_mat, env_elev, env_map, ndvi))

# write.csv(all_env,'C:/git/Biotic-Interactions/data/all_env.csv',row.names=F)
####----Creating an environmental matrix ----####
occuraw <- read.csv("data/2001_2015_bbs_occupancy.csv", header = T) # read in updated bbs data
occumatrix <- occuraw %>%
  filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

route.locs = read.csv('data/latlong_rtes.csv')
all_env = read.csv("data/all_env.csv", header = TRUE)

latlongs = subset(route.locs, select = c('stateroute', 'latitude', 'longitude'))
latlongs$latitude = abs(latlongs$latitude)
route.sp = coordinates(latlongs[,3:2])

#mean data for all variables
envtable <- subset(all_env, select = c('stateroute', 'mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')) 
envtable$mat.mean = as.numeric(envtable$mat.mean)
envtable$map.mean = as.numeric(envtable$map.mean)
envtable$elev.mean = as.numeric(envtable$elev.mean)
envtable$ndvi.mean = as.numeric(envtable$ndvi.mean)
### Calculate metrics of interest
####---- Creating final data frame----####
#For loop to calculate mean & standard dev environmental variables for each unique species
uniq.spp = unique(occumatrix$aou, header = "Species")

# calc mean spp abundance
bbs_raw = read.csv("data/bbs_abun.csv", header = TRUE)
bbs_abun <- bbs_raw %>%
   filter(aou > 2880) %>%
   filter(aou < 3650 | aou > 3810) %>%
   filter(aou < 3900 | aou > 3910) %>%
   filter(aou < 4160 | aou > 4210) %>%
   filter(aou != 7010)

bbs_abun_sum = bbs_abun %>% 
  group_by(stateroute, aou) %>%
  summarise(abun_weights = sum(speciestotal)) %>%
  left_join(occumatrix, ., by = c("stateroute", "aou")) 
  # left_join(., bbs_abun_sum, by = c("stateroute", "aou"))

occumatrix2 = left_join(bbs_abun_sum, envtable, by = "stateroute")
occumatrix2$tempabun_weights = occumatrix2$abun_weights * occumatrix2$mat.mean
occumatrix2$elevabun_weights = occumatrix2$abun_weights * occumatrix2$elev.mean
occumatrix2$mapabun_weights = occumatrix2$abun_weights * occumatrix2$map.mean
occumatrix2$ndviabun_weights = occumatrix2$abun_weights * occumatrix2$ndvi.mean
# removing routes with NA, only eliminates one focal spp
occumatrix2 = na.omit(occumatrix2)

birdsoutput = data.frame(Mean.Temp=integer(), 
   Mean.Elev = integer(), 
   Mean.Precip = integer(), 
   Mean.NDVI = integer(),
   SD.Temp = integer(), 
   SD.Precip = integer(), 
   SD.NDVI = integer(),
   species = integer(), 
   stringsAsFactors=FALSE) 

for (species in uniq.spp) {
  species = as.numeric(species)
  spec.routes <- occumatrix2[(occumatrix2$aou) == species, "stateroute"] #subset routes for each species (i) in tidybirds
  env.sub <- occumatrix2[occumatrix2$stateroute %in% spec.routes, ] #subset routes for each env in tidybirds
  env.spec <- subset(env.sub, aou == species)
  # calc weighted means
  tempmean = (sum(env.spec$tempabun_weights)/sum(env.spec$abun_weights)) # na.rm = TRUE instead of omit
  mapmean = sum(env.spec$mapabun_weights)/sum(env.spec$abun_weights)
  elevmean = sum(env.spec$elevabun_weights)/sum(env.spec$abun_weights)
  ndvimean = sum(env.spec$ndviabun_weights)/sum(env.spec$abun_weights)
  
  envmeans = cbind(tempmean, mapmean, elevmean, ndvimean)
  envmeans = data.frame(envmeans)

  envar = c()
  for(i in 1:nrow(envmeans)){
  # calc weighted sd
    sd_num = sum(env.spec$abun_sum * (env.spec$mat.mean - tempmean)^2) # write a function=, input = env var, delete NAs
    sd_denom = (sum(env.spec$abun_sum > 0) -1) * sum(env.spec$abun_sum) / sum(env.spec$abun_sum > 0) 
    sd = sqrt(sd_num/sd_denom)
    
  tempvar = sqrt(abs(sum(na.omit(env.spec$mat.mean - tempmean)*env.spec$abun_sum)^2/(length(env.spec$abun_mean)-1*sum(na.omit(env.spec$tempabun_sum)))/length(env.spec$abun_mean)))
 #wtd.var(env.spec$tempabun_sum/env.spec$abun_mean))
  mapvar = sqrt(abs((sum(na.omit(envmeans[i, 2])))^2/(length(env.spec$abun_mean)-1*sum(na.omit(env.spec$mapabun_sum)))/length(env.spec$abun_mean)))
  elevvar = sqrt(abs((sum(na.omit(envmeans[i, 3])))^2/(length(env.spec$abun_mean)-1*sum(na.omit(env.spec$elevabun_sum)))/length(env.spec$abun_mean)))
  ndvivar = sqrt(abs((sum(na.omit(envmeans[i, 4])))^2/(length(env.spec$abun_mean)-1*sum(na.omit(env.spec$ndviabun_sum)))/length(env.spec$abun_mean)))
  envar = rbind(envar, c(tempvar, mapvar, elevvar, ndvivar))
  }
  envar = data.frame(envar)
  names(envar) = c("SD.Temp", "SD.Precip", "SD.Elev", "SD.NDVI")

  envmeans$species = species
  
  birdsoutput = rbind(birdsoutput, c(envmeans, envar))
}

# names(birdsoutput) = c("Species", "Mean.Temp", "Mean.Precip", "Mean.Elev", "Mean.NDVI", "SD.Temp", "SD.Precip", "SD.Elev", "SD.NDVI")
### Combine relevant information from each of your two or more datasets using merge()
#(species/occupancy/expected env variables/observed env variables)
occubirds <- left_join(birdsoutput, occumatrix[c("aou", "stateroute","occ")], by = c("species" = "aou"), na.rm = T)
occuenv <- left_join(occubirds, envtable, by = "stateroute", na.rm = T)
# occuenv <- na.omit(occuenv) got rid of too many aous

### Conduct analyses
#Calculating z scores for each environmnetal variable (observed mean - predicted mean/predicted SD)
occuenv$zTemp = (occuenv$mat.mean - occuenv$Mean.Temp) / occuenv$SD.Temp
occuenv$zPrecip = (occuenv$map.mean - occuenv$Mean.Precip) / occuenv$SD.Precip
occuenv$zElev = (occuenv$elev.mean - occuenv$Mean.Elev) / occuenv$SD.Elev
occuenv$zNDVI = (occuenv$ndvi.mean - occuenv$Mean.NDVI) / occuenv$SD.NDVI

write.csv(occuenv, "data/occuenv.csv", row.names= FALSE)






