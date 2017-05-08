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

# write.csv(all_env,'C:/git/Biotic-Interactions/all_env.csv',row.names=F)
####----Creating an environmental matrix ----####
occumatrix <- read.csv("data/2001_2015_bbs_occupancy.csv", header = T) # read in updated bbs data
route.locs = read.csv('data/latlong_rtes.csv')

latlongs = subset(route.locs, select = c('stateroute', 'latitude', 'longitude'))
latlongs$latitude = abs(latlongs$latitude)
route.sp = coordinates(latlongs[,3:2])

#mean data for all variables
envtable <- subset(all_env, select = c('stateroute', 'mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')) 

### Calculate metrics of interest
####---- Creating final data frame----####
#For loop to calculate mean & standard dev environmental variables for each unique species
uniq.spp = unique(occumatrix$Aou, header = "Species")
birdsoutputm = c()
for (species in uniq.spp) {
  spec.routes <- occumatrix[(occumatrix$Aou) == species, "stateroute"] #subset routes for each species (i) in tidybirds
  env.sub <- envtable[envtable$stateroute %in% spec.routes, ] #subset routes for each env in tidybirds

  envmeans = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')], 2, mean))
  envsd = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')], 2, sd))

  birdsoutputm = rbind(birdsoutputm, c(species, envmeans, envsd))
  
}
birdsoutput = data.frame(birdsoutputm)

names(birdsoutput) = c("Species", "Mean.Temp", "Mean.Precip", "Mean.Elev", "Mean.NDVI", "SD.Temp", "SD.Precip", "SD.Elev", "SD.NDVI")
### Combine relevant information from each of your two or more datasets using merge()
#(species/occupancy/expected env variables/observed env variables)
occubirds <- merge(birdsoutput, occumatrix, by.x = "Species", by.y="Aou", na.rm = T)
occuenv <- merge(envtable, occubirds, by = "stateroute", all = F, na.rm = T)
occuenv <- na.omit(occuenv)

### Conduct analyses
#Calculating z scores for each environmnetal variable (observed mean - predicted mean/predicted SD)
occuenv$zTemp = (occuenv$mat.mean - occuenv$Mean.Temp) / occuenv$SD.Temp
occuenv$zPrecip = (occuenv$map.mean - occuenv$Mean.Precip) / occuenv$SD.Precip
occuenv$zElev = (occuenv$elev.mean - occuenv$Mean.Elev) / occuenv$SD.Elev
occuenv$zNDVI = (occuenv$ndvi.mean - occuenv$Mean.NDVI) / occuenv$SD.NDVI

write.csv(occuenv, "occuenv.csv", row.names= FALSE)
