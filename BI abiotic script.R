###This script imports environmental data for BBS routes

library(raster)
library(sp)
library(rgdal)
library(RCurl)

# Define lat/long window for raster data
box = c(-170,-54,24,74) # North America

# Fill in which box you want to use
my.extent = extent(box) 

# Define projection to be used throughout analysis
prj.string <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"

####################################################################################
#### Import route locations and draw sample circles around them

# derived from BBS_occ script
routes = read.csv("latlong_rtes.csv",header =TRUE)
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
#### ----EVI ----#####
setwd('Z:/GIS/EVI')

#url <- "https://search.earthdata.nasa.gov/granules/download.html?project=5135767041&collection=C194001239-LPDAAC_ECS"
#filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE)

#filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
#filenames <- strsplit(filenames, "\r\n")
#filenames = unlist(filenames)

# downloads all NDVI files afrom FTP site and saves to MODIS NDVI folder
if(TRUE){for (filename in filenames) {
  download.file(paste(url, filename, sep = ""), paste(getwd(), "/", filename, sep = ""))
}
}


r.st <- stack(lapply(fnms, function(x) raster(x)))  
img <- raster('Z:/GIS/MODIS NDVI/2000_2016/2000/MOD13A3.A2000122.h01v07.005.2007111071631.hdf')

bn1 <- "MOD13A3.A2000122.h01v07.005.2007111071631.hdf"
b01 <- raster(readGDAL(paste("Z:/GIS/MODIS NDVI/2000_2016/2000/",bn1, sep = ""), silent = TRUE))


#### skip above ####

# need to convert HDFs into TIFFs, EVI is band 2
library(gdalUtils)
fnms <- list.files(path= "Z:/GIS/EVI/NASA_NEW", pattern="*.hdf")
fnms = fnms[2117:2912]
for(i in fnms){
  hdf4_dataset <- paste('Z:/GIS/EVI/NASA_NEW/', i, sep = "")
  i = gsub(".hdf", ".tif", i)
  now = gdal_translate(hdf4_dataset, paste('Z:/GIS/EVI/', i, sep = ""),sd_index=2)
} 
  



#### VRT business
gdal_setInstallation()
valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
if(valid_install)
{
  tiffs <- list.files(path= "Z:/GIS/EVI", pattern="*.tif")
  output.vrt <- paste(tempfile(),".vrt",sep="")
  gdalbuildvrt(gdalfile=tiffs,output.vrt=output.vrt,separate=TRUE,verbose=TRUE)
  gdalinfo(output.vrt)
}



#### from Tracie CC script
e = raster(extent(15.961,69.9,-164.883,-56.953))
tiffs <- list.files(path= "Z:/GIS/EVI", pattern="*.tif")

# create a raster with that extent, and the number of rows and colums to achive a
# similar resolution as you had before, you might have to do some math here....
# as crs, use the same crs as in your rasters before, from the crs slot
s<-raster(e, nrows=2906, ncols=2906, crs=files@crs)

# use this raster to reproject your original raster (since your using the same crs,
# resample should work fine
r1<-resample(r1, s, method="ngb")


files<-paste('Z:/GIS/EVI/',tiffs,sep='')
evimeans<-stack(files)

writeRaster(template, file="MyBigNastyRasty.tif", format="GTiff")
mosaic_rasters(gdalfile=tiffs,dst_dataset="MyBigNastyRasty.tif",of="GTiff")
gdalinfo("MyBigNastyRasty.tif")
for(j in tiffs){
  tif<-readTIFF(j)
}
 
# error @ 844, 1006, 1007, 2065-2067, 2113 "MOD13A3.A2012122.h12v05.006.2015246095954.tif"



setwd('C:/Git/Biotic-Interactions/ENV DATA')
# Read in stack of layers from all 12 months 
tmean <- getData("worldclim", var = "tmean", res = 10)
files<-paste('tmean',1:12,'.bil',sep='')
tmeans<-stack(files)
mat = calc(tmeans, mean)

#### ----Precip ----#####
#read in precip data from world clim, stack data to get 1 MAP value
# Read in stack of layers from all 12 months
prec <- getData("worldclim", var = "prec", res = 10)
pfiles<-paste('prec',1:12,'.bil',sep='')
pmeans<-stack(pfiles)
map = calc(pmeans, mean)

#### ----Elev ----#####
#read in elevation data from world clim
elev <- getData("worldclim", var = "alt", res = 10)
alt_files<-paste('alt_10m_bil', sep='')


setwd('C:/Git/Biotic-Interactions')

# Define the projection of the raster layer (this may be different for different data)
#  See documentation in PROJ4
projection(evi.proj) 

  
#for(env.i in env){
# Extract Data
mat.point = raster::extract(mat, routes)
mat.mean = raster::extract(mat, circs.sp, fun = mean, na.rm=T)
mat.var = raster::extract(mat, circs.sp, fun = var, na.rm=T)
#env.data = rbind(env.data, c(env.i,env.point, env.mean,env.var))
#}
env_mat = data.frame(stateroute = names(circs.sp), mat.point = mat.point, mat.mean = mat.mean, mat.var = mat.var)

env_mat = read.csv("env_mat.csv", header = TRUE)
# Extract Data
elev.point = raster::extract(elev, routes)
elev.mean = raster::extract(elev, circs.sp, fun = mean, na.rm=T)
elev.var = raster::extract(elev, circs.sp, fun = var, na.rm=T)

env_elev = data.frame(stateroute = names(circs.sp), elev.point = elev.point, elev.mean = elev.mean, elev.var = elev.var)
env_elev=read.csv("env_elev.csv", header = TRUE)

# Extract Data
map.point = raster::extract(map, routes)
map.mean = raster::extract(map, circs.sp, fun = mean, na.rm=T)
map.var = raster::extract(map, circs.sp, fun = var, na.rm=T)

env_map = data.frame(stateroute = names(circs.sp), map.point = map.point, map.mean = map.mean, map.var = map.var)
env_map = read.csv("env_map.csv", header = TRUE)

# Extract EVI Data
evi.point = raster::extract(evi.proj, routes.laea)
evi.mean = raster::extract(evi.proj, circs.sp, fun = mean, na.rm=T)
evi.var = raster::extract(evi.proj, circs.sp, fun = var, na.rm=T)

# Put into dataframe
env_evi = data.frame(stateroute = names(circs.sp), evi.point = evi.point, evi.mean = evi.mean, evi.var = evi.var)
env_evi = read.csv("env_evi.csv", header = TRUE)
# toplot = filter(env_evi, evi.mean > 0) 
#env_evi[is.na(env_evi)] <- 0
# merge together
all_env = Reduce(function(x, y) merge(x, y, by = "stateroute"), list(env_mat, env_elev, env_map, env_evi))

# write.csv(all_env,'C:/git/Biotic-Interactions/all_env.csv',row.names=F)
plot(evi.proj)
points(toplot$longitude, toplot$latitude)
####----Creating an environmental matrix ----####
occumatrix <- read.csv("2001_2015_bbs_occupancy.csv", header = T) # read in updated bbs data
route.locs = read.csv('latlong_rtes.csv')

latlongs = subset(route.locs, select = c('stateroute', 'latitude', 'longitude'))
latlongs$latitude = abs(latlongs$latitude)
route.sp = coordinates(latlongs[,3:2])

#mean data for all variables
envtable <- subset(all_env, select = c('stateroute', 'mat.mean', 'map.mean', 'elev.mean', 'evi.mean')) 

### Calculate metrics of interest
####---- Creating final data frame----####
#For loop to calculate mean & standard dev environmental variables for each unique species
uniq.spp = unique(occumatrix$Aou, header = "Species")
birdsoutputm = c()
for (species in uniq.spp) {
  spec.routes <- occumatrix[(occumatrix$Aou) == species, "stateroute"] #subset routes for each species (i) in tidybirds
  env.sub <- envtable[envtable$stateroute %in% spec.routes, ] #subset routes for each env in tidybirds
  envmeans = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'evi.mean')], 2, mean))
  envsd = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'evi.mean')], 2, sd))
  
  birdsoutputm = rbind(birdsoutputm, c(species, envmeans, envsd))
  
}
birdsoutput = data.frame(birdsoutputm)
names(birdsoutput) = c("Species", "Mean.Temp", "Mean.Precip", "Mean.Elev", "Mean.EVI", "SD.Temp", "SD.Precip", "SD.Elev", "SD.EVI")

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
occuenv$zEVI = (occuenv$evi.mean - occuenv$Mean.EVI) / occuenv$SD.EVI

write.csv(occuenv, "occuenv_new.csv", row.names= FALSE)
