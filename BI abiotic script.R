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
#### ----NDVI ----#####
start_yr <- 2001
end_yr <- 2015
min_num_yrs <- 15
richness_w_env <- get_richness_ts_env_data(start_yr, end_yr, min_num_yrs)

<<<<<<< HEAD
=======
#################################################
#This assumes we want all gimms files that are available. It queries for files
#that are available for download and compares against files already in the gimms
#data directory.
#Returns a list of files to download, which may be length 0.
##################################################
get_gimms_download_list=function(gimms_folder = './data/gimms_ndvi/'){
  available_files_download_path=gimms::updateInventory()
  available_files_name=basename(available_files_download_path)
  
  files_present=list.files(gimms_folder)
  #hdr files are created from some of the gimms processing that we don't want to
  #use here.
  files_present=files_present[!grepl('hdr', files_present)]
  
  to_download=available_files_download_path[! available_files_name %in% files_present]
  
  return(to_download)
}

################################################
#Extract values from a single gimms file given a set of coordinates
################################################
#' @importFrom gimms rasterizeGimms
gimms_file_path = "Z:/GIS/MODIS NDVI/2000_2016/GIMMS"
route_locations = routes.laea
extract_gimms_data=function(gimms_file_path, route_locations){
  #Have to load and extract twice. Once for the actual NDVI, once for the quality flag.
  gimmsRaster=gimms::rasterizeGimms(gimms_file_path, flag=FALSE)
  ndvi=raster::extract(gimmsRaster, route_locations)
  gimmsRaster=gimms::rasterizeGimms(gimms_file_path, flag=TRUE)
  flag=raster::extract(gimmsRaster, route_locations)
  
  year=as.numeric(substr(basename(gimms_file_path), 4,5))
  month=substr(basename(gimms_file_path), 6,8)
  day=substr(basename(gimms_file_path), 11,11)
  
  #Convert the a b to the 1st and 15th
  day=ifelse(day=='a',1,15)
  
  #Convert 2 digit year to 4 digit year
  year=ifelse(year>50, year+1900, year+2000)
  
  return(data.frame(year=year, month=month, day=day, ndvi=ndvi, flag=flag, site_id=route_locations@data$site_id, stringsAsFactors = FALSE))
}

################################################
#Extract the NDVI time series for all bbs routes
#from all years of gimms data
################################################

#' @importFrom dplyr bind_rows
#' @importFrom dplyr %>%
#' @importFrom sp coordinates<-
process_gimms_ndvi_bbs=function(gimms_folder = './data/gimms_ndvi/'){
  
  bbs_data <- get_bbs_data()
  route_locations <- unique(dplyr::select(bbs_data, site_id, long, lat))
  coordinates(route_locations) <- c("long", "lat")
  
  gimms_files=list.files(gimms_folder, full.names = TRUE)
  #hdr files are created from some of the gimms processing that we don't want to
  #use here.
  gimms_files=gimms_files[!grepl('hdr', gimms_files)]
  
  gimms_ndvi_bbs=data.frame()
  for(file_path in gimms_files){
    gimms_ndvi_bbs=extract_gimms_data(file_path, route_locations) %>%
      bind_rows(gimms_ndvi_bbs)
  }
  save_provenance(gimms_ndvi_bbs)
  return(gimms_ndvi_bbs)
}


src_sqlite <- function(path, create = FALSE) {
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("RSQLite package required to connect to sqlite db", call. = FALSE)
  }
  
  if (!create && !file.exists(path)) {
    stop("Path does not exist and create = FALSE", call. = FALSE)
  }
  
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  RSQLite::initExtension(con)
  
  src_sql("sqlite", con, path = path)
}

get_bbs_gimms_ndvi = function(
  sqlite_db_file='./data/bbsforecasting.sqlite',
  database = src_sqlite(sqlite_db_file, create = TRUE),
  gimms_folder = './data/gimms_ndvi/'
){
  dir.create(gimms_folder, showWarnings = FALSE, recursive = TRUE)
  
  if('gimms_ndvi_bbs_data' %in% dplyr::src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from gimms_ndvi_bbs_data')), n = Inf))
  } else {
    print('Gimms NDVI bbs data not found, processing from scratch')
    
    files_to_download=get_gimms_download_list()
    if(length(files_to_download)>0){
      print('Downloading GIMMS data')
      downloadGimms(x=files_to_download, dsn=gimms_folder)
    }
    
    gimms_ndvi_bbs_data=process_gimms_ndvi_bbs()
    
    gimms_ndvi_bbs_data=filter_gimms_data(gimms_ndvi_bbs_data)
    
    copy_to(database, gimms_ndvi_bbs_data, temporary = FALSE,
            indexes = list(c('site_id','year','month')))
    
    return(gimms_ndvi_bbs_data)
  }
}
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
#tiffs  = tiffs[1202:2912]
files<-paste('Z:/GIS/EVI/',tiffs,sep='')

evimeans2<-stack(files)
# create a raster with that extent, and the number of rows and colums to achive a
# similar resolution as you had before, you might have to do some math here....
# as crs, use the same crs as in your rasters before, from the crs slot
s<-raster(e, nrows=2906, ncols=2906, crs=files@crs)

# use this raster to reproject your original raster (since your using the same crs,
# resample should work fine
r <- raster(nrow=2906, ncol=2)
r[] <- 1:ncell(r)
s <- raster(nrow=10, ncol=10)
r1<-resample(r1, s, method="bilinear")


files<-paste('Z:/GIS/EVI/',tiffs,sep='')
is.error <- function(x) inherits(x, "try-error")

test = c()
for(file in files){
  w <- raster(file, crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
# evimeans<-stack(files)
  #e = extent(w)
# OutExtent = is.error(try(clip1 <- crop(imported_raster, extent(single))))
  #rbind(test, e)
  plot(w)
}


result = tryCatch({
  evimeans<-stack(files)
}, warning = function(x) {
  suppressWarnings(na.pass(files))
}, finally = {
  cleanup-code
})
  
raster_data <- list.files(pattern='\\.tif$', full.names=TRUE)
w <- raster("MOD13A3.A2000092.h01v07.006.2015137040355.tif", crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
  
for (i in 1:length(raster_data)){
    r <- raster(raster_data[i], crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
    rc = crop(r,e)
    rp <- projectRaster(from = rc, to = w,
                        filename = file.path ("./crop", raster_data[i]),
                        method = "bilinear",
                        format = "raster",
                        overwrite = TRUE)
  }

 
# error @ 844, 1006, 1007, 2065-2067, 2113 "MOD13A3.A2012122.h12v05.006.2015246095954.tif"



>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
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

<<<<<<< HEAD

# NDVI 
gimms_ndvi = read.csv("ENV DATA/gimms_ndvi_bbs_data.csv", header = TRUE)
gimms_agg = gimms_ndvi %>% filter(month == c("may", "jun", "jul")) %>% 
  group_by(site_id, year, month)  %>%  summarise(ndvi.mean=mean(ndvi))
gimms_agg$stateroute = gimms_agg$site_id
ndvi = gimms_agg[,c("stateroute", "ndvi.mean")]
# merge together
all_env = Reduce(function(x, y) merge(x, y, by = "stateroute"), list(env_mat, env_elev, env_map, ndvi))

# write.csv(all_env,'C:/git/Biotic-Interactions/all_env.csv',row.names=F)
=======
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
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
####----Creating an environmental matrix ----####
occumatrix <- read.csv("2001_2015_bbs_occupancy.csv", header = T) # read in updated bbs data
route.locs = read.csv('latlong_rtes.csv')

latlongs = subset(route.locs, select = c('stateroute', 'latitude', 'longitude'))
latlongs$latitude = abs(latlongs$latitude)
route.sp = coordinates(latlongs[,3:2])

#mean data for all variables
<<<<<<< HEAD
envtable <- subset(all_env, select = c('stateroute', 'mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')) 
=======
envtable <- subset(all_env, select = c('stateroute', 'mat.mean', 'map.mean', 'elev.mean', 'evi.mean')) 
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2

### Calculate metrics of interest
####---- Creating final data frame----####
#For loop to calculate mean & standard dev environmental variables for each unique species
uniq.spp = unique(occumatrix$Aou, header = "Species")
birdsoutputm = c()
for (species in uniq.spp) {
  spec.routes <- occumatrix[(occumatrix$Aou) == species, "stateroute"] #subset routes for each species (i) in tidybirds
  env.sub <- envtable[envtable$stateroute %in% spec.routes, ] #subset routes for each env in tidybirds
<<<<<<< HEAD
  envmeans = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')], 2, mean))
  envsd = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'ndvi.mean')], 2, sd))
=======
  envmeans = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'evi.mean')], 2, mean))
  envsd = as.vector(apply(env.sub[, c('mat.mean', 'map.mean', 'elev.mean', 'evi.mean')], 2, sd))
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
  
  birdsoutputm = rbind(birdsoutputm, c(species, envmeans, envsd))
  
}
birdsoutput = data.frame(birdsoutputm)
<<<<<<< HEAD
names(birdsoutput) = c("Species", "Mean.Temp", "Mean.Precip", "Mean.Elev", "Mean.NDVI", "SD.Temp", "SD.Precip", "SD.Elev", "SD.NDVI")
=======
names(birdsoutput) = c("Species", "Mean.Temp", "Mean.Precip", "Mean.Elev", "Mean.EVI", "SD.Temp", "SD.Precip", "SD.Elev", "SD.EVI")
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2

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
<<<<<<< HEAD
occuenv$zNDVI = (occuenv$ndvi.mean - occuenv$Mean.NDVI) / occuenv$SD.NDVI

write.csv(occuenv, "Z:/Snell/occuenv.csv", row.names= FALSE)
=======
occuenv$zEVI = (occuenv$evi.mean - occuenv$Mean.EVI) / occuenv$SD.EVI

write.csv(occuenv, "occuenv_new.csv", row.names= FALSE)
>>>>>>> 8e9689638283c792a876fc521cf5e1ee1ac6a7f2
