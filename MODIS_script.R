###This script imports environmental data for BBS routes

library(raster)

# Define lat/long window for raster data
box = c(-170,-54,24,74) # North America

# Fill in which box you want to use
my.extent = extent(box) 

# Define projection to be used throughout analysis
prj.string <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"

####################################################################################
#### Import route locations and draw sample circles around them
library('sp')
library('rgdal')

# derived from BBS_occ script
routes = read.csv("latlong_rtes.csv",header =TRUE)
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

# Check that circle loactions look right
plot(circs.sp)

#####################################################################################
### Environmental Data Layers

setwd('Z:/GIS/MODIS NDVI')

env.data = raster('Vegetation_Indices_may-aug_2000-2010')


# Define the projection of the raster layer (this may be different for different data)
#  See documentation in PROJ4
projection(env.data) = CRS("+proj=longlat +ellps=WGS84") 


# Make a dataframe to store environmental data
#   Adds the first column atPoint by extracting values found exactly at the route coordinates

# Extract Data
env.point = raster::extract(env.data, routes)
env.mean = raster::extract(env.data, circs.sp, fun = mean, na.rm=T)
env.var = raster::extract(env.data, circs.sp, fun = mean, na.rm=T)

# Put into dataframe
env = data.frame(stateroute = names(circs.sp), env.point = env.point, env.mean = env.mean, env.var = env.var)



write.csv(env,'your_new_data_file.csv',row.names=F)



