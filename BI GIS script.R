#### BI analysis - Use GIS to get species with overlapping ranges, calculate centroids of range, 
#    and generate expected presence of landbird species for pairwise occupancy comparison

#### Subset of BI Data Cleaning Code
library(maps)
library(rgdal)
library(shapefiles)
library(maptools)
library(raster)
library(rgeos)
library(gtools)
library(sp)
library(tidyr)
library(dplyr)

# setwd("C:/git/Biotic-Interactions")
pi_catch = function(pi) {
  if (!is.null(pi)){
    pi_out = tryCatch( #### TRYCATCH
      
      {
        suppressWarnings(pi) 
      },
      warning = function(cond) {
        test = gArea(pi)
      },
      error = function(cond) {
        return(NA)
      })
    return(gArea(pi))
  }
  else area_overlap = NA
}

# read in temporal occupancy data from BI occ script
temp_occ = read.csv("data/bbs_sub1.csv", header=TRUE)
temp_occ$Aou[temp_occ$Aou == 4810] = 4812
# read in lat long data
bbs_routes = read.csv("data/latlong_rtes.csv",header =TRUE)
# read in bird range shps
shapefile_path = 'Z:/GIS/birds/All/All/'
# on mac shapefile_path = '/Volumes/hurlbertlab/GIS/birds/All/All'
all_spp_list = list.files(shapefile_path)



# read in new_spec_weights file created in data cleaning code
new_spec_weights=read.csv("data/new_spec_weights.csv", header=TRUE)
new_spec_weights$focalAOU = as.numeric(new_spec_weights$focalAOU)
new_spec_weights$compAOU = as.numeric(new_spec_weights$CompAOU)
new_spec_weights = unique(new_spec_weights)
# for loop to select a genus_spp from pairwise table, read in shp, subset to permanent habitat, plot focal distribution
filesoutput = c()
# dropping non-intersecting polygons
focal_spp = unique(new_spec_weights$focalcat)

intl_proj = CRS("+proj=longlat +datum=WGS84")
sp_proj = CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km")

####### for loop generating shapefiles and area table for all spp - DO NOT RUN! ######
if(TRUE) {  #Blocking out the for loop below. Need to change to TRUE if you want the loop to run.
  
  for (sp in focal_spp) {
    print(sp)
    t1 = all_spp_list[grep(sp, all_spp_list)]
    t2 = t1[grep('.shp', t1)]
    t3 = strsplit(t2, ".shp")
    test.poly <- readShapePoly(paste(shapefile_path, t3, sep = "")) # reads in species-specific shapefile
    # test.poly <- rgdal::readOGR(dsn = shapefile_path, layer = t3) # reads in species-specific shapefile
    proj4string(test.poly) <- intl_proj
    colors = c("blue", "yellow", "green", "red", "purple")
    # subset to just permanent or breeding residents
    sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
    sporigin = spTransform(sporigin, CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km"))
    plot(sporigin, col = colors, border = NA) 
    gArea(spTransform(sporigin, CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km")))
    
    # list this focal spp competitor
    tmp = filter(new_spec_weights, sp == new_spec_weights$focalcat)
    comp_spp = tmp$compcat
    
    # match competitor sp to focal spp, intersect its range with the focal range,
    # and calcualte the area of overlap between the two species.
    for(co in comp_spp) {          
      c1 = all_spp_list[grep(co, all_spp_list)]
      c2 = c1[grep('.shp', c1)]
      c3 = strsplit(c2, ".shp")
      comp.poly <- readShapePoly(paste(shapefile_path, c3, sep = "")) # reads in species-specific shapefile
      proj4string(comp.poly) <- intl_proj
      corigin = comp.poly[comp.poly@data$SEASONAL == 1|comp.poly@data$SEASONAL == 2|comp.poly@data$SEASONAL ==5,]
      corigin = spTransform(corigin,  CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km"))
      plot(corigin, add = TRUE ,col = colors, border = NA) 
      # intersect from raster package
      sporigin = gBuffer(sporigin, byid=TRUE, width=0)
      corigin = gBuffer(corigin, byid=TRUE, width=0)
      pi = raster::intersect(sporigin, corigin)
      spArea = gArea(sporigin) # in m
      coArea = gArea(corigin)
      area_overlap = pi_catch(pi)
      focalAOU = unique(new_spec_weights[new_spec_weights$focalcat == sp, c('focalAOU')])
      compAOU = unique(new_spec_weights[new_spec_weights$compcat == co, c('CompAOU')])
      filesoutput = rbind(filesoutput, c(sp, focalAOU, co, compAOU, spArea, coArea, area_overlap))
      
    }
  } 
  

  
  filesoutput = data.frame(filesoutput)
  colnames(filesoutput) = c("Focal", "focalAOU","Competitor", "compAOU","FocalArea", "CompArea", "area_overlap")
  # string split to get sci name with spaces
  filesoutput$Focal = gsub('_',' ',filesoutput$Focal)
  filesoutput$Competitor = gsub('_',' ',filesoutput$Competitor)
  write.csv(filesoutput, file = "data/shapefile_areas.csv", row.names = FALSE)
}

######## Calculating centroids for each species - using whole range #####
centroid = c()
new_spec_weights$focalcat = as.character(new_spec_weights$focalcat)
if(FALSE) {for (sp in focal_spp){
  print(sp)
  t1 = all_spp_list[grep(sp, all_spp_list)]
  t2 = t1[grep('.shp', t1)]
  t3 = strsplit(t2, ".shp")
  
  test.poly <- readShapePoly(paste(shapefile_path, t3, sep = "")) # reads in species-specific shapefile
  proj4string(test.poly) <- intl_proj
  colors = c("blue", "yellow", "green", "red", "purple")
  # subset to just permanent or breeding residents
  sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
  sporigin = spTransform(sporigin, CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km"))
  #plot(sporigin, col = colors, border = NA)
  trueCentroid = gCentroid(sporigin)
  coord = coordinates(spTransform(trueCentroid, CRS("+proj=longlat +datum=WGS84")))
  focalAOU = unique(new_spec_weights[new_spec_weights$focalcat == sp, c('focalAOU')])
  centroid = rbind(centroid, c(sp, focalAOU, coord))
}
}
centroid = data.frame(centroid)
names(centroid) = c("Species", "FocalAOU", "Long", "Lat")
centroid$Lat = as.numeric(paste(centroid$Lat))
centroid$Long = as.numeric(paste(centroid$Long))
write.csv(centroid,"data/centroid.csv",row.names=FALSE)

######## Calculating expected presences for each species using whole range #####
bbs_routes$latitude = abs(bbs_routes$latitude)
bbs_routes = dplyr::select(bbs_routes, -c(statenum, route))
routes = unique(bbs_routes$stateroute)

coordinates(bbs_routes) <- c("longitude", "latitude")
proj4string(bbs_routes) <- intl_proj

expect_pres = c()

file_names = dir('sp_routes/')
setwd("sp_routes/")

if(FALSE) {for (sp in focal_spp){
  print(sp)
  
  focalAOU = subset(new_spec_weights, focalcat == sp)
  spAOU = unique(focalAOU$focalAOU)
  
  t1 = all_spp_list[grep(sp, all_spp_list)]
  t2 = t1[grep('.shp', t1)]
  t3 = strsplit(t2, ".shp")
  
  test.poly <- readShapePoly(paste(shapefile_path, t3, sep = "")) # reads in species-specific shapefile
  proj4string(test.poly) <- intl_proj
  sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
  
  plot(sporigin)
  # source: http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
  routes_inside <- bbs_routes[!is.na(sp::over(bbs_routes,as(sporigin,"SpatialPolygons"))),]
  plot(routes_inside, add = T)
  routes_inside = data.frame(routes_inside)
  routes_inside = cbind(routes_inside, spAOU)
  expect_pres=rbind(expect_pres, routes_inside)
}
}
setwd("C:/Git/Biotic-Interactions")
expect_pres = data.frame(expect_pres)
  
expect_pres = dplyr::select(expect_pres, -optional)

# write.csv(expect_pres,"C:/Git/Biotic-Interactions/data/expect_pres.csv",row.names=FALSE)

######## PDF of each species BBS occurrences ########
focalcompsub = read.csv("data/focalcompsub.csv", header=TRUE)
# merge in lat/long
latlongs = read.csv('data/routes 1996-2010 consecutive.csv', header = T)
plotdata_all = merge(focalcompsub, latlongs, by = "stateroute") 

subfocalspecies = unique(focalcompsub$FocalAOU)
# Making pdf of ranges for each focal spp
pdf('Plots_RangeMaps.pdf', height = 8, width = 10)
par(mfrow = c(3, 4))

for(sp in subfocalspecies){ 
  print(sp)
  plotsub = plotdata_all[plotdata_all$FocalAOU == sp,]
  map("state") 
  points(plotsub$Longi, plotsub$Lati, col = 3,  pch = 20)
  title(main = (unique(plotdata_all$FocalSciName[plotdata_all$FocalAOU == sp])))
}

dev.off() 



###### Talk plot 

sp = "Pipilo_maculatus"

focalAOU = subset(new_spec_weights, focalcat == sp)
spAOU = unique(focalAOU$focalAOU)

#file_names = dir('sp_routes/')


t1 = all_spp_list[grep(sp, all_spp_list)]
t2 = t1[grep('.shp', t1)]
t3 = strsplit(t2, ".shp")

test.poly <- readShapePoly(paste(shapefile_path, t3, sep = "")) # reads in species-specific shapefile
proj4string(test.poly) <- intl_proj
sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]


# source: http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
routes_inside <- bbs_routes[!is.na(sp::over(bbs_routes,as(sporigin,"SpatialPolygons"))),]
# plot(routes_inside, add = T)

setwd("C:/git/Biotic-Interactions")
r = read.csv("C:/git/Biotic-Interactions/all_expected_pres.csv", header =TRUE)
r = subset(r,Focal == "Spotted Towhee")

occ_loc = merge(r[,c("stateroute", "FocalOcc", "FocalAbundance","Focal")],routes_inside, by = "stateroute")
occ_sub = subset(occ_loc, Focal = "Spotted Towhee")

colfunc <- colorRampPalette(c("blue","red"))
colfunc(16)
occ_sub$Col <- colfunc(16)[as.numeric(cut(occ_sub$FocalOcc,16))]
occ_sub  = occ_sub %>%
  arrange(FocalOcc)

occ_sub$lab = signif(occ_sub$FocalOcc, 2)

setwd("Z:/GIS/geography")
usa = readShapePoly("Z:/GIS/geography/continent.shp")
proj4string(usa) <- intl_proj
usa = spTransform(usa, CRS("+proj=longlat +datum=WGS84"))

pdf('Figures/spottedtow_plot.pdf', height = 8, width = 11)
par(mfrow = c(1, 1), mar = c(6, 6, 1, 1), mgp = c(4, 1, 0), 
    cex.axis = 1.5, cex.lab = 2, las = 1)

plot(usa, col=rgb(0, 1, 0, 0.2)) #extent(usa, -173.4401,-57.92498,15.11517,77.79797)

plot(sporigin, add=TRUE, col = "gray", border = "gray")
setwd("C:/git/Biotic-Interactions")
points(occ_sub$longitude, occ_sub$latitude, col = unique(occ_sub$Col), pch=20, cex = 0.3)
dev.off()

colfunc <- colorRampPalette(c("gold","red"))
colfunc(10)
occ_sub$Cola <- colfunc(10)[as.numeric(cut(occ_sub$FocalAbundance,10))]
occ_sub  = occ_sub %>%
  arrange(FocalAbundance)

occ_abun$lab = signif(occ_abun$FocalAbundance, 2)
points(occ_sub$longitude, occ_sub$latitude, col = occ_sub$Cola, pch=20, cex = 0.3)

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.5, y = seq(0,1,l=5), labels = seq(1,0,l=5))
rasterImage(legend_image, 0, 0, 1,1)
