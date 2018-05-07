
setwd("C:/git/core-transient")

library(lme4)
library(plyr) # for core-transient functions
library(ggplot2)
library(tidyr)
library(maps)
library(gridExtra)
library(RColorBrewer)
library(sp)
library(rgdal)
library(raster)
library(dplyr)
library(digest)


source('scripts/R-scripts/core-transient_functions.R')

#############################################

# Specify here the datasetIDs and then run the code below.
dataformattingtable = read.csv('data_formatting_table.csv', header = T) 

datasetIDs = dataformattingtable$dataset_ID[dataformattingtable$format_flag == 1]

# BBS (dataset 1) will be analyzed separately for now.
datasetIDs = datasetIDs[!datasetIDs %in% c(1)]

# rbind in new BBS data ## bbs_be_lat has lots of rows we don't use
lat_scale_rich = read.csv("output/tabular_data/lat_scale_rich.csv", header = TRUE)
lat_scale_bbs = subset(lat_scale_rich, datasetID == 1)
lat_scale_bbs$spRich = NA
lat_scale_bbs$nTime = NA
lat_scale_bbs$scale = 1
lat_scale_bbs = lat_scale_bbs[,c("datasetID","site","spRich","nTime","meanAbundance", "scale")]

# Plotting summary results across datasets for Core-Transient analysis
summ = read.csv('output/tabular_data/core-transient_summary.csv', header=T) #87 unique ids

allsummaries = read.csv("output/tabular_data/allsummaries.csv", header = TRUE) #39 unique ids


# BBS (dataset 1) will be analyzed separately for now.
summ2 = summ[!summ$datasetID %in% allsummaries$datasetID,] #now 48 unique ids
summ2$spRich = summ2$spRichTotal
summ2$scale = 1
summ2 = summ2[,c("datasetID","site","spRich","nTime","meanAbundance", "scale")]

scalesum = rbind(allsummaries, summ2) #87 datasets


all_assemblages = rbind(scalesum, lat_scale_bbs)
median(table(all_assemblages$datasetID))




