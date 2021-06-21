##****************************************************************************##
# fire severity and weigthed endemism overlap for invertebrates
# Aaron Greenville
# June 2020
##****************************************************************************##

## >> load libraries ####
#library(phyloregion)
library(rgdal)
library(raster)

## >> load endemism results (polygons) ####
load("data/m1_point_0-5.rds")
load("data/m1_poly_0-5.rds")


## >> Load in fire severity raster (re-classed) and get unique classes (raster) ####
fire_severity <- raster("data/AUS_GEEBAM_Fire_Severity_NIAFED20200224/AUS_GEEBAM_Fire_Severity_NIAFED20200224.tif") #severity5_eqar250_native_paa.tif


## >> Extract fire severity values in each endemism polygon ####

point.overlap <- raster::extract(fire_severity, m1)
