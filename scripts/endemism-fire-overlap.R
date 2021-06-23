##****************************************************************************##
# fire severity and weighted endemism overlap for invertebrates
# Aaron Greenville
# June 2020
##****************************************************************************##

## >> load libraries ####
#library(phyloregion)
library(rgdal)
library(raster)
library(tidyr)
library(ggplot2)

## >> load endemism results (polygons) ####

load("data/m1_point_1.rds")

load("data/m1_point_0-5.rds")
load("data/m1_poly_0-5.rds")

aust <- readOGR("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\GIS_files\\World boundaries\\Aust\\AustraliaAdmin.shp")
plot(aust)

## >> Load in fire severity raster (re-classed) and get unique classes (raster) ####
# fire_severity <- raster("data/AUS_GEEBAM_Fire_Severity_NIAFED20200224/AUS_GEEBAM_Fire_Severity_NIAFED20200224.tif",
#                         RAT = TRUE) #severity5_eqar250_native_paa.tif

# From Payal: Note this is the GEEBAM severity layer that has been reprojected to 250m. sq 
# at Albers equal area and extent to cover islands and offland territories,
# then clipped to NVIS native vegetation and finally clipped to the
# Preliminary Analysis Area

fire_severity <- raster("data/severity5_eqar250_native_paa.tif",
                        RAT = TRUE) #

analysis.area <- readOGR("data/Preliminary_Analysis_Areas/prelim_analysis_areas_dissolve.shp")
plot(analysis.area)

analysis.area.utm <- spTransform(analysis.area,
                               crs(fire_severity))

# fire.data_dbf <- foreign::read.dbf("data/AUS_GEEBAM_Fire_Severity_NIAFED20200224/AUS_GEEBAM_Fire_Severity_NIAFED20200224.tif.vat.dbf")
# 
# fire_severity.gbam <- reclassify(fire_severity, data.frame(fire.data_dbf$Value, fire.data_dbf$GEEBAMcls))


#aggregate from 40x40 resolution to 252x252 (factor = 6.3)

# current res
# proj = "+proj=lcc +lat_1=-28 +lat_2=-36 +lat_0=-32 +lon_0=135 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
# fire_severity.gbam.lam <- projectRaster(fire_severity.gbam,
#             crs = proj)
# 
# save(fire_severity.gbam.lam, "data/GEEBAM_UTS.rds")

# fire_severity.gbam.aggregate <- aggregate(fire_severity.gbam, fact=6.3)
# res(fire_severity.gbam.aggregate)

## >> break endemism value into quantiles and assign to new field ####

m1$corrected_endemism_quantile <- cut(m1$corrected_endemism , breaks=quantile(m1$corrected_endemism),
                                         labels=1:4, include.lowest=TRUE)

# check
Corrected.weighted.endemism.q <- tm_shape(m1) +
  tm_polygons("corrected_endemism_quantile", 
              #style="pretty", 
              title="Points: Corrected weighted endemism",
              lwd = 0.5,
              palette="OrRd")+ #"-RdBu"
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)


## >> Dissolve polys based on corrected_endemism_quantile ####
m1.dissolve <- aggregate(m1, by = "corrected_endemism_quantile")

  
# check
  tm_shape(m1.dissolve) +
    tm_polygons("corrected_endemism_quantile", 
                #style="pretty", 
                title="Points: Corrected weighted endemism",
                lwd = 0.5,
                palette="OrRd")+ #"-RdBu"
    tm_shape(aust_WGS84)+
    tm_polygons("PLACENAME", 
                alpha = 0,
                legend.show = F,
                lwd = 2, border.col = 1)
  

## >> Extract fire severity values in each endemism polygon ####

plot(fire_severity)  

# convert to same coordinate system  
m1.dissolve.utm <- spTransform(m1.dissolve,
                            crs(fire_severity))

# aust.utm <- spTransform(aust,
#                                crs(fire_severity))

# crop to mainland extent
fire.sev.crop <- crop(fire_severity, extent(analysis.area.utm))  
m1.dissolve.utm.crop <- crop(m1.dissolve.utm, extent(analysis.area.utm))

# extract fire severity for each endemism polygon  
extract.fire <- exactextractr::exact_extract(fire.sev.crop, m1.dissolve.utm.crop, include_area = TRUE) 

# # remove na's
# extract.fire.na <- lapply(extract.fire, na.omit)

# add endemism col to list elements
#extract.fire.na<- Map(cbind,extract.fire.na, Endemism = c(1:4) )

extract.fire.all<- Map(cbind,extract.fire, Endemism = c(1:4) )

# collapse list in df  
#extract.fire.df<-do.call("rbind",extract.fire.na )

extract.fire.all.df<-do.call("rbind",extract.fire.all )
head(extract.fire.all.df)

# calc area per cell

endemism.burnt <- extract.fire.all.df %>% 
  mutate(area.adj = area*coverage_fraction) %>%
  group_by(Endemism, value) %>%
  replace_na(list(value = 0)) %>%
  summarise(burnt.area.km = sum(area.adj)/1000000) %>%
  mutate(fire.severity = recode(value, 
                                '1' = "No data",
                                '2' = "Unburnt",
                                '3' = "Low and moderate",
                                '4' = "High",
                                '5' = "Very high",
                                '0' = "Unburnt")) %>%
  group_by(Endemism, fire.severity) %>%
  summarise(burnt.area.km = sum(burnt.area.km)) 
                        
  
ggplot(endemism.burnt, aes(Endemism, burnt.area.km, fill = fire.severity )) +
  geom_col(position = "dodge" )









