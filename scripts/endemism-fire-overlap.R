##****************************************************************************##
# fire severity and weighted endemism overlap for invertebrates
# - point data
# Aaron Greenville
# June 2020
##****************************************************************************##

## >> load libraries ####
#library(phyloregion)
library(rgdal)
library(raster)
library(tidyverse)
library(tmap)
#library(ggplot2)

## >> load endemism results (polygons) ####
# load("data/m1_poly_0-5.rds")
m1 <- m1.poly

## >> load endemism results (points) ####
#load("data/m1_point_1.rds")
load("data/m1_point_0-5.rds")


aust <- readOGR("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\GIS_files\\World boundaries\\Aust\\AustraliaAdmin.shp")
plot(aust)



## >> Load in fire severity raster (re-classed) and get unique classes (raster) ####
# From Payal: Note this is the GEEBAM severity layer that has been reprojected to 250m. sq 
# at Albers equal area and extent to cover islands and offland territories,
# then clipped to NVIS native vegetation and finally clipped to the
# Preliminary Analysis Area

fire_severity <- raster("data/severity5_eqar250_native_paa.tif",
                        RAT = TRUE) #
aust.utm <- spTransform(aust,
                               crs(fire_severity))
analysis.area <- readOGR("data/Preliminary_Analysis_Areas/prelim_analysis_areas_dissolve.shp")
plot(analysis.area)

analysis.area.utm <- spTransform(analysis.area,
                               crs(fire_severity))

plot(aust.utm)
plot(analysis.area.utm, add=TRUE, col = "lightgrey")
plot(aust.utm, add = TRUE)


## >> break endemism value into quantiles and assign to new field ####

# m1$corrected_endemism_quantile <- cut(m1$corrected_endemism , breaks=quantile(m1$corrected_endemism),
#                                          labels=1:4, include.lowest=TRUE)

# break into percentiles
m1$corrected_endemism_quantile <- cut(m1$corrected_endemism , breaks=quantile(m1$corrected_endemism,
                                                                              probs = seq(0, 1, by = 0.20)),
                                      labels=c("0-20%","20-40%", "40-60%", "60-80%", "80-100%"), 
                                      include.lowest=TRUE)


# check
Corrected.weighted.endemism.q <- tm_shape(m1) +
  tm_polygons("corrected_endemism_quantile", 
              #style="pretty", 
              title="Points: Corrected weighted endemism",
              lwd = 0.5,
              palette="OrRd")+ #"-RdBu"
  tm_shape(aust)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)



## >> Extract fire severity values in each endemism polygon ####

plot(fire_severity)  

# convert to same coordinate system  
m1.utm <- spTransform(m1,crs(fire_severity))

# crop to mainland extent
fire.sev.crop <- crop(fire_severity, (analysis.area.utm))  
m1.utm.crop <- crop(m1.utm, (analysis.area.utm))

# replace fire severity class 5 (Very high) with 4 (high) for mapping below
fire.sev.crop.recl <- reclassify(fire.sev.crop, c(4,5,4))

# extract fire severity for each endemism polygon. Note still has all classes for flexibility  
extract.fire <- exactextractr::exact_extract(fire.sev.crop, m1.utm.crop, include_area = TRUE,
                                             include_cols  = "corrected_endemism_quantile") 

# collapse list in df  
extract.fire.all.df<-do.call("rbind",extract.fire )
head(extract.fire.all.df)

## >> calc area per endemism and fire severity ####

endemism.burnt <- extract.fire.all.df %>% 
  mutate(area.adj = area*coverage_fraction) %>%
  group_by(corrected_endemism_quantile, value) %>%
  replace_na(list(value = 0)) %>%
  summarise(burnt.area.km = sum(area.adj)/1000000) %>%
  # commented out for change in fire severity classes
  # mutate(fire.severity = factor(recode(value, 
  #                               '1' = "No data",
  #                               '2' = "Unburnt",
  #                               '3' = "Low and moderate",
  #                               '4' = "High",
  #                               '5' = "Very high",
  #                               '0' = "Unburnt"),
  #                               levels = c("No data",
  #                                          "Unburnt",
  #                                          "Low and moderate",
  #                                          "High",
  #                                          "Very high"))) %>%
  mutate(fire.severity = factor(recode(value, 
                                     '1' = "No data",
                                     '2' = "Unburnt",
                                     '3' = "Low and moderate",
                                     '4' = "High and very high",
                                     '5' = "High and very high",
                                     '0' = "Unburnt"),
                              levels = c("No data",
                                         "Unburnt",
                                         "Low and moderate",
                                         "High and very high"))) %>%
  group_by(corrected_endemism_quantile, fire.severity) %>%
  summarise(burnt.area.km = sum(burnt.area.km)) %>%
  group_by(corrected_endemism_quantile) %>% mutate(percent = burnt.area.km/sum(burnt.area.km)*100)
                        
# plotting overall results  
ggplot(endemism.burnt, aes(factor(corrected_endemism_quantile), burnt.area.km, fill = fire.severity )) +
  geom_col(position = "dodge" ) 

# Plotting burnt areas only

burnt.endemism.area.plot <- endemism.burnt %>%
  filter(!fire.severity %in% c("Unburnt","No data")) %>%
  ggplot(aes(factor(corrected_endemism_quantile), burnt.area.km, fill = fire.severity )) +
  geom_col(position = "dodge" , color = "black") +
  scale_fill_brewer(palette = "OrRd", "Fire severity")+
  scale_x_discrete(name = "Endemism percentile") +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression("Area burnt "~(km^2)))+
  theme_classic()


# cowplot::save_plot("output/endemism_area_recl-0-5-plot.png", burnt.endemism.area.plot, base_height = 6, base_width = 8)

## Map of inputs ####

input.map <- tm_shape(aust.utm) +
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)+
 tm_shape(m1.utm.crop) +
   tm_polygons("corrected_endemism_quantile", 
               #style="quantile", 
               title="Percentile: Corrected weighted endemism",
               lwd = 0.5,
               palette="YlGnBu",
               alpha = 0.5,
               style = "cat")+
 tm_shape(aust.utm) +
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)+
  tm_shape(analysis.area.utm)+
  tm_polygons( "Id",
              alpha = 0,
              legend.show = T,
              lwd = 1, border.col = 4,
              title = "Preliminary analysis area",
              label = "")+
  tm_shape(fire.sev.crop.recl)+
    tm_raster("severity5_eqar250_native_paa",
              title = "Fire severity",
              style = "cat",
               labels = c("No data",
                          "Unburnt",
                          "Low and moderate",
                          "High and very high"))
   
  
#tmap_arrange(Corrected.weighted.endemism, input.map, ncol=2 )

# tmap_save(input.map, filename = "output/tmap_inputMap_recl-0-5.png")


## << Moran I test for spatial autocorr in endemism

library(spdep)

# extract the center of each polygon
coo <- coordinates(m1)

# Search radius to include all neighboring polygon (0 - 200km)
S.dist  <-  dnearneigh(coo, 0, 200000)  

#identify all neighboring polygons for each polygon in the dataset.
lw <- nb2listw(S.dist, style="W",zero.policy=T) 

# Run the MC simulation
MI  <-  moran.mc(m1$corrected_endemism, lw, nsim=9999,zero.policy=T)

