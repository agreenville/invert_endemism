##****************************************************************************##
# fire severity and weighted endemism overlap for invertebrates
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

load("data/m1_point_1.rds")

# load("data/m1_point_0-5.rds")
# load("data/m1_poly_0-5.rds")

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

m1$corrected_endemism_quantile <- cut(m1$corrected_endemism , breaks=quantile(m1$corrected_endemism),
                                         labels=1:4, include.lowest=TRUE)

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


## >> Dissolve polys based on corrected_endemism_quantile ####
m1.dissolve <- aggregate(m1, by = "corrected_endemism_quantile")

  
# check
  tm_shape(m1.dissolve) +
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
m1.dissolve.utm <- spTransform(m1.dissolve,
                            crs(fire_severity))

m1.utm <- spTransform(m1,crs(fire_severity))

# crop to mainland extent
fire.sev.crop <- crop(fire_severity, (analysis.area.utm))  
m1.utm.crop <- crop(m1.utm, (analysis.area.utm))


# extract fire severity for each endemism polygon  
extract.fire <- exactextractr::exact_extract(fire.sev.crop, m1.utm.crop, include_area = TRUE) 

# add endemism col to list elements
extract.fire.all<- Map(cbind,extract.fire, Endemism = c(1:4) )

# collapse list in df  

extract.fire.all.df<-do.call("rbind",extract.fire.all )
head(extract.fire.all.df)

## >> calc area per endemism and fire severity ####

endemism.burnt <- extract.fire.all.df %>% 
  mutate(area.adj = area*coverage_fraction) %>%
  group_by(Endemism, value) %>%
  replace_na(list(value = 0)) %>%
  summarise(burnt.area.km = sum(area.adj)/1000000) %>%
  mutate(fire.severity = factor(recode(value, 
                                '1' = "No data",
                                '2' = "Unburnt",
                                '3' = "Low and moderate",
                                '4' = "High",
                                '5' = "Very high",
                                '0' = "Unburnt"),
                                levels = c("No data",
                                           "Unburnt",
                                           "Low and moderate",
                                           "High",
                                           "Very high"))) %>%
  group_by(Endemism, fire.severity) %>%
  summarise(burnt.area.km = sum(burnt.area.km)) 
                        
# plotting overall results  
ggplot(endemism.burnt, aes(factor(Endemism), burnt.area.km, fill = fire.severity )) +
  geom_col(position = "dodge" ) +
  scale_x_discrete(labels=c("1" = "Q1", "2" = "Q2",
                            "3" = "Q3", "4"="Q4"))

# Plotting burnt areas only

burnt.endemism.area.plot <- endemism.burnt %>%
  filter(!fire.severity %in% c("Unburnt","No data")) %>%
  ggplot(aes(factor(Endemism), burnt.area.km, fill = fire.severity )) +
  geom_col(position = "dodge" , color = "black") +
  scale_fill_brewer(palette = "OrRd", "Fire severity")+
  scale_x_discrete(name = "Endemism quantile", labels=c("1" = "Q1", "2" = "Q2",
                            "3" = "Q3", "4"="Q4")) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression("Area burnt "~(km^2)))+
  theme_classic()

# cowplot::save_plot("output/endemism_area_plot.png", burnt.endemism.area.plot, base_height = 6, base_width = 8)

## Map of inputs ####

input.map <- tm_shape(aust.utm) +
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)+
tm_shape(m1.utm.crop) +
   tm_polygons("corrected_endemism_quantile", 
               #style="quantile", 
               title="Quantile: Corrected weighted endemism",
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
  tm_shape(fire.sev.crop)+
    tm_raster("severity5_eqar250_native_paa",
              title = "Fire severity",
              style = "cat",
               labels = c("No data",
                          "Unburnt",
                          "Low and moderate",
                          "High",
                          "Very high"))
   
  
#tmap_arrange(Corrected.weighted.endemism, input.map, ncol=2 )

# tmap_save(input.map, filename = "output/tmap_inputMap_1.png")




