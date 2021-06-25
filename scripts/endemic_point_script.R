##****************************************************************************##
# Calc weighed endemism for invertebrates point data
# Aaron Greenville
# March 2020
##****************************************************************************##

## >> Load packages ####
library(phyloregion)
library(rgdal)
library(raster)
library(tmap)



## >> load data ####

aust = readOGR("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\GIS_files\\World boundaries\\Aust\\AustraliaAdmin.shp")
plot(aust)

#change datum for aust to match inverts and add
aust_WGS84 <- spTransform(aust,
                          crs("+proj=longlat +datum=WGS84"))

# Notes from Payal:
# 1.	The data now has 58821 unique species as 
# per field spfile. Please note that this is the
# field to be used for unique species identification. 
# The number of unique species by scientificName will be 
# less because there are duplicates. Spfile instead indicates unique 
# IDs by each scientificName-class-family combination.
# 2.	I now have the data in WGS84 
# data was masked using 1km resolution WGS84 mask to Australia

invert.point <- read.csv("data/data_ALAnonALA_wgs84_corrected.csv",
                          header = TRUE)
  
# old data  read.csv("polygons_Aaron/data_ALAnonALA_wgs84.csv",
# header = TRUE)

head(invert.point)



## calc community matrix ####
comm.point <- points2comm(invert.point, lon = "longitude",
            lat = "latitude",
            species = "spfile",
            res = 0.5, # size of grid in decimal degrees 0.1 = ~11km. 0.2 min I can go
            trace=1,
            shp.grids = NULL)  

head(comm.point)



## calc weighed endemism ####
Endm.invert <- weighted_endemism(comm.point$comm_dat)
head(Endm.invert)

# join results back to spatial community data
m1 <- merge(comm.point$poly_shp, data.frame(grids=names(Endm.invert),
                                            WE=Endm.invert), by="grids")
m1 <- m1[!is.na(m1@data$WE),]



## calc corrected weighted endemism #### 
# (weighted endemism tally per cell divided by the species richness of that cell)

m1$corrected_endemism <- m1$WE/m1$richness

# Set coord system for inverts
crs(comm.point$poly_shp) <- "+proj=longlat +datum=WGS84"
crs(m1) <- "+proj=longlat +datum=WGS84"



## Plotting ####

#save plot
# png(file="output/endemism_inverts_point_updated-0-5.png",
#    width=1000, height=600)

# plot  c(bottom, left, top, right)
par(mfrow = c(2, 2), mar=c(0, 5.5, 0, 0), oma=c(0, 0, 0, 0), xpd=TRUE)
  plot_swatch(m1, values = m1$WE, k=10, leg = 5 ,border = NA,
            key_label = "Weighted endemism",
            pos = "bottomleft")#

  plot(aust_WGS84, add=TRUE)

  #species richness
  plot_swatch(m1, values = m1$richness, k=10, leg = 5 ,
            border = NA, key_label = "Species richness")
  plot(aust_WGS84, add=TRUE)

  #corrected weighted endemism
  plot_swatch(m1, values = m1$corrected_endemism, k=10, 
            leg = 5 ,border = NA, key_label = "Corrected weighted endemism",
            breaks = "quantile")

  plot(aust_WGS84, add=TRUE)
par(mfrow = c(1, 1), mar=c(5, 4, 4, 2) + 0.1, xpd=FALSE)

#dev.off()

# library(ggplot2)
# 
# m1.2 <- m1
# m1.2@data$id <- rownames(m1.2@data)
# gg.data <- fortify(m1.2, region = "id" )
# 
# gg.data.df <- merge(gg.data, m1.2@data, by = "id")
# 
# breaks <- quantile(gg.data.df$corrected_endemism)
# 
# aust_WGS84.sf <- st_as_sf(aust_WGS84)
# 
# p0 <- ggplot() +
#   geom_polygon(data = gg.data.df, aes(x = long, y = lat, group = group, fill = corrected_endemism)) +
#   #geom_path(color = "white", size = 0.2) +
#   scale_fill_fermenter(breaks = breaks, palette = "BrBG", direction=1) + 
#   geom_sf(aust_WGS84.sf, aes(x = geometry ,group = PLACENAME, fill = NA )) +
#   coord_equal() +
#   theme(panel.background=element_blank())+
#   theme(panel.background= element_rect(color="black")) +
#   theme(axis.title = element_blank(), axis.text = element_blank()) +
#   labs(title = "Corrected weighted endemism")
# p0


#fire_sf <- st_read("data/NIAFED_v20200623")



## Nicer plots ####

library(tmap)

spRh <- tm_shape(m1) +
  tm_polygons("richness", 
              style="quantile", 
              title="Species Richness",
              lwd = 0.5,
              palette="YlGnBu")+
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)

weighted.endemism <- tm_shape(m1) +
  tm_polygons("WE", 
              style="quantile", 
              title="Weighted endemism",
              lwd = 0.5,
              palette="YlGnBu")+
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)

Corrected.weighted.endemism <- tm_shape(m1) +
  tm_polygons("corrected_endemism", 
              style="quantile", 
              title="Corrected weighted endemism",
              lwd = 0.5,
              palette="YlGnBu") + # palette="YlGnBu"
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1) #+
  # tm_shape(fire_sf)+
  # tm_polygons(
  #             alpha = 0.5,
  #             legend.show = T,
  #             lwd = 2, border.col = 1)


point.maps <- tmap_arrange(spRh, weighted.endemism, Corrected.weighted.endemism)




## Saving plot and results ####

# tmap_save(point.maps, filename = "output/tmap_invert_point_endemism_0-5.png")
# 
# tmap_save(Corrected.weighted.endemism, filename = "output/tmap_invert_point_Wendemism_0-5.png")

### Save out endemism file
save(m1, file = "data/m1_point_0-5.rds")


## END ####
