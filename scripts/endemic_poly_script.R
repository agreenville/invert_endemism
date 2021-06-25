##****************************************************************************##
# Calc weighed endemism for invertebrates shapefiles
# Aaron Greenville
##****************************************************************************##

## >> load packages ####
library(phyloregion)
library(rgdal)
library(raster)


## >> load data ####

aust = readOGR("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\GIS_files\\World boundaries\\Aust\\AustraliaAdmin.shp")
plot(aust)

## Species with polygons saved a list of Spatial Polygons Data Frames
output_dir <- "data/"
species_polys <- readRDS(file.path(output_dir, "species_ahullEOOspdf.rds"))

# Combine the list of sf objects into a single object
combinedShp <- do.call(what = rbind, args=species_polys)

# add species names
## Associated species names
all_species <- names(species_polys)
combinedShp$dummy <-all_species

##  calc community matrix ####
comm.poly <- polys2comm(dat = combinedShp, species = "dummy", trace=1, res = 0.5)
head(comm.poly)


## calc weighed endemism ####
Endm.invert.poly <- weighted_endemism(comm.poly$comm_dat)
head(Endm.invert.poly)

# join results back to spatial community data
m1.poly <- merge(comm.poly$poly_shp, data.frame(grids=names(Endm.invert.poly), WE=Endm.invert.poly), by="grids")
m1.poly <- m1.poly[!is.na(m1.poly@data$WE),]

## calc corrected weighted endemism #### 
# (weighted endemism tally per cell divided by the species richness of that cell)

m1.poly$corrected_endemism <- m1.poly$WE/m1.poly$richness



## Plotting ####

# save plot
# png(file="output/endemism_inverts_poly_update-0-5.png",
#    width=1000, height=600)

# plot  c(bottom, left, top, rightt)
par(mfrow = c(2, 2), mar=c(0, 0, 0, 0))
  plot_swatch(m1.poly, values = m1.poly$WE, k=10, leg = 5 ,border = NA,key_label = "Weighted endemism",
            pos = "bottomleft")#


  #change datum for aust to match inverts and add
  aust_WGS84 <- spTransform(aust,
                          crs(comm.poly$poly_shp))

  plot(aust_WGS84, add=TRUE)

  #species richness
  plot_swatch(m1.poly, values = m1.poly$richness, k=10, leg = 5 ,border = NA, key_label = "Species richness")
  plot(aust_WGS84, add=TRUE)

  #corrected weighted endemism
  plot_swatch(m1.poly, values = m1.poly$corrected_endemism, k=10, leg = 5 ,border = NA, key_label = "Corrected weighted endemism")
  plot(aust_WGS84, add=TRUE)
par(mfrow = c(1, 1), mar=c(5, 4, 4, 2) + 0.1)

#dev.off()


## Nicer plots ###

library(tmap)

spRh.poly <- tm_shape(m1.poly) +
  tm_polygons("richness", 
              style="quantile", 
              title="Species Richness",
              palette="YlGnBu")+
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)

weighted.endemism.poly <- tm_shape(m1.poly) +
  tm_polygons("WE", 
              style="quantile", 
              title="Weighted endemism",
              palette="YlGnBu")+
  tm_shape(aust_WGS84)+
  tm_polygons("PLACENAME", 
              alpha = 0,
              legend.show = F,
              lwd = 2, border.col = 1)

Corrected.weighted.endemism.poly <- tm_shape(m1.poly) +
  tm_polygons("corrected_endemism", 
              style="quantile", 
              title="Polygons: Corrected weighted endemism",
              lwd = 0.5,
              palette="YlGnBu")+
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




maps.poly <- tmap_arrange(spRh.poly, weighted.endemism.poly, Corrected.weighted.endemism.poly, ncol =1 )

# tmap_save(maps.poly, filename = "output/tmap_invert_poly_endemism_1.png",
#           height = 7, width = 7)

maps.poly.point <- tmap_arrange( Corrected.weighted.endemism, Corrected.weighted.endemism.poly, ncol =1 )

# tmap_save(maps.poly.point, filename = "output/tmap_invert_poly-point_endemism_0-5.png",
#           height = 7, width = 7)


## Save out endemism file ####
save(m1.poly, file = "data/m1_poly_0.5.rds")



