################################################################################
# Calc weigthed endemism for invertebrates shapfiles
# Aaron Greenville
################################################################################

library(phyloregion)
library(rgdal)
library(raster)


#### load data ####

aust = readOGR("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\GIS_files\\World boundaries\\Aust\\AustraliaAdmin.shp")
plot(aust)

# invert shapefiles
baseDir <- file.path("C:\\Users\\aaron\\Dropbox (Sydney Uni)\\Projects\\2019-20_fires\\nesp_fire_invertebrates\\polygons_Aaron\\shapefiles")
filenames <- list.files(baseDir, pattern = ".shp")
filepaths <- paste(baseDir,"/", filenames, sep='')

# Read each shapefile and return a list of sf objects
listOfShp <- lapply(filepaths, readOGR, verbose = FALSE)

# Look to make sure they're all in the same CRS
unique(sapply(listOfShp, crs))

# Combine the list of sf objects into a single object
combinedShp <- do.call(what = rbind, args=listOfShp)

# add species names
combinedShp$dummy <-filenames

# calc community matrix
comm <- polys2comm(dat = combinedShp, species = "dummy", trace=1, res = 0.5)
head(comm)

# calc weigthed endemism
Endm.invert <- weighted_endemism(comm$comm_dat)
head(Endm.invert)

# join results back to spatial community data
m1 <- merge(comm$poly_shp, data.frame(grids=names(Endm.invert), WE=Endm.invert), by="grids")
m1 <- m1[!is.na(m1@data$WE),]

# calc corrected weighted endemism, 
# (weighted endemism tally per cell divided by the species richness of that cell)

m1$corrected_endemism <- m1$WE/m1$richness


# save plot
# png(file="output/endemism_inverts_5.png",
#    width=1000, height=600)

# plot  c(bottom, left, top, rightt)
par(mfrow = c(2, 2), mar=c(0, 0, 0, 0))
plot_swatch(m1, values = m1$WE, k=10, leg = 5 ,border = NA,key_label = "Weighted endemism",
            pos = "bottomleft")#


#change datum for aust to match inverts and add
aust_WGS84 <- spTransform(aust,
                          crs(comm$poly_shp))

plot(aust_WGS84, add=TRUE)

#species richness
plot_swatch(m1, values = m1$richness, k=10, leg = 5 ,border = NA, key_label = "Species richness")
plot(aust_WGS84, add=TRUE)

#corrected weighted endemism
plot_swatch(m1, values = m1$corrected_endemism, k=10, leg = 5 ,border = NA, key_label = "Corrected weighted endemism")
plot(aust_WGS84, add=TRUE)


par(mfrow = c(1, 1), mar=c(5, 4, 4, 2) + 0.1)

#dev.off()
