library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(parallel)
library(biwavelet)
library(geodist)
source('Scripts/functions.R', echo=FALSE)

### map 
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
lbox <- st_bbox(c(xmin = -30, xmax = -90, ymin = -57, ymax = 13), crs = 4326) %>% st_as_sfc()
r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
inBBB <- suppressMessages(unlist(st_intersects(lbox, pol)))

r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(r0))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(r0)


## batches
system("mount_smbfs //slisovsk:X3cp-8PO@smb.isipd.dmawi.de/projects/bioing/user/slisovsk /Users/slisovski/Documents/bioing/slisovsk")

batch_dir <- "~/Documents/bioing/slisovsk/SeasonalChange/Results/Batches/"
phen_dir  <- "~/Documents/bioing/slisovsk/SeasonalChange/Results/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("Results/dateSequence.rda") 

batch = 1027
