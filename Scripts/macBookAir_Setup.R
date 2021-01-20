library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(zoo)
library(biwavelet)
library(MESS)
library(zoo)
library(parallel)
library(geodist)
source('Scripts/functions.R', echo=FALSE)

### map 
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)

r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- rasterToPolygons(r0)


## batches
# system("mount_smbfs //slisovsk:X3cp-8PO@smb.isipd.dmawi.de/projects/bioing/user/slisovsk /Users/slisovski/Documents/bioing/slisovsk")

batch_dir <- "~/Documents/bioing/slisovsk/Batches/"
phen_dir  <- "~/Documents/bioing/slisovsk/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("Results/dateSequence.rda") 

# plot(land)
# plot(pol[batchID,])

load("Scripts/isNull.rda")
batch <- batchID[batchID%in%isNull][1]
