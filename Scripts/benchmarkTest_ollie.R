args  <- commandArgs(trailingOnly = T)

#### libraries
.libPaths(c(.libPaths(), "/home/ollie/slisovsk/RServerLibs/3.6"))

suppressPackageStartupMessages({
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

  ### google drive
  library(googledrive)
  library(httr)
  drive_auth(email = "simeon.lisovski@gmail.com", cache = "/work/ollie/slisovsk/Projects/ArcticPhenology/gargle-oauth")
  set_config(config(ssl_verifypeer = 0L))
  ###  
  
source('/work/ollie/slisovsk/Projects/SeasonalChange/Scripts/functions.R', echo=FALSE)
})


### map 
land <- read_sf("Data/GeoDat/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
lbox <- st_bbox(c(xmin = 108, xmax = 130, ymin = 55, ymax = 80), crs = 4326) %>% st_as_sfc()
r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
inBBB <- suppressMessages(unlist(st_intersects(lbox, pol)))

### batches
batch_dir <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/Batches/"
phen_dir  <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("/work/ollie/slisovsk/Projects/SeasonalChange/Results/dateSequence.rda") 


### Pixel analysis
batch <- batchID[which(batchID%in%inBBB)][as.numeric(args[1])]

if(!file.exists(paste0(phen_dir, "phenBatch_", batch, ".rda"))) {

  load(paste0(batch_dir, "Batch_", batch, ".rda"))
  
  inBatch_sf = st_as_sf(outBatch$crds[,c(1:2)], coords = c("x", "y"), crs = 4326, agr = "constant")
  inBatch_sf$onLand <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
    !outBatch$crds[,ncol(outBatch$crds)] &
    outBatch$crds[,3]
  inBatch_sf$inBatch <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
    !outBatch$crds[,ncol(outBatch$crds)]

  # a <- Sys.time()
  pxlPhen <- mclapply(which(inBatch_sf$onLand), evalPxl, mc.cores = 15)
  # Sys.time() - a
  
  tt <- sapply(pxlPhen, function(x) is.data.frame(x))
  pxlPhen[[which(!tt)[1]]]
  
  phenA <- tryCatch(abind::abind(lapply(1:ncol(pxlPhen[[1]]), function(l) matrix(unlist(sapply(pxlPhen, function(x) x[,l])), ncol = nrow(pxlPhen[[1]]), byrow = T)), along = 3),
                    error = function(e) NULL)
  
  phenR <- rasterFromXYZ(cbind(outBatch$crds[outBatch$crds[,3],1:2],
                               (apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
                                  !outBatch$crds[,ncol(outBatch$crds)])[outBatch$crds[,3]]), crs = CRS("+proj=longlat"))
  
  # r_test <- phenR; r_test[] <- NA
  # r_test[phenR[]==1][1:500] <-  apply(phenA, c(1,3), median, na.rm = T)[,9]
  # plot(r_test)
  
  names(phenR) <- paste0("batch_", batch)
  phenOut <- list(dat = phenA, raster = phenR)
  save(phenOut, file = paste0(phen_dir, "phenBatch_", batch, ".rda"))

  drive_upload(paste0(phen_dir, "phenBatch_", batch, ".rda"),
               path = "SeasonalChange/phenBatchesRDA/", name = paste0("phenBatch_", batch, ".rda"), verbose = FALSE, overwrite = T)
}