library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(parallel)

## land mask
land   <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
# lbox   <- st_bbox(c(xmin = -20, xmax = 50, ymin = 20, ymax = 75), crs = 4326) %>% st_as_sfc()
r0     <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol    <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
# inEurope <- suppressMessages(unlist(st_intersects(lbox, pol)))


## batches
phen_dir  <- "/Users/slisovsk/Google Drive/SeasonalChange/phenBatchesRDA/"
batches   <- list.files(phen_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
# euBatch   <- batchID[which(batchID%in%inEurope)]

## emptyRaster
r0    <- raster("Results/phenRaster.tif")
r_ind <- r0; r_ind[] <- 1:length(r_ind[])

medDat <- do.call("rbind", mclapply(batchID, function(x) {
  
  load(paste0(phen_dir, batches[which(batchID==x)]))
  
  ind <- unlist(raster::extract(r_ind, as(pol[x,], "Spatial")))
  out <- cbind(ind, matrix(ncol = 6, nrow = length(ind)))
  
  if(!is.null(phenOut$dat)) {
    out[phenOut$raster[]==1,2:7] <- apply(phenOut$dat, c(1,3), median, na.rm = T)
  }

  out

}, mc.cores = 5, mc.preschedule = T))


datR <- do.call("stack", mclapply(2:7, function(x) {
  rOut <- r0; rOut[medDat[,1]] <- medDat[,x]
  rOut
}, mc.cores = 5))


lbox <- st_bbox(c(xmin = -30, xmax = -90, ymin = -57, ymax = 13), crs = 4326) %>% st_as_sfc()
datRC <- crop(datR, as(lbox, "Spatial"))
landC <- land %>% st_intersection(lbox)

# 1: per1, 2: sig1, 3: per2, 4: sig2, 5: per3, 6: sig3, 7: max, 8: amp, 9: area, 10: st1, 11: st2, 12: en2, 13: en1 
plot(landC, col = "grey90", border = NA)
plot(datRC[[6]], add = T)

plot(pol[batchID], add = T)
# crds <- st_coordinates(st_centroid(pol))
text(crds[,1], crds[,2], 1:nrow(crds), cex = 0.5)
