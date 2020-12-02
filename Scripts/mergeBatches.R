library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(parallel)

## land mask
lake <- read_sf("Data/GeoDat/ne_50m_lakes/ne_50m_lakes.shp")
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp")
lbox <- st_bbox(land, crs = 4326) %>% st_as_sfc()

r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(r0))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(r0)

## batches
phen_dir  <- "/Volumes/bioing/user/slisovsk/SeasonalChange/Results/phenBatches/"
batches   <- list.files(phen_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))

## emptyRaster
r0    <- raster("Results/phenRaster.tif")
r_ind <- r0; r_ind[] <- 1:length(r_ind[])

ampDat <- do.call("rbind", mclapply(batchID, function(x) {

  # plot(r0)
  # plot(pol[x,], add = T)
  
  load(paste0(phen_dir, batches[which(batchID==x)]))
  
  out <- cbind(unlist(raster::extract(r_ind, pol[x,])), NA)
  
  if(!is.null(phenOut$dat)) {
    ## period, max, amp, area, q10gub, q50gub, q10sen, q50sen
    out[phenOut$raster[]==1,2] <- apply(phenOut$dat[,,3], 1, median, na.rm = T)
  } else {
    out[phenOut$raster[]==1,2] <- 0
  }
  
  out
  
}))

ampMed <- r0; ampMed[ampDat[,1]] <- ampDat[,2]

plot(land$geometry, col = "grey90", border = NA)
plot(ampMed, add = T)
