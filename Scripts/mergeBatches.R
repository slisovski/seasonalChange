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
# phen_dir  <- "/Users/slisovski/Google Drive/SeasonalChange/phenBatchesRDA/"
phen_dir  <- "/Users/slisovsk/Google Drive/SeasonalChange/phenBatchesRDA/"
batches   <- list.files(phen_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
# lbox <- st_bbox(c(xmin = 80, xmax = 170, ymin = 50, ymax = 80), crs = 4326) %>% st_as_sfc()
# inBBB <- suppressMessages(unlist(st_intersects(lbox, pol)))


isNull <- unlist(mclapply(batchID, function(x) {
  load(paste0(phen_dir, batches[which(batchID==x)]))
  if(is.null(phenOut$dat)) x 
  }, mc.cores = 4))


save(isNull, file = "Scripts//isNull.rda")


plot(land, col = "grey60", border = NA)
plot(pol[batchID[batchID%in%isNull],], col = adjustcolor("orange", 0.1), add = T)
crds <- suppressWarnings(st_coordinates(st_centroid(pol)))
text(crds[batchID[batchID%in%isNull],1], crds[batchID[batchID%in%isNull],2], isNull, cex = 0.5)


# plot(landC, col = "grey90", border = NA)
# plot(pol[which(batchID%in%batchID[inBBB])], add = T, col = adjustcolor("transparent", alpha.f = 0.5))

# plot(pol[batchID,])
# plot(pol[batchID[batchID%in%inBBB],], col = "orange", add = T)

## emptyRaster
r0    <- raster("Results/phenRaster.tif")
r_ind <- r0; r_ind[] <- 1:length(r_ind[])

medDat <- do.call("rbind", mclapply(batchID[batchID%in%isNull], function(x) {
  
  load(paste0(phen_dir, batches[which(batchID==x)]))
  
  # r_test <- phenOut$raster; r_test[] <- NA
  # r_test[phenOut$raster[]==1][1:500] <-  apply(phenOut$dat, c(1,3), median, na.rm = T)[,9]
  # plot(r_test)
  
  ind <- unlist(raster::extract(r_ind, as(pol[x,], "Spatial")))
  out <- cbind(ind, matrix(ncol = 13, nrow = length(ind)))
  
  if(!is.null(phenOut$dat)) {
      out[which(phenOut$raster[]==1),2:14] <- apply(phenOut$dat, c(1,3), median, na.rm = T)
  }

  out

}, mc.cores = 4))


datR <- do.call("stack", mclapply(c(3,8), function(x) {
  rOut <- r0; rOut[medDat[,1]] <- medDat[,x]
  rOut
}, mc.cores = 4))
# plot(datR[[1]])


# datRC <- crop(datR, as(lbox %>% st_buffer(5), "Spatial"))
# landC <- land %>% st_intersection(lbox %>% st_buffer(5))

# 1: per1, 2: sig1, 3: per2, 4: sig2, 5: per3, 6: sig3, 7: max, 8: amp, 9: area, 10: st1, 11: st2, 12: en2, 13: en1 
plot(land, col = "grey60", border = NA)
plot(pol[batchID[batchID%in%isNull],], col = adjustcolor("orange", 0.1), add = T)

crds <- st_coordinates(st_centroid(pol))
text(crds[batchID,1], crds[batchID,2], batchID, cex = 0.5)
plot(datR[[2]], add = T, legend = T)

