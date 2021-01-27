library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(abind)
library(parallel)


## land mask
land   <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
r0     <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol    <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)


## batches
# phen_dir  <- "/Users/slisovski/Google Drive/SeasonalChange/phenBatchesRDA/"
phen_dir  <- "/Users/slisovski/Google Drive/SeasonalChange/phenBatchesRDA/"
batches   <- list.files(phen_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))


# isNull <- unlist(mclapply(batchID, function(x) {
#   load(paste0(phen_dir, batches[which(batchID==x)]))
#   if(is.null(phenOut$dat)) x 
#   }, mc.cores = 4))
# 
# save(isNull, file = "Scripts//isNull.rda")

## emptyRaster
r0    <- raster("Results/phenRaster.tif")
r_ind <- r0; r_ind[] <- 1:length(r_ind[])

for(i in 1:13) {
  
  medDat <- do.call("rbind", lapply(batchID, function(x) {
  
      load(paste0(phen_dir, batches[which(batchID==x)]))
      
      # r_test <- phenOut$raster; r_test[] <- NA
      # r_test[phenOut$raster[]==1] <-  apply(phenOut$dat, c(1,3), median, na.rm = T)[,9]
      # plot(r_test)
      
      ind <- unlist(raster::extract(r_ind, as(pol[x,], "Spatial")))
    
      if(!is.null(phenOut$dat) & length(ind[which(phenOut$raster[]==1)])>1) {
        cbind(ind[which(phenOut$raster[]==1)], phenOut$dat[,,i])
      }

  }))

  
  rOut <- do.call("stack", lapply(2:41, function(y) {
    rOut <- r0; rOut[medDat[,1]] <- medDat[,y]
    rOut
  }))
  
  names(rOut) <- paste0("Y", 1981:2020)
  save(rOut, file = paste0('~/Desktop/tmp/', "_", i, 'test.rda'))
  
}
  

listfls <- list.files("/Volumes/projects/bioing/user/slisovsk/", pattern = "*raster.rda")
 id     <- as.numeric(sapply(strsplit(listfls, "_"), function(x) x[[1]]))

seasList <- lapply(order(id), function(x) {
  load(paste0("/Volumes/projects/bioing/user/slisovsk/", listfls[x]))
  rOut
}) 
 

nms <- c("per1", "perSig1", "per2", "perSig2", "seasMax", "seasAmp", "seasArea", "gup10", "gup50", "gup90", "sen90", "sen50", "sen10")
length(seasList)

names(seasList) <- nms
save(seasList, file = "~/Google Drive/SeasonalChange/seasList.rda")


  

# for(y in 2:41) {
#     rOut <- r0; rOut[medDat[,1,i]] <- medDat[,y,i]
#     writeRaster(rOut,paste0('~/Desktop/tmp/', i, "_", y, 'test.tif'), options=c('TFW=YES'))
# }
# 
# }
# 
# 
# seasonalList <- lapply(1:13, function(x) {
#   
#   
#   
#   
#   rOut <- do.call("stack", lapply(2:41, function(y) {
#     rOut <- r0; rOut[medDat[,1,1]] <- medDat[,y,x]
#     rOut
#   }))
#   names(rOut) <- paste0("Y", 1981:2020)
#   rOut
#   
# })
# 
# 
# 
# 
# datR <- do.call("stack", mclapply(c(2:8), function(x) {
#   rOut <- r0; rOut[medDat[,1]] <- medDat[,x]
#   rOut
# }, mc.cores = 4))
# # plot(datR[[1]])
# 
# 
# # lbox <- st_bbox(c(xmin = 108, xmax = 135, ymin = 63, ymax = 75), crs = 4326) %>% st_as_sfc()
# # 
# # 
# # datRC <- crop(datR, as(lbox, "Spatial"))
# # landC <- land %>% st_intersection(lbox)
# # 
# # plot(landC, col = "grey60", border = NA)
# # plot(datRC[[5]], add = T, legend = T)
# 
# landP <- land %>% st_transform(CRS("+proj=moll"))
# datRP <- projectRaster(datR, crs = CRS("+proj=moll"))
# 
# # 1: per1, 2: sig1, 3: per2, 4: sig2, 5: per3, 6: sig3, 7: max, 8: amp, 9: area, 10: st1, 11: st2, 12: en2, 13: en1 
# plot(land, col = "grey60", border = NA)
# plot(pol[batchID,], col = adjustcolor("orange", 0.1), add = T)
# 
# crds <- st_coordinates(st_centroid(pol))
# text(crds[batchID,1], crds[batchID,2], batchID, cex = 0.5)
# plot(datR[[1]], add = T, legend = T, col = rev(viridis::inferno(100)))

