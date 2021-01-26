library(doParallel)
registerDoParallel(cl <- makeCluster(4))

### map 
land   <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
lbox   <- st_bbox(c(xmin = -20, xmax = 50, ymin = 20, ymax = 75), crs = 4326) %>% st_as_sfc()
r0     <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol    <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
inEurope <- suppressMessages(unlist(st_intersects(lbox, pol)))


## batches
phen_dir  <- "/Users/slisovsk/Google Drive/SeasonalChange/phenBatchesRDA/"
batches   <- list.files(phen_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
euBatch   <- batchID[which(batchID%in%inEurope)]


r0    <- raster("Results/phenRaster.tif")
r_ind <- r0; r_ind[] <- 1:length(r_ind[])


for(b in euBatch) {
  
  b <- 479
  load(paste0(phen_dir, "phenBatch_", b, ".rda"))
  
  
  res <- foreach(j = seq_len(dim(phenOut$dat)[1]), .combine = "rbind") %:% 
    foreach(i = seq_len(dim(phenOut$dat)[3]), .combine = "c") %dopar% {
      dat <- phenOut$dat[j,,i]
      if(sum(is.na(dat))<length(dat)*0.3) {
        mod <- lm(dat~c(1981:2020))
        summary(mod)$coefficients[2,1]
      } else NA
    }
    
  rTest      <- phenOut$raster
  ind        <- phenOut$raster[]==1
  rTest[ind] <- res[,5]
  plot(rTest)
  hist(rTest, breaks = 30)
  
  
  
    
stopCluster(cl)
  
  
  
  
  
  
  
  
  rTest      <- phenOut$raster
  ind        <- phenOut$raster[]==1
  rTest[ind] <- ch
  plot(rTest)
  hist(rTest, breaks = 30)
}