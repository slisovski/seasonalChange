library(sf)
library(raster)
library(rgdal)
library(rgeos)

## land mask
lake <- read_sf("Data/GeoDat/ne_50m_lakes/ne_50m_lakes.shp")
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp")
lbox <- st_bbox(land, crs = 4326) %>% st_as_sfc()

r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(r0))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(r0)


## batches
batch_dir <- "/Users/slisovsk/Desktop/batches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("Results/dateSequence.rda") 


## run
plot(land$geometry, col = adjustcolor("grey90", alpha.f = 0.5))
plot(pol, add = T)
# plot(lake$geometry, add =T)

centr <- do.call("rbind", lapply(1:length(pol), function(x) coordinates(gCentroid(pol[x,]))))
text(centr[batchID,1], centr[batchID,2], batchID, cex = 0.4)

# for(batch in batchID) {
  
  batch <- 613

  load(paste0(batch_dir, "batch_", batch, ".rda"))
  
  onLand  <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) & 
    !outBatch$crds[,ncol(outBatch$crds)] & 
    outBatch$crds[,3]
  
  inBatch <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) & 
    !outBatch$crds[,ncol(outBatch$crds)]
  
  st_crds  <- st_sfc(lapply(1:nrow(outBatch$crds), function(x) st_point(as.numeric(outBatch$crds[x,1:2]))), crs = 4326)
  
  # test <- rasterFromXYZ(cbind(outBatch$crds[,1:2], outBatch$dat[,420,1]))
  # plot(test, breaks = seq(0,0.8, length = 100), col = terrain.colors(100))
  # plot(land$geometry, add = T, col = NA)
  # plot(lake$geometry, add = T)
  # points(outBatch$crds[inBatch,1:2], cex = 0.3, pch = 16, col = "black")
  # points(outBatch$crds[onLand,1:2], cex = 0.3, pch = 16, col = "red")

  coordinates_dt <- optiRum::CJ.dt(data.table::as.data.table(outBatch$crds[inBatch,1:2]),
                                   data.table::as.data.table(outBatch$crds[inBatch,1:2]))
  
  distM <- matrix(with(coordinates_dt, spatialrisk::haversine(y, x, i.y, i.x)), 
                  ncol = sum(inBatch), nrow = sum(inBatch))

  
  pxlPhen <- mclapply(which(onLand), function(pxl) {
    
      inbfr  <- which(inBatch)[which(distM[which(inBatch)==pxl,]<15000)]
      
      # points(outBatch$crds[pxl, 1], outBatch$crds[pxl, 2])
      # points(outBatch$crds[inbfr, 1:2], pch = 16, cex = 0.5, col = "magenta")
      
      dst    <- distM[pxl,inbfr]/1000
      weigth <- approx(c(0, max(dnorm(dst, 0, 5))), c(0,1), dnorm(dst, 0, 5))$y
      
      dat    <- outBatch$dat[inbfr,,]
      evi    <- ifelse(!is.na(dat[,,2]) | dat[,,1] < -0.1 , NA, dat[,,1])
      
  
      ### segments
      med <- apply(evi, 2, median , na.rm = T)
      # matplot(dates, t(evi), type= "o", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16)
      
      if(sum(!is.na(med))>(length(med)*0.75)) {    
          
          q50 <- quantile(evi, prob = 0.5, na.rm = T)
          
          peaks   <- FindPeaks(med)
          pars    <- list(rel_amp_frac = 0.25, rel_peak_frac = 0.1, min_seg_amplitude = 0.15)
          
          segs    <- do.call("rbind", GetSegs(peaks, med, pars))
          segs    <- segs[order(segs[,1]),]
          
          # plot(dates, med, type = "o")
          # apply(segs, 1, function(x) rect(dates[x[1]], 0, dates[x[3]], 1, col = adjustcolor("grey90", alpha.f = 0.2), border =  "grey10"))
          # abline(v = dates[segs[,2]], col = "orange", lty = 2)
          
          segL    <- lapply(1:nrow(segs), function(x) evi[,segs[x,1]:segs[x,3]])  
          dtsS    <- lapply(1:nrow(segs), function(x) dates[segs[x,1]:segs[x,3]])  
          
          phen <- do.call("rbind", lapply(1:nrow(segs), function(s) {
          
                # matplot(dtsS[[s]], t(segL[[s]]), type= "o", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16)
                tmpDat <- data.frame(x = rep(dtsS[[s]], each = nrow(segL[[s]])),
                                     y = c(segL[[s]]),
                                     w = rep(weigth, each = ncol(segL[[s]])))    
            
                spl <- with(tmpDat[!is.na(tmpDat$y),], smooth.spline(x, y, spar = 0.3, w = w))
                xSmooth <- predict(spl, as.numeric(dtsS[[s]]))$y
        
                spline  <- data.frame(ind = 1:length(dtsS[[s]]), d = dtsS[[s]], evi = xSmooth, 
                                      gup = ifelse(dtsS[[s]]<=dates[segs[s,2]], TRUE, FALSE))
                
                # lines(dtsS[[s]], xSmooth, lwd = 2, col = "orange")
                # abline(v = dates[segs[s,2]])
                # abline(h = q50)
                
                year   <- median(as.numeric(format(dtsS[[s]], "%Y")))
                
                max    <- as.numeric(format(dates[segs[s,2]], "%j"))
                maxInd <- which(spline$d==dates[segs[s,2]])
                amp    <- diff(range(xSmooth))
        
                  
                q10gup <- with(spline[spline$gup,],  rev(ind)[GetThresh(min(evi) + diff(range(evi))*0.1, rev(evi), gup = F, first_greater = T)])
                q10sen <- with(spline[!spline$gup,], rev(ind)[GetThresh(min(evi) + diff(range(evi))*0.1, rev(evi), gup = T, first_greater = F)])
                
                q50gup <- with(spline[spline$gup,],  rev(ind)[GetThresh(q50, rev(evi), gup = F, first_greater = T)])
                q50sen <- with(spline[!spline$gup,], rev(ind)[GetThresh(q50, rev(evi), gup = T, first_greater = F)])
                
                area   <- with(spline[spline$ind>=q10gup & spline$ind<=q10sen,],
                               MESS:::auc(ind, evi, type = 'spline'))
                
                c(year, max, amp, area, abs(c(q10gup, q50gup, q10sen, q50sen)-maxInd))

        }))
    
          phenOut <- merge(data.frame(year = 1981:2020), data.frame(year = phen[,1],
                                                                    max  = phen[,2],
                                                                    amp  = phen[,3],
                                                                    area = phen[,4],
                                                                    q10gup = phen[,5],
                                                                    q50gub = phen[,6],
                                                                    q10sen = phen[,7],
                                                                    q50sen = phen[,8]), all.x = T)
  
      }  else {
        
         data.frame(year = 1981:2020, 
                    max  = NA, 
                    amp  = NA, area = NA, 
                    q10gup = NA, q50gub = NA, 
                    q10sen = NA, q50sen = NA)
      }
  

      
  }, mc.cores = 5)
  

mx <- sapply(pxlPhen, function(x) median(x[,7], na.rm = T))
  
test0   <- rasterFromXYZ(cbind(outBatch$crds[,1:2], outBatch$dat[,420,1]))
test0[] <- NA
test0[which(onLand)] <- mx
  
plot(test0)
  