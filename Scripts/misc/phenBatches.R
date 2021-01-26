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
lbox <- st_bbox(c(xmin = -20, xmax = 50, ymin = 20, ymax = 75), crs = 4326) %>% st_as_sfc()
r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
inEurope <- suppressMessages(unlist(st_intersects(lbox, pol)))


r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(r0))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(r0)


## batches
batch_dir <- "/Volumes/bioing/user/slisovsk/SeasonalChange/Results/Batches/"
phen_dir  <- "/Volumes/bioing/user/slisovsk/SeasonalChange/Results/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("Results/dateSequence.rda") 


## run
# plot(land$geometry, col = adjustcolor("grey90", alpha.f = 0.5))
# plot(pol, add = T)
# plot(lake$geometry, add =T)
# centr <- do.call("rbind", lapply(1:length(pol), function(x) coordinates(gCentroid(pol[x,]))))
# text(centr[batchID,1], centr[batchID,2], batchID, cex = 0.4)


for(batch in batchID[sample(1:length(batchID))]) {

  # batch = 1060
  
if(!file.exists(paste0(phen_dir, "phenBatch_", batch, ".rda"))) { 
  
load(paste0(batch_dir, "Batch_", batch, ".rda"))

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

distM <- suppressMessages(geodist(outBatch$crds[inBatch,1:2], measure = "cheap"))

pxlPhen <- mclapply(which(onLand), function(pxl) {
  
  # pxl <- which(onLand)[3601]
  
  inbfr  <- which(inBatch)[which(distM[which(inBatch)==pxl,]<15000)]

  # points(outBatch$crds[pxl, 1], outBatch$crds[pxl, 2])
  # points(outBatch$crds[inbfr, 1:2], pch = 16, cex = 0.5, col = "magenta")

  dst    <- distM[which(inBatch)==pxl, which(inBatch)%in%inbfr]/1000
  weigth <- approx(c(0, max(dnorm(dst, 0, 4))), c(0,1), dnorm(dst, 0, 4))$y

  dat    <- outBatch$dat[inbfr,,]
  if(outBatch$crds[pxl, 2]>=0) {
    evi    <- ifelse(!is.na(dat[,,2]) | dat[,,1] < -0.1 , NA, dat[,,1])
  }  else evi    <- ifelse(dat[,,1] < -0.1 , NA, dat[,,1])

  # matplot(dates, t(evi), type= "o", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16)
  
  med <- apply(evi, 2, median , na.rm = T)
  
  wt <- biwavelet::wt(cbind(1:length(med), zoo::na.approx(med)))
  power  <- log2(wt$power.corr)
  per    <- wt$period[which.max(apply(power, 1, median, na.rm = T))]

  if(round(per,0)%in%c(42:62)) {
      period = 1
  } else {
      if(round(wt$period[which.max(apply(power, 1, median))],0)%in%c(16:36)) {
        period = 0.5
    } else period = NA
  }
    
  if((diff(quantile(med, prob = c(0.025,0.975), na.rm = T))>0.15) & !is.na(period)) { ## diff ts > 0.1 OR period == NA (no annual or bi-annual seasonality)
    
    datCurve <- merge(data.frame(year = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%Y")),
                                 week = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%U"))),
                      data.frame(year = as.numeric(format(dates, "%Y")), week = as.numeric(format(dates, "%U")),
                                 date = dates,
                                 evi  = med, t(evi)), all.x = T)
      
    if(period == 1) {
      fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = 52, Mx = zoo::na.approx(datCurve$evi) - median(datCurve$evi, na.rm = T), sd = 0.01)
      curve <- fit0$par[1]*cos(pi*((1:length(datCurve$evi))/(length(datCurve$evi)/((length(datCurve$evi)/52)*2))) +
                                 (pi+fit0$par[2])) +  mean(med, na.rm=T)
    } else{
      fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = 26, Mx = zoo::na.approx(datCurve$evi) - median(datCurve$evi, na.rm = T), sd = 0.01)
      curve <- fit0$par[1]*cos(pi*((1:length(datCurve$evi))/(length(datCurve$evi)/((length(datCurve$evi)/26)*2))) +
                                 (pi+fit0$par[2])) +  mean(med, na.rm=T)
    }
    
   # plot( datCurve$evi, type = "o", pch = 16, cex = 0.5, col = "darkgreen")
   # par(new = T)
   # plot(curve, type = "l", col = "orange")
    
    mins <- unique(c(1, which(diff(sign(diff(-curve)))==-2)+1, length(curve)))
    maxs <- which(diff(sign(diff( curve)))==-2)+1
    
    segL <- apply(do.call("rbind", lapply(maxs, function(x) mins[c(1:length(mins))[order(abs(x-mins))][1:2]]+c(-10,10))), 1, 
                  function(x) datCurve[ifelse(x[1]<1, 1, x[1]):ifelse(x[2]>length(curve), length(curve), x[2]),])
    
    q50 <- median(sapply(segL, function(x) quantile(x[,4], prob = 0.5, na.rm = T)), na.rm = T)
    
    phen0 <- do.call("rbind", lapply(segL, function(s) {
      
      # s <- segL[[77]]
      # matplot(s$date, s[,-c(1:4)], type= "l", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16)
      # lines(s$date, s[,-c(1:4)][,weigth==1], lwd = 4, col = "orange")

      dateSeg <- zoo::na.approx(s$date)
      year    <- median(as.numeric(format(s$date, "%Y")), na.rm = T)
      
      if(length(dateSeg)>20) {
      
        
      tmpDat <- subset(data.frame(x = rep(s$date, ncol(s)-4),
                                  y = unlist(t(c(s[,-c(1:4)]))),
                                  w = rep(weigth, each = nrow(s))), !is.na(y))
      
      spl     <- smooth.spline(x = tmpDat$x, y = tmpDat$y, spar = 0.3, w = tmpDat$w)
      xSmooth <- predict(spl, dateSeg)$y
      peaks   <- FindPeaks(xSmooth)
      
      # lines(dateSeg, xSmooth, lwd = 6, col = "orange")
      # abline(v = dateSeg[peaks], lty = 2, col = "grey80")
      
      pars    <- list(rel_amp_frac = 0.15, rel_peak_frac = 0.1, min_seg_amplitude = 0.1)
      segs    <- tryCatch(do.call("rbind", GetSegs(peaks, xSmooth, pars)), error = function(e) NULL)
      
      # apply(segs, 1, function(x) rect(dateSeg[x[1]], 0, dateSeg[x[3]], 1, col = adjustcolor("grey80", alpha.f = 0.2), border =  "grey10"))
      
      if(!is.null(segs)) {      
                 
            if(nrow(segs)>1) {
              
              max_in  <- min(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
              max_out <- max(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
              
              seqRan  <- c(min(segs[order(xSmooth[segs[,2]], decreasing = T),1]),
                           max(segs[order(xSmooth[segs[,2]], decreasing = T),3]))
                
            } else {
              
              max_in <- max_out  <- segs[2]
              seqRan <- segs[c(1,3)]
              
            }
        
            max    <- as.numeric(format(as.POSIXct(dateSeg, origin = "1970-01-01")[max_in], "%j"))
            amp    <- diff(range(xSmooth[seqRan[1]:seqRan[2]]))
            
            
            q10gup <- with(data.frame(t = dateSeg, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
              # abline(v = q10gup)
            q50gup <- with(data.frame(t = dateSeg, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(q50, rev(y), gup = F, first_greater = T)])
              # abline(v = q50gup)
            
            q10sen <- with(data.frame(t = dateSeg, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
              # abline(v = q10sen)
            q50sen <- with(data.frame(t = dateSeg, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(q50, rev(y), gup = T, first_greater = F)])
              # abline(v = q50sen)
            
            area   <- with(data.frame(t = 1:length(dateSeg), y = xSmooth)[seqRan[1]:seqRan[2],],
                           MESS:::auc(t, y, type = 'spline'))
            
            
            out <- c(year = year, max = max, amp = round(amp, 2), area = round(area,2), 
                     as.numeric(format(as.POSIXct(c(q10gup, q50gup, q10sen, q50sen), origin = "1970-01-01"), "%j")))
              out[5:6] <- ifelse(out[5:6]>=max, out[5:6] - 365, out[5:6])
              out[7:8] <- ifelse(out[7:8]<=max, 365 + out[7:8], out[7:8])
            
            out
      } else {
        c(year = year, max = NA, amp = NA, area = NA, rep(NA, 4))
      }
      
      } else {
        c(year = year, max = NA, amp = NA, area = NA, rep(NA, 4))
      }
    }))
    
    if(period==1) {
      phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
    } else {
      pk   <- round(median(unlist(aggregate(phen0[,3], by = list(phen0[,1]), FUN = function(a) which.max(a))$x), na.rm = T), 0)
      phen <- do.call("rbind", lapply(split(as.data.frame(phen0), round(phen0[,1],0)), function(a) {
        a[1,]
      }))
    }

    merge(data.frame(year = 1981:2020), data.frame(year = phen[,1],
                                                   per  = round(per,0),
                                                   max  = phen[,2],
                                                   amp  = phen[,3],
                                                   area = phen[,4],
                                                   q10gup = phen[,5],
                                                   q50gub = phen[,6],
                                                   q10sen = phen[,7],
                                                   q50sen = phen[,8]), all.x = T)[,-1]
    
    
    
  } else { #### diff ts range < 0.1
    
    data.frame(year = 1981:2020, 
               per = NA,
               max  = NA, 
               amp  = NA, area = NA, 
               q10gup = NA, q50gub = NA, 
               q10sen = NA, q50sen = NA)[,-1]
    
  }
  
}, mc.cores = max(5, detectCores()-10))

## period, max, amp, area, q10gub, q50gub, q10sen, q50sen
phenA <- tryCatch(abind::abind(lapply(1:ncol(pxlPhen[[1]]), function(l) matrix(unlist(sapply(pxlPhen, function(x) x[,l])), ncol = nrow(pxlPhen[[1]]), byrow = T)), along = 3),
                  error = function(e) NULL)

phenR <- rasterFromXYZ(cbind(outBatch$crds[outBatch$crds[,3],1:2], 
                             (apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
                                !outBatch$crds[,ncol(outBatch$crds)])[outBatch$crds[,3]]), crs = CRS("+proj=longlat"))
names(phenR) <- paste0("batch_", batch)

phenOut <- list(dat = phenA, raster = phenR)

save(phenOut, file = paste0(phen_dir, "phenBatch_", batch, ".rda"))

# rTest <- phenOut$raster
# rTest[rTest[]==1] <- apply(phenOut$dat[,,1], 1, median, na.rm = T)
# plot(rTest)

} ## file exists

}