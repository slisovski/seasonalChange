args  <- commandArgs(trailingOnly = T)

.libPaths(c(.libPaths(), "/home/ollie/slisovsk/RServerLibs/3.6"))

library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(zoo)
library(biwavelet)
library(parallel)
library(MESS)
library(zoo)
library(geodist)

source('/work/ollie/slisovsk/Projects/SeasonalChange/Scripts/functions.R', echo=FALSE)

### google drive
library(googledrive)
library(httr)
drive_auth(email = "simeon.lisovski@gmail.com", cache = "/work/ollie/slisovsk/Projects/ArcticPhenology/gargle-oauth")
set_config(config(ssl_verifypeer = 0L))
###


## batches
batch_dir <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/Batches/"
phen_dir  <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("/work/ollie/slisovsk/Projects/SeasonalChange/Results/dateSequence.rda") 

batch <- batchID[as.numeric(args[1])]

load(paste0(batch_dir, "Batch_", batch, ".rda"))

onLand  <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
  !outBatch$crds[,ncol(outBatch$crds)] &
  outBatch$crds[,3]

inBatch <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
  !outBatch$crds[,ncol(outBatch$crds)]

distM   <- suppressMessages(geodist(outBatch$crds[inBatch,1:2], measure = "cheap"))

pxlPhen <- mclapply(which(onLand), function(pxl) {
  
  inbfr  <- which(inBatch)[which(distM[which(inBatch)==pxl,]<15000)]
  dst    <- distM[which(inBatch)==pxl, which(inBatch)%in%inbfr]/1000
  weigth <- approx(c(0, max(dnorm(dst, 0, 4))), c(0,1), dnorm(dst, 0, 4))$y
  
  dat    <- outBatch$dat[inbfr,,]
  if(outBatch$crds[pxl, 2]>=0) {
    evi    <- ifelse(!is.na(dat[,,2]) | dat[,,1] < -0.1 , NA, dat[,,1])
  }  else evi    <- ifelse(dat[,,1] < -0.1 , NA, dat[,,1])
  
  med <- apply(evi, 2, median , na.rm = T)
  
  if(sum(is.na(med))<length(med)*0.33) {
    
    wt     <- wt(cbind(1:length(med), na.approx(med)))
    power  <- log2(wt$power.corr)
    per    <- wt$period[which.max(apply(power, 1, median, na.rm = T))]
    pow    <- wt$power[which.max(apply(power, 1, median, na.rm = T))]
    
    
    if(round(per,0)%in%c(42:62)) {
      period = 1
    } else {
      if(round(wt$period[which.max(apply(power, 1, median))],0)%in%c(16:36)) {
        period = 0.5
      } else period = NA
    }
    
      
    datCurve <- merge(data.frame(year = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%Y")),
                                 week = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%U"))),
                      data.frame(year = as.numeric(format(dates, "%Y")), week = as.numeric(format(dates, "%U")),
                                 date = dates,
                                 evi  = med, pow = apply(power, 2, median, na.rm = T), t(evi)), all.x = T)
      
    if(period == 1) {
        fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = 52, Mx = na.approx(datCurve$evi) - median(datCurve$evi, na.rm = T), sd = 0.01)
        curve <- fit0$par[1]*cos(pi*((1:length(datCurve$evi))/(length(datCurve$evi)/((length(datCurve$evi)/52)*2))) +
                                   (pi+fit0$par[2])) +  mean(med, na.rm=T)
    } else {
        fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = 26, Mx = na.approx(datCurve$evi) - median(datCurve$evi, na.rm = T), sd = 0.01)
        curve <- fit0$par[1]*cos(pi*((1:length(datCurve$evi))/(length(datCurve$evi)/((length(datCurve$evi)/26)*2))) +
                                   (pi+fit0$par[2])) +  mean(med, na.rm=T)
    }
      
    mins <- unique(c(1, which(diff(sign(diff(-curve)))==-2)+1, length(curve)))
    maxs <- which(diff(sign(diff( curve)))==-2)+1
      
      
    if((diff(quantile(med, prob = c(0.025,0.975), na.rm = T))>0.15) & !is.na(period) & length(mins)>20 & length(maxs)>20) {
      
          segL <- apply(do.call("rbind", lapply(maxs, function(x) mins[c(1:length(mins))[order(abs(x-mins))][1:2]]+c(-10,10))), 1,
                        function(x) datCurve[ifelse(x[1]<1, 1, x[1]):ifelse(x[2]>length(curve), length(curve), x[2]),])
          
          q50 <- median(sapply(segL, function(x) quantile(x[,4], prob = 0.5, na.rm = T)), na.rm = T)
          
          phen0 <- do.call("rbind", lapply(segL, function(s) {
            
            # s <- segL[[1]]
            # matplot(s$date, s[,-c(1:5)], type= "l", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16)
            # lines(s$date, s[,-c(1:5)][,weigth==1], lwd = 4, col = "orange")
            
            dateSeg <- na.approx(s$date)
            year    <- median(as.numeric(format(s$date, "%Y")), na.rm = T)
            
            if(length(dateSeg)>20) {
              
              tmpDat <- subset(data.frame(x = rep(s$date, ncol(s)-5),
                                          y = unlist(t(c(s[,-c(1:5)]))),
                                          w = rep(weigth, each = nrow(s))), !is.na(y))
              
              spl     <- smooth.spline(x = tmpDat$x, y = tmpDat$y, spar = 0.3, w = tmpDat$w)
              xSmooth <- predict(spl, dateSeg)$y
              peaks   <- FindPeaks(xSmooth)
              
              # lines(dateSeg, xSmooth, lwd = 6, col = "orange")
              # abline(v = dateSeg[peaks], lty = 2, col = "grey80")
              
              pars    <- list(rel_amp_frac = 0.15, rel_peak_frac = 0.1, min_seg_amplitude = 0.15)
              segs    <- tryCatch(do.call("rbind", GetSegs(peaks, xSmooth, pars)), error = function(e) NULL)
              pSeg    <- median(s$pow, na.rm = T)
              
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
                               auc(t, y, type = 'spline'))
                
                
                out <- c(year = year, max = max, amp = round(amp, 2), area = round(area,2),
                         as.numeric(format(as.POSIXct(c(q10gup, q50gup, q10sen, q50sen), origin = "1970-01-01"), "%j")))
                out[5:6] <- ifelse(out[5:6]>=max, out[5:6] - 365, out[5:6])
                out[7:8] <- ifelse(out[7:8]<=max, 365 + out[7:8], out[7:8])
                
                out
                
              } else {
                c(year = year, per = per, max = max(xSmooth), pow = pSeg, amp = diff(range(xSmooth)), area = NA, rep(NA, 4))
              }
              
            } else {
              c(year = year, per = per, 
                             max = as.numeric(quantile(s[,-c(1:5)][,weigth==1], prob = 0.975 , na.rm = T)), 
                             pow = median(s$pow, na.rm = T),
                             amp = as.numeric(diff(quantile(s[,-c(1:5)][,weigth==1], prob = c(0.025, 0.975)) , na.rm = T)), area = NA, rep(NA, 4))
            }
            
      }))
      
      }  else {
        phen0 = cbind(year = 1981:2020, per  = per,
                      max  = as.numeric(diff(quantile(evi[weigth==1,], prob = c(0.025, 0.975) , na.rm = T))), 
                      pow  = pow,
                      amp  = as.numeric(quantile(evi[weigth==1,], prob = 0.975, na.rm = T)), 
                      area = NA, matrix(NA, ncol = 4, nrow = length(1981:2020)))
    } 
      
    
  } else {
    phen0 = cbind(year = 1981:2020, per  = NA,
                                    max  = as.numeric(diff(quantile(evi[weigth==1,], prob = c(0.025, 0.975) , na.rm = T))), 
                                    pow  = NA,
                                    amp  = as.numeric(quantile(evi[weigth==1,], prob = 0.975, na.rm = T)), 
                                    area = NA, matrix(NA, ncol = 4, nrow = length(1981:2020)))
  }
  
  
  if(period==1) {
    phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
  } else {
    pk   <- round(median(unlist(aggregate(phen0[,3], by = list(phen0[,1]), FUN = function(a) which.max(a))$x), na.rm = T), 0)
    phen <- do.call("rbind", lapply(split(as.data.frame(phen0), round(phen0[,1],0)), function(a) a[1,]))
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
  
}, mc.cores = 4)


phenA <- tryCatch(abind::abind(lapply(1:ncol(pxlPhen[[1]]), function(l) matrix(unlist(sapply(pxlPhen, function(x) x[,l])), ncol = nrow(pxlPhen[[1]]), byrow = T)), along = 3),
                  error = function(e) NULL)

phenR <- rasterFromXYZ(cbind(outBatch$crds[outBatch$crds[,3],1:2],
                             (apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
                                !outBatch$crds[,ncol(outBatch$crds)])[outBatch$crds[,3]]), crs = CRS("+proj=longlat"))


rTest <- phenR
if(!is.null(phenA)) rTest[rTest[]==1] <- apply(phenA[,,2], 1, median, na.rm = T)

png(paste0("/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/phenPNG_", batch,".png"), width = 400, height = 400, res = 100)
plot(rTest)
dev.off()

drive_upload(paste0("/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/phenPNG_", batch,".png"),
             path = "SeasonalChange/phenBatches/", name = paste0("phenPNG_", batch, ".png"), verbose = FALSE, overwrite = T)


names(phenR) <- paste0("batch_", batch)
phenOut <- list(dat = phenA, raster = phenR)
save(phenOut, file = paste0(phen_dir, "phenBatch_", batch, ".rda"))
