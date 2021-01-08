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
  
evalPxl <- function(pxl) {
    
    # pxl   <- which(inBatch_sf$onLand)[1]
  
    dst   <- subset(data.frame(ind  = 1:sum(inBatch_sf$inBatch), 
                               dist = suppressMessages(geodist(st_coordinates(inBatch_sf$geometry[inBatch_sf$inBatch]), 
                                                               st_coordinates(inBatch_sf$geometry[pxl,]), measure = "cheap"))/1000), dist<15)
    
    # plot(inBatch_sf$geometry, pch = 16, cex = 0.6)
    # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch),], pch = 16, col = "orange", add = T)
    # plot(inBatch_sf$geometry[which(inBatch_sf$onLand),], pch = 16, cex = 0.5, col = "darkgreen", add = T)
    # plot(inBatch_sf$geometry[pxl,], pch = 16, col = "red", add =T, cex = 0.85)
    # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch)[dst[,1]],], pch = "x", add =T, cex = 0.85)
    
    weigth <- approx(c(0, max(dnorm(dst[,2], 0, 4))), c(0,1), dnorm(dst[,2], 0, 4), rule = 2)$y
    
    dat    <- outBatch$dat[which(inBatch_sf$inBatch)[dst[,1]],,]
    evi    <- ifelse(dat[,,1] <= 0 , NA, dat[,,1])
    
    med <- apply(evi, 2, median, na.rm = T)
    
    if(sum(is.na(med))<length(med)*0.15) {
      
      wt   <- wt(cbind(1:length(med), na.approx(med)))
      pwL  <- apply(log2(abs(wt$power/wt$sigma2)), 1, median, na.rm = T)
      sig  <- apply(wt$signif, 1, median, na.rm = T)
      
      wt.pks <- FindPeaks(pwL)
      wt.s   <- cbind(pwL[wt.pks], round(wt$period[wt.pks],0), (sig>=1)[wt.pks])[order(pwL[wt.pks], decreasing = T),]
      wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[1:2,]  
      
      datCurve <- merge(data.frame(year = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%Y")),
                                   week = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%U"))),
                        data.frame(year = as.numeric(format(dates, "%Y")), week = as.numeric(format(dates, "%U")),
                                   date = dates, id = 1:length(dates),
                                   evi  = na.approx(med), t(evi)), all.x = T)
      
      if(any(wt.sig[,3]==1)) {
        
        fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = wt.sig[1,2], Mx = datCurve$evi - median(med, na.rm = T), sd = 0.01)
        curve <- fit0$par[1]*cos(pi*((1:nrow(datCurve))/(nrow(datCurve)/((nrow(datCurve)/wt.sig[1,2])*2))) +
                                   (pi+fit0$par[2])) +  mean(med, na.rm=T)
        
        mins <- unique(c(1, which(diff(sign(diff(-curve)))==-2)+1, length(curve)))
        maxs <- which(diff(sign(diff(curve)))==-2)+1
        
      } else {
        mins <- maxs <- NA
      }
      
      if((diff(quantile(med, prob = c(0.025,0.975), na.rm = T))>0.15) & length(mins)>20 & length(maxs)>20) {
        
        segL <- apply(do.call("rbind", lapply(maxs, function(x) mins[c(1:length(mins))[order(abs(x-mins))][1:2]]+c(-10,10))), 1,
                      function(x) datCurve[ifelse(x[1]<1, 1, x[1]):ifelse(x[2]>length(curve), length(curve), x[2]),])
        
        q50 <- median(sapply(segL, function(x) quantile(x[,5], prob = 0.5, na.rm = T)), na.rm = T)
        
        phen0 <- do.call("rbind", lapply(segL, function(s) {
          
          dateSeg <- na.approx(s$date)
          year    <- median(as.numeric(format(s$date, "%Y")), na.rm = T)
          
          if(length(dateSeg)>20 & sum(is.na(s[,-c(1:5)]))<length(unlist(c(s[,-c(1:5)])))*0.35) {
            
            pwL  <- apply(log2(abs(wt$power/wt$sigma2))[,s$id], 1, median, na.rm = T)
            sig  <- apply(wt$signif[,s$id], 1, median, na.rm = T)
            
            wt.pks <- FindPeaks(pwL)
            wt.s   <- cbind(pwL[wt.pks], wt$period[wt.pks], (sig>=1)[wt.pks])
            wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[order(c(wt.s[wt.s[,3]==1,1], rep(0,3)), decreasing = T),][1:2,]
            
            # plot(wt$period, pwL, type = "o")
            # points(wt$period, pwL, col = ifelse(sig>=1, "red", "grey90"), pch = 16)
            
            
            tmpDat <- subset(data.frame(x = rep(s$date, ncol(s)-5),
                                        y = unlist(t(c(s[,-c(1:5)]))),
                                        w = rep(weigth, each = nrow(s))), !is.na(y))
         
            dtsSm   <- seq(min(dateSeg), max(dateSeg), by = 24*60*60)
                       
            spl     <- smooth.spline(x = tmpDat$x, y = tmpDat$y, spar = 0.3, w = tmpDat$w)
            xSmooth <- predict(spl, dtsSm)$y
            peaks   <- FindPeaks(xSmooth)
               
              # plot(tmpDat$x, tmpDat$y, pch = 16, cex = 0.5)
              # lines(dtsSm, xSmooth, type= "l", lwd = 6, col = "orange")
              # abline(v = dtsSm[peaks], lty = 2, col = "grey80")

            pars    <- list(rel_amp_frac = 0.15, rel_peak_frac = NULL, min_seg_amplitude = 0.1)
            segs    <- tryCatch(do.call("rbind", GetSegs(peaks, xSmooth, pars)), error = function(e) NULL)
             
              # apply(segs, 1, function(x) rect(dtsSm[x[1]], 0, dtsSm[x[3]], 1, col = adjustcolor("grey80", alpha.f = 0.2), border =  "grey10"))
       
            if(!is.null(segs)) {

              if(nrow(segs)>1) {

                  max_in   <- min(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
                  max_out  <- max(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])

                  seqRan   <- c(min(segs[order(xSmooth[segs[,2]], decreasing = T),1]),
                                max(segs[order(xSmooth[segs[,2]], decreasing = T),3]))
                          
              } else {
                max_in <- max_out  <- segs[2]
                seqRan <- segs[c(1,3)]
              }
              
              
              max    <- as.numeric(format(as.POSIXct(dtsSm, origin = "1970-01-01")[max_in], "%j"))
              # abline(v = as.POSIXct(dtsSm, origin = "1970-01-01")[max_in], lty = 3, col = "red")
              amp    <- diff(range(xSmooth[seqRan[1]:seqRan[2]]))
              
              q10gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
              # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,], min(y) + diff(range(y))*0.1), v = q10gup, col = "blue")
              q50gup <- suppressWarnings(with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(q50, rev(y), gup = F, first_greater = T)]))
              # abline(h = q50, v = q50gup, col = "cyan")
              q90gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
              # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,], max(y) - diff(range(y))*0.1), v = q90gup, col = "firebrick")
              
              q90sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
              # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],], max(y) - diff(range(y))*0.1), v = q90sen, col = "blue")
              q50sen <- suppressWarnings(with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(q50, rev(y), gup = T, first_greater = F)]))
              # abline(h = q50, v = q50sen, col = "cyan")
              q10sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
              # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],], min(y) + diff(range(y))*0.1), v = q10sen, col = "firebrick")
              
              area   <- with(data.frame(t = 1:length(dtsSm), y = xSmooth)[seqRan[1]:seqRan[2],],
                             MESS::auc(t, y, type = 'spline'))
              
              
              out <- c(year = year,                                  #1
                       per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),  #2
                       sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),  #3
                       per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),  #4
                       per2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),  #5
                       max = max,                                    #6
                       amp = round(amp, 2),                          #7
                       area = round(area,2),                         #8
                       as.numeric(format(as.POSIXct(c(q10gup, q50gup, q90gup, q90sen, q50sen, q10sen), origin = "1970-01-01"), "%j"))) #9, 10, 11, 12, 13, 14
              
              # out[9:11] <- ifelse(out[9:11]>=max, out[9:11] - 365, out[9:11])
              # out[12:14] <- ifelse(out[12:14]<=max, 365 + out[12:14], out[12:14])
              
              out
              
            } else  {
              c(year = year,
                per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
                sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
                per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
                sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
                max = NA,
                amp = NA, area = NA, rep(NA, 6))
            }
            
            
            
          } else {
            
            c(year = year,
              per1 = NA,
              sig1 = NA,
              per2 = NA,
              sig2 = NA,
              max  = NA,
              amp  = NA, area = NA, rep(NA, 6))
            
          }
          
        }))
        
      }  else {
        
        phen0 =  cbind(year = 1981:2020,                      
                       per1 = NA,
                       sig1 = NA,
                       per2 = NA,
                       sig2 = NA,
                       max  = NA,
                       amp  = NA,
                       area = NA, rep(NA, 6))
      }
      
    } else {
      
      phen0 =  cbind(year = 1981:2020,                      
                     per1 = NA,
                     sig1 = NA,
                     per2 = NA,
                     sig2 = NA,
                     max  = NA,
                     amp  = NA,
                     area = NA, rep(NA, 6))
    }
    
    
    if(exists("wt.sig") & any(wt.sig[,3]==1) & wt.sig[1,2]%in%c(42:62)) {
      phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
    } else {
      phen <- do.call("rbind", lapply(split(as.data.frame(phen0), round(phen0[,1],0)), function(a) a[1,]))
    }
    
    merge(data.frame(year = 1981:2020), as.data.frame(phen), all.x = T)[,-1]
    
  }
})


### map 
land <- read_sf("Data/GeoDat/ne_50m_land/ne_50m_land.shp") %>% st_geometry()
lbox <- st_bbox(c(xmin = -30, xmax = -90, ymin = -57, ymax = 13), crs = 4326) %>% st_as_sfc()
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
batch <- batchID[batchID%in%inBBB][as.numeric(args[1])]

# if(!file.exists(paste0(phen_dir, "phenBatch_", batch, ".rda"))) {

  load(paste0(batch_dir, "Batch_", batch, ".rda"))
  
  inBatch_sf = st_as_sf(outBatch$crds[,c(1:2)], coords = c("x", "y"), agr = "constant") %>% st_set_crs("+proj=longlat")
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

# }