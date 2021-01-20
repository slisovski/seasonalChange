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


### map 
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_geometry()
lbox <- st_bbox(c(xmin = -30, xmax = -90, ymin = -57, ymax = 13), crs = 4326) %>% st_as_sfc()
r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
pol <- st_as_sfc(rasterToPolygons(r0)) %>% st_set_crs(4326)
inBBB <- suppressMessages(unlist(st_intersects(lbox, pol)))


## batches
batch_dir <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/Batches/"
phen_dir  <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/"
batches   <- list.files(batch_dir, pattern = ".rda")
batchID   <- sapply(strsplit(batches, "_"), function(x) as.numeric(substring(x[[2]], 1, nchar(x[[2]]) -4)))
load("/work/ollie/slisovsk/Projects/SeasonalChange/Results/dateSequence.rda") 

batch <- batchID[batchID%in%inBBB][as.numeric(args[1])]

# if(!file.exists(paste0(phen_dir, "phenBatch_", batch, ".rda"))) {

load(paste0(batch_dir, "Batch_", batch, ".rda"))

inBatch_sf = st_as_sf(outBatch$crds[,1:2], coords = c("x", "y"), crs = 4326, agr = "constant")
  inBatch_sf$onLand <- outBatch$crds$land &
    !outBatch$crds[,ncol(outBatch$crds)] &
    outBatch$crds[,3]
  inBatch_sf$inBatch <- outBatch$crds$land &
    !outBatch$crds[,ncol(outBatch$crds)]
  
  # plot(inBatch_sf$geometry, pch = 16, cex = 0.5)
  # plot(inBatch_sf$geometry[inBatch_sf$inBatch,], pch = 16, cex = 0.5, col = "orange", add = T)
  # plot(inBatch_sf$geometry[inBatch_sf$onLand,], pch = 16, cex = 0.5, col = "green", add = T)

pxlPhen <- mclapply(which(inBatch_sf$onLand), function(pxl) {
  
  # pxl <- which(inBatch_sf$onLand)[28]
  
  dst   <- suppressMessages(geodist(st_coordinates(inBatch_sf$geometry[inBatch_sf$inBatch]), st_coordinates(inBatch_sf$geometry[pxl,]), measure = "cheap"))/1000
  inbfr <- which(dst<15)

  # plot(inBatch_sf$geometry[pxl,], pch = 16, col = "cyan", add =T)
  # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch)[inbfr],], pch = 16, col = "magenta", add =T)

  weigth <- approx(c(0, max(dnorm(dst, 0, 4))), c(0,1), dnorm(dst[inbfr], 0, 4), rule = 2)$y

  dat    <- outBatch$dat[which(inBatch_sf$inBatch)[inbfr],,]
  evi    <- ifelse(dat[,,1] <= 0 , NA, dat[,,1])

  med <- apply(evi, 2, median, na.rm = T)
  
  if(sum(is.na(med))<length(med)*0.2) {
    
    wt   <- wt(cbind(1:length(med), na.approx(med)))
    pwL  <- apply(log2(abs(wt$power/wt$sigma2)), 1, median, na.rm = T)
    sig  <- apply(wt$signif, 1, median, na.rm = T)
    
    wt.pks <- FindPeaks(pwL)
    wt.s   <- cbind(pwL[wt.pks], wt$period[wt.pks], (sig>=1)[wt.pks])[order(pwL[wt.pks], decreasing = T),]
    wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[1:2,]  

    
    if(round(wt.sig[1,2],0)%in%c(42:62)) {
      period = 52
    } else {
      if(round(round(wt.sig[1,2],0),0)%in%c(16:36)) {
        period = 26
      } else period = NA
    }
    
    datCurve <- merge(data.frame(year = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%Y")),
                                   week = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%U"))),
                      data.frame(year = as.numeric(format(dates, "%Y")), week = as.numeric(format(dates, "%U")),
                                   date = dates, id = 1:length(dates),
                                   evi  = med, t(evi)), all.x = T)
    
    if(!is.na(period)) {
      
    fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = period, Mx = na.approx(datCurve$evi) - median(med, na.rm = T), sd = 0.01)
    curve <- fit0$par[1]*cos(pi*((1:nrow(datCurve))/(nrow(datCurve)/((nrow(datCurve)/period)*2))) +
                                 (pi+fit0$par[2])) +  mean(med, na.rm=T)
      
    mins <- unique(c(1, which(diff(sign(diff(-curve)))==-2)+1, length(curve)))
    maxs <- which(diff(sign(diff(curve)))==-2)+1
 
    } else {
      maxs <- mins <- NA
    }
    
    if((diff(quantile(med, prob = c(0.025,0.975), na.rm = T))>0.15) & !is.na(period) & length(mins)>20 & length(maxs)>20) {

          segL <- apply(do.call("rbind", lapply(maxs, function(x) mins[c(1:length(mins))[order(abs(x-mins))][1:2]]+c(-10,10))), 1,
                        function(x) datCurve[ifelse(x[1]<1, 1, x[1]):ifelse(x[2]>length(curve), length(curve), x[2]),])

          q50 <- median(sapply(segL, function(x) quantile(x[,5], prob = 0.5, na.rm = T)), na.rm = T)

          phen0 <- do.call("rbind", mclapply(segL, function(s) {

            # s <- segL[[77]]
            # matplot(s$date, s[,-c(1:5)], type= "l", lty = 1, lwd = 1, col = adjustcolor("darkgreen", alpha.f = 0.5), pch = 16, xaxt = "n")
            # axis(1, at = as.numeric(seq(min(s$date, na.rm = T), max(s$date, na.rm = T), length = 10)), labels = format(seq(min(s$date, na.rm = T), max(s$date, na.rm = T), length = 10), "%Y-%m-%d"))

            dateSeg <- na.approx(s$date)
            year    <- median(as.numeric(format(s$date, "%Y")), na.rm = T)

            if(length(dateSeg)>20 & sum(is.na(s[,-c(1:5)][,weigth==1]))<nrow(s)*0.33) {
                
                pwL  <- apply(log2(abs(wt$power/wt$sigma2))[,s$id], 1, median, na.rm = T)
                sig  <- apply(wt$signif[,s$id], 1, median, na.rm = T)
    
                wt.pks <- FindPeaks(pwL)
                wt.s   <- cbind(pwL[wt.pks], wt$period[wt.pks], (sig>=1)[wt.pks])
                wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[order(c(wt.s[wt.s[,3]==1,1], rep(0,3)), decreasing = T),][1:2,]
    
                # plot(wt$period, pwL, type = "o")
                # points(wt$period, pwL, col = ifelse(sig>=1, "red", "grey90"), pch = 16)


    #           tmpDat <- subset(data.frame(x = rep(s$date, ncol(s)-5),
    #                                       y = unlist(t(c(s[,-c(1:5)]))),
    #                                       w = rep(weigth, each = nrow(s))), !is.na(y))
    #           
    #           dtsSm   <- seq(min(dateSeg), max(dateSeg), by = 24*60*60)
    #           
    #           spl     <- smooth.spline(x = tmpDat$x, y = tmpDat$y, spar = 0.3, w = tmpDat$w)
    #           xSmooth <- predict(spl, dtsSm)$y
    #           peaks   <- FindPeaks(xSmooth)
    #           
    #           # lines(dtsSm, xSmooth, type= "l", lwd = 6, col = "orange")
    #           # abline(v = dtsSm[peaks], lty = 2, col = "grey80")
    #           
    #           pars    <- list(rel_amp_frac = 0, rel_peak_frac = 0, min_seg_amplitude = 0.15)
    #           segs    <- tryCatch(do.call("rbind", GetSegs(peaks, xSmooth, pars)), error = function(e) NULL)
    # 
    #           # apply(segs, 1, function(x) rect(dtsSm[x[1]], 0, dtsSm[x[3]], 1, col = adjustcolor("grey80", alpha.f = 0.2), border =  "grey10"))
    #           
    #           if(!is.null(segs)) {
    #             
    #             if(nrow(segs)>1) {
    #               
    #               max_in   <- min(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
    #               max_out  <- max(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
    #               
    #               seqRan   <- c(min(segs[order(xSmooth[segs[,2]], decreasing = T),1]),
    #                             max(segs[order(xSmooth[segs[,2]], decreasing = T),3]))
    #               
    #             } else {
    #               
    #               max_in <- max_out  <- segs[2]
    #               seqRan <- segs[c(1,3)]
    #               
    #             }
    #             
    #             max    <- as.numeric(format(as.POSIXct(dtsSm, origin = "1970-01-01")[max_in], "%j"))
    #             # abline(v = as.POSIXct(dtsSm, origin = "1970-01-01")[max_in], lty = 3, col = "red")
    #             amp    <- diff(range(xSmooth[seqRan[1]:seqRan[2]]))
    #             
    #             q10gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
    #             # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,], min(y) + diff(range(y))*0.1), v = q10gup, col = "blue")
    #             q50gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(q50, rev(y), gup = F, first_greater = T)])
    #             # abline(h = q50, v = q50gup, col = "cyan")
    #             q90gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
    #             # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,], max(y) - diff(range(y))*0.1), v = q90gup, col = "firebrick")
    #             
    #             q90sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
    #             # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],], max(y) - diff(range(y))*0.1), v = q90sen, col = "blue")
    #             q50sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(q50, rev(y), gup = T, first_greater = F)])
    #             # abline(h = q50, v = q50sen, col = "cyan")
    #             q10sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
    #             # abline(h = with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],], min(y) + diff(range(y))*0.1), v = q10sen, col = "firebrick")
    #             
    #             area   <- with(data.frame(t = 1:length(dtsSm), y = xSmooth)[seqRan[1]:seqRan[2],],
    #                            MESS::auc(t, y, type = 'spline'))
    #             
    #             
    #             out <- c(year = year,                                  #1
    #                      per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),  #2
    #                      sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),  #3
    #                      per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),  #4
    #                      sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),  #5
    #                      per3 = ifelse(wt.sig[3,3], wt.sig[3,2], NA),  #6
    #                      sig3 = ifelse(wt.sig[3,3], wt.sig[3,1], NA),  #7
    #                      max = max,                                    #8
    #                      amp = round(amp, 2),                          #9
    #                      area = round(area,2),                         #10
    #                      as.numeric(format(as.POSIXct(c(q10gup, q50gup, q90gup, q90sen, q50sen, q10sen), origin = "1970-01-01"), "%j"))) #11, 12, 13, 14
    #             out[11:13] <- ifelse(out[11:13]>=max, out[11:13] - 365, out[11:13])
    #             out[14:16] <- ifelse(out[14:16]<=max, 365 + out[14:16], out[14:16])

              out <- c(year = year,                                
                       per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
                       sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
                       per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
                       sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
                       max = round(quantile(s[,-c(1:5)][,weigth==1], probs = 0.98, na.rm = T),2), 
                       amp = round(diff(quantile(s[,-c(1:5)][,weigth==1], probs = c(0.02, 0.98), na.rm = T)),2))
              
              
              out

            } else {
    #             c(year = year, 
    #               per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
    #               sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
    #               per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
    #               sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
    #               per3 = ifelse(wt.sig[3,3], wt.sig[3,2], NA),
    #               sig3 = ifelse(wt.sig[3,3], wt.sig[3,1], NA),
    #               max = NA, 
    #               amp = NA, area = NA, rep(NA, 6))
                   c(year = year,
                     per1 = NA,
                     sig1 = NA,
                     per2 = NA,
                     sig2 = NA,
                     max  = NA,
                     amp  = NA)
              }
    #           
    #         } else {
    #           c(year = year,                      
    #             per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
    #             sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
    #             per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
    #             sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
    #             per3 = ifelse(wt.sig[3,3], wt.sig[3,2], NA),
    #             sig3 = ifelse(wt.sig[3,3], wt.sig[3,1], NA),
    #             max = NA,
    #             amp = NA, area = NA, rep(NA, 6))
    #         }
    #         
      }))
  
      }  else {
    #     phen0 = cbind(year = 1981:2020, 
    #                   per1 = NA,
    #                   sig1 = NA,
    #                   per2 = NA,
    #                   sig2 = NA,
    #                   per3 = NA,
    #                   sig3 = NA,
    #                   max  = NA, 
    #                   amp  = NA, 
    #                   area = NA, matrix(NA, ncol = 6, nrow = length(1981:2020)))
        phen0 =  cbind(year = 1981:2020,                      
                       per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
                       sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
                       per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
                       sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
                       max  = NA,
                       amp = NA)
    }
    
    
  } else {
    # phen0 = cbind(year = 1981:2020, 
    #               per1 = NA,
    #               sig1 = NA,
    #               per2 = NA,
    #               sig2 = NA,
    #               per3 = NA,
    #               sig3 = NA,
    #               max  = NA, 
    #               amp  = NA, 
    #               area = NA, matrix(NA, ncol = 6, nrow = length(1981:2020)))
    
    phen0 =  cbind(year = 1981:2020,                      
                   per1 = NA,
                   sig1 = NA,
                   per2 = NA,
                   sig2 = NA,
                   max  = NA,
                   amp  = NA)
  }
  
  
  if(!is.na(period) & period==52) {
    phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
  } else {
    phen <- do.call("rbind", lapply(split(as.data.frame(phen0), round(phen0[,1],0)), function(a) a[1,]))
  }
  
  merge(data.frame(year = 1981:2020), as.data.frame(phen), all.x = T)[,-1]
  
}, mc.cores = 15)

# tt <- sapply(pxlPhen, function(x) is.data.frame(x))
# pxlPhen[[which(!tt)[1]]]

phenA <- tryCatch(abind::abind(lapply(1:ncol(pxlPhen[[1]]), function(l) matrix(unlist(sapply(pxlPhen, function(x) x[,l])), ncol = nrow(pxlPhen[[1]]), byrow = T)), along = 3),
                  error = function(e) NULL)

phenR <- rasterFromXYZ(cbind(outBatch$crds[outBatch$crds[,3],1:2],
                             (apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
                                !outBatch$crds[,ncol(outBatch$crds)])[outBatch$crds[,3]]), crs = CRS("+proj=longlat"))


# rTest <- phenR
# if(!is.null(phenA)) rTest[rTest[]==1] <- apply(phenA[,,13], 1, median, na.rm = T)

# png(paste0("/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/phenPNG_", batch,".png"), width = 400, height = 400, res = 100)
# plot(rTest)
# dev.off()
# 
# drive_upload(paste0("/work/ollie/slisovsk/Projects/SeasonalChange/Results/phenBatches/phenPNG_", batch,".png"),
#              path = "SeasonalChange/phenBatches/", name = paste0("phenPNG_", batch, ".png"), verbose = FALSE, overwrite = T)


names(phenR) <- paste0("batch_", batch)
phenOut <- list(dat = phenA, raster = phenR)
save(phenOut, file = paste0(phen_dir, "phenBatch_", batch, ".rda"))

drive_upload(paste0(phen_dir, "phenBatch_", batch, ".rda"),
             path = "SeasonalChange/phenBatchesRDA/", name = paste0("phenBatch_", batch, ".rda"), verbose = FALSE, overwrite = T)

} ## file exists