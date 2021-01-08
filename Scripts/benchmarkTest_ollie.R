#### benchmark test
.libPaths(c(.libPaths(), "/home/ollie/slisovsk/RServerLibs/3.6"))

options(warn=-1)
library(sf)
library(dplyr)
library(raster)
library(rgdal)
library(rgeos)
library(zoo)
library(biwavelet)
library(MESS)
library(zoo)
library(geodist)

source('/work/ollie/slisovsk/Projects/SeasonalChange/Scripts/functions.R', echo=FALSE)


batch_dir <- "/work/ollie/slisovsk/Projects/SeasonalChange/Results/Batches/"
load("/work/ollie/slisovsk/Projects/SeasonalChange/Results/dateSequence.rda") 

batch <- 1027
load(paste0(batch_dir, "Batch_", batch, ".rda"))

inBatch_sf = st_as_sf(outBatch$crds[,c(1:2)], coords = c("x", "y"), agr = "constant") %>% st_set_crs("+proj=longlat")
inBatch_sf$onLand <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
  !outBatch$crds[,ncol(outBatch$crds)] &
  outBatch$crds[,3]
inBatch_sf$inBatch <- apply(outBatch$crds[,-c(1:3, ncol(outBatch$crds))], 1, any) &
  !outBatch$crds[,ncol(outBatch$crds)]

evalPxl <- function(pxl) {
  
  # pxl   <- which(inBatch_sf$onLand)[28]
  dst   <- subset(data.frame(ind  = 1:sum(inBatch_sf$inBatch), 
                             dist = suppressMessages(geodist(st_coordinates(inBatch_sf$geometry[inBatch_sf$inBatch]), 
                                                             st_coordinates(inBatch_sf$geometry[pxl,]), measure = "cheap"))/1000), dist<15)
  
  # plot(inBatch_sf$geometry, pch = 16, cex = 0.6)
  # plot(inBatch_sf$geometry[pxl,], pch = 16, col = "cyan", add = T)
  # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch)[dst[,1]],], pch = 16, col = "magenta", add =T, cex = 0.55)
  
  weigth <- approx(c(0, max(dnorm(dst[,2], 0, 4))), c(0,1), dnorm(dst[,2], 0, 4), rule = 2)$y
  
  dat    <- outBatch$dat[which(inBatch_sf$inBatch)[which(dst<15)],,]
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
        
        if(length(dateSeg)>20 & sum(is.na(s[,-c(1:5)][,weigth==1]))<nrow(s)*0.33) {
          
          pwL  <- apply(log2(abs(wt$power/wt$sigma2))[,s$id], 1, median, na.rm = T)
          sig  <- apply(wt$signif[,s$id], 1, median, na.rm = T)
          
          wt.pks <- FindPeaks(pwL)
          wt.s   <- cbind(pwL[wt.pks], wt$period[wt.pks], (sig>=1)[wt.pks])
          wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[order(c(wt.s[wt.s[,3]==1,1], rep(0,3)), decreasing = T),][1:2,]
          
          out <- c(year = year,                                
                   per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
                   sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
                   per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
                   sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
                   max = round(quantile(s[,-c(1:5)][,weigth==1], probs = 0.98, na.rm = T),2), 
                   amp = round(diff(quantile(s[,-c(1:5)][,weigth==1], probs = c(0.02, 0.98), na.rm = T)),2))
          
          out
          
        
          } else {
          
          c(year = year,
            per1 = NA,
            sig1 = NA,
            per2 = NA,
            sig2 = NA,
            max  = NA,
            amp  = NA)
          
        }

      }))
      
    
      }  else {

      phen0 =  cbind(year = 1981:2020,                      
                     per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),
                     sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),
                     per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),
                     sig2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),
                     max  = NA,
                     amp = NA)
    }
    
    
  
    } else {

    phen0 =  cbind(year = 1981:2020,                      
                   per1 = NA,
                   sig1 = NA,
                   per2 = NA,
                   sig2 = NA,
                   max  = NA,
                   amp  = NA)
  }
  
  
  if(exists("wt.sig") & any(wt.sig[,3]==1) & wt.sig[1,2]%in%c(42:62)) {
    phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
  } else {
    phen <- do.call("rbind", lapply(split(as.data.frame(phen0), round(phen0[,1],0)), function(a) a[1,]))
  }
  
  merge(data.frame(year = 1981:2020), as.data.frame(phen), all.x = T)[,-1]
  
}

### mclapply
library(parallel)
pxlPhen <- mclapply(which(inBatch_sf$onLand), evalPxl, mc.cores = parallel::detectCores())

