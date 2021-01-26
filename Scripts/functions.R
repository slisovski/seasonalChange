FindPeaks <- function(x, mag_order=T){
  # Function to identify peaks in time series x (or troughs if x=-x), supports "flat top" peaks
  # if mag_order is TRUE, peaks are returned in order of increasing magnitude (of x)
  d <- diff(x)
  d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec
  peaks <- unlist(gregexpr("12", paste(d_code, collapse=""))) # no match is -1
  if(peaks[1] == -1) peaks <- NULL
  flat_peaks <- unlist(gregexpr("10+2", paste(d_code, collapse=""))) # no match is -1
  if(flat_peaks[1] == -1) flat_peaks <- NULL
  d_code_rle <- rle(d_code)
  flat_peaks <- flat_peaks + round(d_code_rle$l[match(flat_peaks, cumsum(d_code_rle$l)) + 1] / 2)
  # all_peaks <- c(ifelse(peaks[1] == -1, NULL, peaks + 1), ifelse(flat_peaks[1] == -1, NULL, flat_peaks + 1))
  peaks <- sort(c(peaks + 1, flat_peaks + 1))
  if(mag_order) return(peaks[order(x[peaks])])
  return(peaks)
}

GetSegs <- function(peaks, x, pars, peak=NA){
  # identifies valid increasing-decreasing segments in x subject to the parameters in pars
  # returns a list of segments: c(start, peak, end). DON'T call directly w/ peak!=NA
  # NOTE: returned segments will not necessarily be in order, and may not completely partition x
  
  # ensure that peaks are in increasing order of x's magnitude
  tmp_peaks <- peaks[order(x[peaks])] # so we only have to sort once if they're in the wrong order
  if(!identical(tmp_peaks, peaks)) peaks <- tmp_peaks
  
  # if no peak is specified, we start at the beginning
  if(is.na(peak)) peak <- peaks[1]
  
  # get the next largest peak; will be NA if this peak is the highest (last one to do)
  next_highest_peak <- peaks[which(peaks == peak) + 1]
  
  # check if we're doing C5-style relative amplitude and peak identification
  # if(!is.na(pars$rel_amp_frac) & !is.na(pars$rel_peak_frac)){
  #   global_max <- max(x, na.rm=T)
  #   seg_thresh <- (global_max - min(x, na.rm=T)) * pars$rel_amp_frac
  #   peak_thresh <- global_max * pars$rel_peak_frac
  # }else{
  #   seg_thresh <- pars$min_seg_amplitude
  #   peak_thresh <- 0
  # }
  
  # we could have any combination of rel_amp_frac, rel_peak_frac, and min_seg_amplitude specified
  # initialize seg_thresh and peak_thresh to zero
  # determine the "global max/min", if peak_frac is specified, set it, if amp_frac is specified, set it
  # if min_seg_amplitude is set, choose the max of that and amp_frac
  seg_thresh <- peak_thresh <- 0
  global_max <- max(x, na.rm=T)
  global_min <- min(x, na.rm=T)
  if(!is.na(pars$rel_amp_frac)) seg_thresh <- (global_max - global_min) * pars$rel_amp_frac
  #if(!is.na(pars$rel_peak_frac)) peak_thresh <- global_max * pars$rel_peak_frac
  if(!is.na(pars$min_seg_amplitude)) seg_thresh <- max(pars$min_seg_amplitude, seg_thresh)
  
  # checks if the period preceding the peak covers enough amplitude
  # search before the peak up to the maximum of: previous peak, the head of x, or the peak - max_increase_length
  previous_peaks <- peaks[peaks - peak < 0]
  previous_peak <- NA
  if(length(previous_peaks) > 0) previous_peak <- max(previous_peaks)
  search_start <- max(1, peak - pars$max_increase_length, previous_peak, na.rm=T)
  search_end <- peak
  # get the index of the closest minimum value within the search window
  # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
  # in the event of repeated minimum values, we take the closest one here
  inc_min_ind <- max(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
  seg_amp <- x[peak] - x[inc_min_ind] # get the increasing segment amplitude
  # if(seg_amp > pars$min_seg_amplitude){
  if((seg_amp >= seg_thresh) & (x[peak] >= peak_thresh)){
    # check for a valid decreasing segment
    next_peaks <- peaks[peaks - peak > 0]
    next_peak <- NA
    if(length(next_peaks) > 0) next_peak <- min(next_peaks)
    # search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
    search_start <- peak
    search_end <- min(length(x), peak + pars$max_decrease_length, next_peak, na.rm=T)
    # get the index of the closest minimum value within the search window
    # NOTE: see above note about finding troughs instead
    dec_min_ind <- min(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
    seg_amp <- x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
    # if(seg_amp > pars$min_seg_amplitude){
    if(seg_amp >= seg_thresh){
      # we found a valid segment, store it as a list with a single vector: c(start, peak, end)
      tmp_seg <- list(c(inc_min_ind, peak, dec_min_ind))
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(c(tmp_seg, GetSegs(peaks, x, pars, peak=next_highest_peak)))
      }else{
        # that was the last peak, and it was valid
        return(tmp_seg) # covers the case where there's only one valid peak
      }
    }else{
      # increase was valid, but decrease was not
      peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(GetSegs(peaks, x, pars, peak=next_highest_peak))
      }else{
        # that was the last peak, and it was invalid
        return(NULL)
      }
    }
  }else{
    # increase segment not valid
    peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
    # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
    if(!is.na(next_highest_peak)){
      return(GetSegs(peaks, x, pars, peak=next_highest_peak))
    }else{
      # that was the last peak, and it was invalid
      return(NULL)
    }
  }
}

GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
  # returns the index of the first/last value  of x that is greater/less than the value of thresh.
  # If gup is False (greendown) then it returns the first/last value of x that is less/greater than
  # the value of thresh. first/last and greater/less determined by first_greater
  # NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round the threshold and each
  # of the evi values to 6 decimal places to compensate
  
  if(gup){
    if(first_greater){
      return(min(which(round(x, 6) >= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) <= round(thresh_value, 6))))
    }
  }else{
    if(first_greater){
      return(min(which(round(x, 6) <= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) >= round(thresh_value, 6))))
    }
  }
}

lsCos <- function(params, f, Mx, sd = 0.001) {
  fit  <- params[1]*cos(pi*((1: length(Mx))/(length(Mx)/((length(Mx)/f)*2))) + (pi+params[2]))
  -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
}


evalPxl <- function(pxl) {
  
  # pxl   <- which(inBatch_sf$onLand)[200]

  dst   <- subset(data.frame(ind  = 1:sum(inBatch_sf$inBatch),
                             dist = suppressMessages(geodist(st_coordinates(inBatch_sf$geometry[inBatch_sf$inBatch]),
                                                             st_coordinates(inBatch_sf$geometry[pxl,]), measure = "cheap"))/1000), dist<15)

  if(nrow(dst)>=4) {

  # plot(inBatch_sf$geometry, pch = 16, cex = 0.6, col = "grey90")
  # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch),], pch = 16, col = "orange", add = T)
  # plot(inBatch_sf$geometry[which(inBatch_sf$onLand),], pch = 16, cex = 0.5, col = "darkgreen", add = T)
  # plot(inBatch_sf$geometry[pxl,], pch = 16, col = "red", add =T, cex = 0.85)
  # plot(inBatch_sf$geometry[which(inBatch_sf$inBatch)[dst[,1]],], pch = "x", add =T, cex = 0.85)

  weigth <- approx(c(0, max(dnorm(dst[,2], 0, 4))), c(0,1), dnorm(dst[,2], 0, 4), rule = 2)$y
  
  dat    <- outBatch$dat[which(inBatch_sf$inBatch)[dst[,1]],,]
  evi    <- ifelse(dat[,,1] <= 0 , NA, dat[,,1])

  med <- apply(evi, 2, median, na.rm = T)
  sno <- ifelse(apply(dat[,,2], 2, median, na.rm = T)<1 | is.na(apply(dat[,,2], 2, median, na.rm = T)), TRUE, FALSE)
  
  if(any(!sno)) {
    ind <- which(diff(sno)!=0)

    spl <- split(sno, cut(1:length(sno), unique(c(0, ind, length(sno)+1)), labels = FALSE))
    msk <- unlist(sapply(spl, function(x) {
      if(length(x)<6) {
        x[] <- TRUE
      } else {
        if(all(!x)) x[unique(c(1:3, (length(x)-2):length(x)))] <- TRUE
      }
      x
    }))
    
    evi <- t(mapply(function(x) {
      if(any(!is.na(ifelse(msk, evi[x,], NA)))) {
        na.approx(ifelse(msk, evi[x,], NA), rule = 2)
      } else evi[x,]
    }, x = 1:nrow(evi)))
  } else {
    msk <- rep(TRUE, ncol(evi))
  }

  
  if(sum(is.na(med[msk]))<length(med[msk])*0.15) {
    
    wt   <- wt(cbind(1:length(med), na.approx(ifelse(!msk, 0, med), rule= 2)))
    pwL  <- apply(log2(abs(wt$power/wt$sigma2)), 1, median, na.rm = T)
    sig  <- apply(wt$signif, 1, quantile, probs = 0.75, na.rm = T)
    
    if(any(sig>=1)) {
    
    wt.pks <- FindPeaks(pwL)
    wt.s   <- cbind(pwL[wt.pks], round(wt$period[wt.pks],0), (sig>=1)[wt.pks])[order(pwL[wt.pks], decreasing = T),]
    wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[1:2,]  
    
    datCurve <- merge(data.frame(year = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%Y")),
                                 week = as.numeric(format(seq(min(dates), max(dates), by = "week"), "%U"))),
                      data.frame(year = as.numeric(format(dates, "%Y")), week = as.numeric(format(dates, "%U")),
                                 date = dates, id = 1:length(dates),
                                 evi  = na.approx(med, rule = 2), mask = msk, t(evi)), all.x = T)
    
    if(any(wt.sig[,3]==1) & any(wt.sig[,2]%in%c(45:60))) {
      
      fit0  <- optim(fn = lsCos, par = c(a = 1, b = 0), f = 52, Mx = na.approx(datCurve$evi - median(med, na.rm = T), rule =  3), sd = 0.01)
      curve <- fit0$par[1]*cos(pi*((1:nrow(datCurve))/(nrow(datCurve)/((nrow(datCurve)/52)*2))) +
                                 (pi+fit0$par[2])) +  mean(med, na.rm=T)
      
      mins <- sort(unique(c(1, FindPeaks(-curve), length(curve))))
      maxs <- which(diff(sign(diff(curve)))==-2)+1
      
    } else {
      mins <- maxs <- NA
    }
    } else {
      mins <- maxs <- NA
    }
    
    # with(datCurve[1:500,], plot(date, evi, type = "o"))
    # abline(v = datCurve$date[mins])
    
    if(diff(quantile(med[msk], prob = c(0.025,0.975), na.rm = T))>0.75 & length(mins)>20 & length(maxs)>20) {
      
      segL <- apply(do.call("rbind", lapply(maxs, function(x) mins[c(1:length(mins))[order(abs(x-mins))][1:2]]+c(-10,10))), 1,
                    function(x) datCurve[ifelse(x[1]<1, 1, x[1]):ifelse(x[2]>length(curve), length(curve), x[2]),])
      
      q50 <- median(med[msk], na.rm = T)
      
      phen0 <- do.call("rbind", lapply(segL, function(s) {
        
        # s <- segL[[2]]
        
        dateSeg <- na.approx(s$date, rule = 2)
        year    <- median(as.numeric(format(s$date, "%Y")), na.rm = T)
        
        # plot(rep(dateSeg, ncol(s)-6), unlist(c(s[,-c(1:6)])), pch = 16, cex = 0.5)
        
        if(length(dateSeg)>20 & sum(is.na(s[s$mask,-c(1:6)]))<length(unlist(c(s[s$mask,-c(1:6)])))*0.4) {
          
          pwL  <- apply(log2(abs(wt$power/wt$sigma2))[,s$id], 1, median, na.rm = T)
          sig  <- apply(wt$signif[,s$id], 1, median, na.rm = T)
          
          wt.pks <- FindPeaks(pwL)
          wt.s   <- cbind(pwL[wt.pks], wt$period[wt.pks], (sig>=1)[wt.pks])
          wt.sig <- rbind(wt.s[wt.s[,3]==1,], matrix(0, ncol = 3, nrow = 3))[order(c(wt.s[wt.s[,3]==1,1], rep(0,3)), decreasing = T),][1:2,]
          
          tmpDat <- subset(data.frame(x = rep(dateSeg, ncol(s)-6),
                                      m = rep(s$mask, ncol(s)-6),
                                      y = unlist(c(s[,-c(1:6)])),
                                      w = rep(weigth, each = nrow(s))), !is.na(y))
  
          
          dtsSm   <- seq(min(dateSeg), max(dateSeg), by = 24*60*60)
          
          spl     <- smooth.spline(x = tmpDat$x, y = tmpDat$y, spar = 0.3, w = tmpDat$w)
          XSeg    <- predict(spl, dateSeg)$y
          xSmooth <- predict(spl, dtsSm)$y
          
          # lines(dtsSm, xSmooth, type= "l", lwd = 6, col = "orange")
          
          peaks   <- FindPeaks(XSeg)
            peaks <- peaks[peaks%in%which(s$mask)]
          
          pars    <- list(rel_amp_frac = 0.15, rel_peak_frac = NULL, min_seg_amplitude = 0.075)
          segs0   <- tryCatch(do.call("rbind", GetSegs(peaks, ifelse(s$mask, XSeg, -1), pars)), error = function(e) NULL)
          
        
          if(!is.null(segs0)) {
            
            segs <- t(apply(segs0, 1, function(x) {
              sapply(dateSeg[x], function(y) which.min(abs(y-as.numeric(dtsSm))))
            }))
            
            if(nrow(segs)>1) {
              
              max_in   <- min(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
              max_out  <- max(segs[order(xSmooth[segs[,2]], decreasing = T),2][1:2])
              
              seqRan   <- c(min(segs[order(xSmooth[segs[,2]], decreasing = T),1]),
                            max(segs[order(xSmooth[segs[,2]], decreasing = T),3]))
              
            } else {
              max_in <- max_out  <- segs[2]
              seqRan <- segs[c(1,3)]
            }
            
            maxDate <- dtsSm[max_in]
            max     <- as.numeric(format(as.POSIXct(dtsSm, origin = "1970-01-01")[max_in], "%j"))
            amp    <- diff(range(xSmooth[seqRan[1]:seqRan[2]]))
            
            q10gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
            q50gup <- suppressWarnings(with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(q50, rev(y), gup = F, first_greater = T)]))
            q90gup <- with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = F, first_greater = T)])
              
            q90sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(max(y) - diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
            q50sen <- suppressWarnings(with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(q50, rev(y), gup = T, first_greater = F)]))
            q10sen <- with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  rev(t)[GetThresh(min(y) + diff(range(y))*0.1, rev(y), gup = T, first_greater = F)])
              
            area   <- with(data.frame(t = 1:length(dtsSm), y = xSmooth)[seqRan[1]:seqRan[2],],
                           MESS::auc(t, y, type = 'spline'))
            
            
            # plot(tmpDat$x, tmpDat$y, pch = 16, cex = 0.5)
            # lines(dtsSm, xSmooth, type= "l", lwd = 6, col = "orange")
            # abline(v = dtsSm[max_in], lty = 3, col = "red", lwd = 3)
            # apply(segs, 1, function(x) rect(dtsSm[x[1]], 0, dtsSm[x[3]], 1, col = adjustcolor("grey80", alpha.f = 0.2), border =  "grey10"))
            # points(q10gup, with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  min(y) + diff(range(y))*0.1), pch = 23, bg = "white", lwd = 2, cex = 3)
            # points(q50gup, q50, pch = 21, bg = "white", lwd = 2, cex = 3)
            # points(q90gup, with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  max(y) - diff(range(y))*0.1), pch = 22, bg = "white", lwd = 2, cex = 3)
            # points(q90sen, with(data.frame(t = dtsSm, y = xSmooth)[seqRan[1]:max_in,],  max(y) - diff(range(y))*0.1), pch = 23, bg = "cornflowerblue", lwd = 2, cex = 3)
            # points(q50sen, q50, pch = 21, bg = "cornflowerblue", lwd = 2, cex = 3)
            # points(q10sen,with(data.frame(t = dtsSm, y = xSmooth)[max_out:seqRan[2],],  min(y) + diff(range(y))*0.1), pch = 22, bg = "cornflowerblue", lwd = 2, cex = 3)

            out <- c(year = year,                                  #1
                     per1 = ifelse(wt.sig[1,3], wt.sig[1,2], NA),  #2
                     sig1 = ifelse(wt.sig[1,3], wt.sig[1,1], NA),  #3
                     per2 = ifelse(wt.sig[2,3], wt.sig[2,2], NA),  #4
                     per2 = ifelse(wt.sig[2,3], wt.sig[2,1], NA),  #5
                     max = max,                                    #6
                     amp = round(amp, 2),                          #7
                     area = round(area,2),                         #8
                     as.POSIXlt(maxDate, origin = "1970-01-01")$yday + (c(q10gup, q50gup, q90gup, q90sen, q50sen, q10sen) -  maxDate)/24/60/60) #9, 10, 11, 12, 13, 14
            
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
      
      wt.sig <- matrix(0, ncol = 3, nrow = 2)
      
      phen0 =  cbind(year = 1981:2020,                      
                     per1 = NA,
                     sig1 = NA,
                     per2 = NA,
                     sig2 = NA,
                     max  = NA,
                     amp  = NA,
                     area = NA, matrix(ncol = 6, nrow = length(1981:2020)))
    }
    
  } else {
    
    wt.sig <- matrix(0, ncol = 3, nrow = 2)
    
    phen0 =  cbind(year = 1981:2020,                      
                   per1 = NA,
                   sig1 = NA,
                   per2 = NA,
                   sig2 = NA,
                   max  = NA,
                   amp  = NA,
                   area = NA, matrix(NA, ncol = 6, nrow = length(1981:2020)))
  }
  
  } else {
    
    wt.sig <- matrix(0, ncol = 3, nrow = 2)
    
    phen0 =  cbind(year = 1981:2020,                      
                   per1 = NA,
                   sig1 = NA,
                   per2 = NA,
                   sig2 = NA,
                   max  = NA,
                   amp  = NA,
                   area = NA, matrix(NA, ncol = 6, nrow = length(1981:2020)))
    
  }
  

  phen <- aggregate(phen0, by = list(phen0[,1]), median)[,-1]
  merge(data.frame(year = 1981:2020), as.data.frame(phen), all.x = T)[,-1]
  
}
