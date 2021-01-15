library(sf)
library(raster)
library(ncdf4)
library(rgdal)
library(parallel)
library(mcapply)

## land mask
lake <- read_sf("Data/GeoDat/ne_50m_lakes/ne_50m_lakes.shp") %>% st_union()
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp") %>% st_union()
lbox <- st_bbox(land, crs = 4326) %>% st_as_sfc()

###############
## batch ######
###############

rPol <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(rPol))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(rPol)

###############
## VHP data ###
###############

pathVHP <- "/bioing/user/slisovsk/VHP_SM_SMN/"

fls  <- list.files(pathVHP)
dts  <- cbind(as.numeric(substr(sapply(strsplit(fls, "[.]"), function(x) x[[5]]), 2, 5)),
              as.numeric(substr(sapply(strsplit(fls, "[.]"), function(x) x[[5]]), 6, 8)))

dSeq <- aggregate(seq(as.POSIXct("1980-01-01"), as.POSIXct("2020-12-31"), by = "day"), by = list(
  format(seq(as.POSIXct("1980-01-01"), as.POSIXct("2020-12-31"), by = "day"), "%Y"),
  format(seq(as.POSIXct("1980-01-01"), as.POSIXct("2020-12-31"), by = "day"), "%U")
), FUN = median)

dates <- dSeq$x[apply(dts, 1, function(x) which(x[1]==as.numeric(dSeq$Group.1) & x[2]==as.numeric(dSeq$Group.2)))]

## result Raster
r0 <- raster(paste0(pathVHP, fls[[1]]))
# plot(r0)
r0[] <- NA 
names(r0) <- "phenRaster_empty"

# writeRaster(r0,'Results/phenRaster.tif',options=c('TFW=YES'))

##############
## Snow data ##
###############

# fls.gz     <- list.files("/bioing/user/slisovsk/24km/", pattern = ".asc.gz", recursive = T,  full.names = T)[-c(1:693)]
# datesSnow  <- as.Date(as.POSIXct(unlist(lapply(strsplit(fls.gz, "ims"), function(x) strsplit(x[[2]], "_24km")))[c(TRUE, FALSE)], format = "%Y%j"))
# 
# weeks      <- cbind(as.numeric(format(datesSnow, "%Y")), as.numeric(format(datesSnow, "%U")))
# 
# weekStack <- do.call("stack", lapply(unique(weeks[,2]), function(x) {
# 
#   rS <- do.call("stack", mclapply(which(weeks[,2]==x), function(w) {
#     tab0 <- readLines(fls.gz[w])
#     ind  <- unlist(suppressWarnings(parallel::mclapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x))), mc.cores = 2)))
#     tab  <- tab0[-which(ind)]
# 
#     z = do.call("rbind", parallel::mclapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]), mc.cores = 2))
#     r0 <- raster(z[nrow(z):1,])
#     r0[] <- ifelse(r0[]==4, 1, NA)
#     r0
#   }, mc.cores = 10))
# 
#   tt <- rS[[1]]; tt[] <- apply(rS[], 1, function(x) sum(x, na.rm = T)/length(x))
#   extent(tt) <- c(-12126597.0, -12126597.0 + 1024*23684.997, -12126840.0, -12126597.0 + 1024*23684.997)
#   proj4string(tt) <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"
# 
#   tt
# }))
# 
# names(weekStack) <- paste0("V_", unique(weeks[,2]))
save(weekStack, file = "Results/weekStack_snow.rda")

load("Results/weekStack_snow.rda")

####################
## create batches ##
####################

plot(land, col = adjustcolor("grey90", alpha.f = 0.5), add = F)
plot(pol, add = T)


for(p in 1224:length(pol)) {

  plot(pol[p,], add = T, col = adjustcolor("orange", alpha.f = 0.25))
  
  if(!file.exists(paste0("/bioing/user/slisovsk/Batches/Batch_", p, ".rda"))) {
  
  ### VHP within batch extent
  tmp <- raster(paste0(pathVHP, fls[1]))
  polE <- st_as_sfc(pol[p,]) %>% st_set_crs(4326) %>% st_transform(3857) %>% st_buffer(25*1000) %>%
    st_intersection(lbox %>% st_transform(3857)) %>%
    st_transform(4326) %>% as("Spatial")
  
  if(!is.null(intersect(extent(tmp), pol[p,]))) {
    
    ### VHP on land
    r_out    <- crop(tmp, polE)
    crds     <- coordinates(r_out)
    st_crds  <- st_sfc(lapply(1:nrow(crds), function(x) st_point(crds[x,])), crs = 4326)
    
    onLand   <- as.vector(suppressMessages(st_intersects(st_crds, land, sparse = F)))
    
    # plot(st_crds, pch = 16, col = "grey80", cex = 0.7)
    # plot(st_crds[onLand], pch= 16, cex = 0.4, add = T)
    
    snow     <- extract(weekStack, st_coordinates(st_crds %>% st_transform(proj4string(weekStack))))
    
    if(any(onLand)) {
      
      inLake   <- as.vector(suppressMessages(st_intersects(st_crds, lake, sparse = F)))
      
      inPol    <- suppressMessages(st_intersects(st_crds, st_as_sfc(pol[p,]) %>% st_set_crs(4326), sparse = F))
      outBatch <- list(crds = data.frame(crds, indBatch = inPol, land = onLand, lake = inLake))
      
      ## loop over dates
      parOut <- do.call("cbind", mclapply(1:length(fls), function(f) {
        
        crp <- crop(raster(paste0(pathVHP, fls[f])), polE)
        
        # if(any(crds[,2]>0)) {
        #   snow <- rasterize(xy, r_snow, field = t(z[,,which.min(abs(snowD - dates[f]))]))
        #   sOut <- extract(snow, project(crds, prj))
        # } else sOut <- rep(NA, nrow(crds))
        
        crp[]
        
      }, mc.cores = 15))
      
      outEviSnow <- abind::abind(parOut, mapply(function(w) snow[,w], w = as.numeric(format(dates, "%U"))), along = 3)
      
      outBatch$dat <- outEviSnow
      save(outBatch, file = paste0("/bioing/user/slisovsk/Batches/Batch_", p, ".rda"))
      
    }
    
  }
  }
}

# test <- rasterFromXYZ(cbind(outBatch$crds[,1:2], outBatch$dat[,100,1]))
# plot(test, breaks = seq(0,1, length = 100))