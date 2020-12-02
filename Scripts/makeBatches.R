library(sf)
library(raster)
library(ncdf4)
library(rgdal)
library(parallel)
library(pbmcapply)

## land mask
lake <- read_sf("Data/GeoDat/ne_50m_lakes/ne_50m_lakes.shp")
land <- read_sf("Data/GeoDat/ne_50m_land//ne_50m_land.shp")
lbox <- st_bbox(land, crs = 4326) %>% st_as_sfc()

###############
## batch ######
###############

r0 <- raster(extent(st_bbox(land)[c(1,3,2,4)]), res = 5)
# plot(rasterToPolygons(r0))
# plot(land$geometry, add = T)
pol <- rasterToPolygons(r0)

###############
## VHP data ###
###############

pathVHP <- "/Volumes/bioing/user/slisovsk/VHP_SM_SMN/"

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

writeRaster(r0,'Results/phenRaster.tif',options=c('TFW=YES'))

##############
## Snow data ##
###############
dat <- nc_open("Data/nhsce_v01r01_19661004_20201102.nc")

## "days since 1966-10-03"
t <- ncvar_get(dat, "time")
z <- ncvar_get(dat, "snow_cover_extent")
ylat <- ncvar_get(dat, "longitude")
xlon <- ncvar_get(dat, "latitude")

### Projection and raster extent
lon <- raster(ylat)
lat <- raster(xlon)

start <- as.POSIXct("1966-10-03", "GMT")
snowD <- start + ((t-6)*24*60*60)

## prj
prj <- "+proj=stere +lat_0=90 +lon_0=-0 +a=6371200 +b=6371200 +units=m +no_defs"
xy  <- project(cbind(lon[], lat[]), prj) 

r_snow <- raster(xmn = min(xy[,1]), xmx = max(xy[,2]), ymn = min(xy[,1]), ymx = max(xy[,2]), 
                 ncol = ncol(lon), nrow = ncol(lon), crs = CRS(prj))

nc_close(dat)

####################
## create batches ##
####################

plot(land$geometry, col = adjustcolor("grey90", alpha.f = 0.5))
plot(rasterToPolygons(r0), add = T)

for(p in 1:length(pol)) {
  
  plot(pol[p,], add = T, col = adjustcolor("orange", alpha.f = 0.25))
  
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
    
    onLand   <- suppressMessages(st_intersects(st_crds, land, sparse = F))
    
    if(any(onLand)) {
      
      inLake   <- apply(suppressMessages(st_intersects(st_crds, lake$geometry, sparse = F)), 1, any)
      
      inPol    <- suppressMessages(st_intersects(st_crds, st_as_sfc(pol[p,]) %>% st_set_crs(4326), sparse = F))
      outBatch <- list(crds = data.frame(crds, indBatch = inPol, land = onLand, lake = inLake))
      
      ## loop over dates
      parOut <- pbmclapply(1:length(fls), function(f) {
        
        crp <- crop(raster(paste0(pathVHP, fls[f])), polE)
        
        if(any(crds[,2]<0)) {
          snow <- rasterize(xy, r_snow, field = t(z[,,which.min(abs(snowD - dates[f]))]))
          sOut <- extract(snow, project(crds, prj))
        } else sOut <- rep(NA, nrow(crds))
        
        array(c(crp[], ifelse(is.na(sOut) | sOut<1, NA, 1)), dim = c(nrow(crds),1,2))
        
      }, mc.cores = 5)
      
      outBatch$dat <- abind::abind(parOut, along = 2)
      save(outBatch, file = paste0("Results/Batches/Batch_", p, ".rda"))
      
    }
    
  }
  
  
}

# test <- rasterFromXYZ(cbind(outBatch$crds[,1:2], outBatch$dat[,100,1]))
# plot(test, breaks = seq(0,1, length = 100))