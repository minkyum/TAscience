rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(RColorBrewer)
library(viridis)

##############################
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))
# tile <- 'h10v03'; year <- 2012


sinu_crs <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

mcd12q1_path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01'
sstr <- paste('*A*.',tile,'.006*.hdf',sep='')
file <- list.files(path=mcd12q1_path,pattern=glob2rx(sstr),full.names=T)
sds <- get_subdatasets(file)
lct <- raster(sds[1])

## panel grid
cxx <- 10
cyy <- 10

panel_grid <- matrix(1:length(lct),dim(lct)[1]/cxx,dim(lct)[2]/cyy)
panel_grid <- raster(panel_grid,crs=sinu_crs,
                     xmn=lct@extent@xmin,xmx=lct@extent@xmax,
                     ymn=lct@extent@ymin,ymx=lct@extent@ymax)


####
path <- paste('/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/MOLA/MYD21A1D.006')
sstr <- paste('*A',year,'*.',tile,'.006*.hdf',sep='')
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T,recursive=T)

# if(vv==1){
#   file <- file
# }else if(vv==2){
#   fileSpr <- file[as.numeric(substr(file,91,92)) > 2 & as.numeric(substr(file,91,92)) < 6]
#   file <- fileSpr
# }else{
#   fileSum <- file[as.numeric(substr(file,91,92)) > 6 & as.numeric(substr(file,91,92)) < 10]
#   file <- fileSum
# }

# print(var[vv])
print(length(file))

#
lstMat <- matrix(NA,240*240,length(file))
for(i in 1:length(file)){
  sds <- get_subdatasets(file[i])
  rast <- raster(sds[1])
  rast <- projectRaster(rast,panel_grid)

  lstMat[,i] <- values(rast)

  if(i%%10==0) print(i)
}

#########
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/met/')
save(lstMat,file,panel_grid,
     file=paste0('lst_',tile,'_',year,'.rda'))

