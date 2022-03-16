rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(RColorBrewer)
library(viridis)

############################################################
# From bash code
args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,3)) # 1-378
pp <- as.numeric(substr(args[3],4,4)) # 2-8
# ss <- 1; pp <- 2

print(ss)
print(pp)

##############################
meta <- read.csv('/projectnb/modislc/users/mkmoon/TAscience/Data/sites_fluxnet2015_Tier1.csv')
    
# Get site coordiante
source('~/R/Codes/latlong2MODIStile.R', echo=TRUE)

sam <- latlong2tile(meta$lat[ss],meta$lon[ss])
hh <- sprintf('%02d',sam[1])
vv <- sprintf('%02d',sam[2])
tile <- paste0('h',hh,'v',vv)

print(tile)

path <- '/projectnb/modislc/users/mkmoon/LCD_C6/product/2001'
sstr <- paste0('*',tile,'*')
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
sds <- get_subdatasets(file)
lct <- raster(sds[1])

pixNum <- setValues(lct,1:length(lct))

geog_crs = CRS("+proj=longlat +datum=WGS84")
site <- data.frame(meta$lon[ss],meta$lat[ss])
colnames(site) <- c('lon','lat')
xy   <- site[,c(1,2)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
bb   <- spTransform(bb,crs(pixNum))

pixN1 <- unlist(extract(pixNum,bb))
pixN3 <- unlist(extract(pixNum,bb,buffer=750))
pixN5 <- unlist(extract(pixNum,bb,buffer=1250))

###############################################
p1b1 <- matrix(NA,15,1)
p3b3 <- matrix(NA,15,length(pixN3))
p5b5 <- matrix(NA,15,length(pixN5))
for(i in 1:15){
  yy <- 2000 + i
  # doy_offset <- as.integer(as.Date(paste(yy, "-1-1", sep="")) - as.Date("1970-1-1"))
  
  path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/product/',yy,sep='')
  sstr <- paste('*',yy,'*',tile,'*.hdf',sep='')
  file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  temp <- get_subdatasets(file)
  
  temp <- raster(temp[pp])
  # temp[temp > 32760] <- NA
  
  p1b1[i,] <- values(temp)[pixN1] #- doy_offset
  p3b3[i,] <- values(temp)[pixN3] #- doy_offset
  p5b5[i,] <- values(temp)[pixN5] #- doy_offset
 
  print(i)
}

####################
var <- c('15gu','50gu','peak','90gu','90gd','50gd','15gd')

outFolder <- '/projectnb/modislc/users/mkmoon/TAscience/Data/forJ/'
outDir <- paste0(outFolder,var[(pp-1)],'/')
if (!dir.exists(outDir)) {dir.create(outDir)}

setwd(outDir)
save(p1b1,p3b3,p5b5,
     file=paste(meta$siteID[ss],'_',var[(pp-1)],'.rda',sep=''))



# ########################
# path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/forJ/'
# sstr <- paste('*15gd.rda',sep='')
# file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T,recursive=T)
# load(file)


