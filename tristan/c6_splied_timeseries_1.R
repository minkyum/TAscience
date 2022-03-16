library(sp)
library(raster)
library(gdalUtils)
library(rgdal)
library(doMC)
library(doParallel)


args <- commandArgs()
print(args)

xx <- as.numeric(args[3])


tile = 'h10v03'
year = c(2019, 2020, 2021)
mcd12a4_path <- paste('/projectnb/modislc/projects/sat/data/e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006')
search_str <- paste('*43A4.A',year[1],'*.',tile,'.006*',sep='')
files <- list.files(path=mcd12a4_path, pattern=glob2rx(search_str),full.names=T,include.dirs=F,recursive=T)

test <- raster(get_subdatasets(files[1])[[1]])

x <- -127.6731 + xx*0.005
y <- 52.54944
points <- cbind(x,y)
v <- SpatialPoints(points, proj4string = CRS("+proj=longlat +datum=WGS84"))
y <- spTransform(v, CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"))

numPix <- length(test)
imgNum <- setValues(test, 1:numPix)
z <- extract(imgNum, y)
sam <- z #make lat/lon coord sample((2400*2400),30)


## Extract splined EVI 
ncol <- 2400
nrow <- 2400
nbands <- 365
cnt <- ncol*nrow*nbands

data <- readBin('/projectnb/modislc/projects/sat/data/spline/h10v03/tg_v1.h10v03.2019.evi2.bip',
                what="integer",n=cnt,size=2,endian="little")
data <- array(data,c(nbands, ncol, nrow))
data <- aperm(data, c(3,2,1)) #for transposing
data <- brick(data)

splined <- matrix(NA,365,1)
for(i in 1:365){
  rast <- data[[i]]
  
  splined[i,] <- rast[sam]
}

registerDoMC()
evi_mm <- matrix(NA,365,1)
for(i in 1:365){
# evi_mm <- foreach(i=1:365,.combine=rbind) %dopar% {
  
  sds <- get_subdatasets(files[i])
  
  b1 <- raster(sds[[8]])[sam]
  b2 <- raster(sds[[9]])[sam]
  
  evi_mm[i,1] <- 2.5*((b2-b1)/(b2 + (2.4*b1) +1)) # calculate EVI2 

  # print(i)
}


setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v9/ts/')
png(file=paste('c6_v9_',xx,'.png',sep=''),res=300,units='in',width=13.5,height=6.5)
par(oma=c(2,2,2,2),mar=c(3,3,1,0))
plot(splined[,1]/10000,ylim=c(0,0.8),
     xlab='Date',ylab='EVI2',lwd=1.5,
     cex.lab=1.3,type='l')
points(evi_mm,col='red',cex=1.5)
dev.off()
