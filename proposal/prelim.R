
#####################################################
# Get NBAR (MCD43A4 C6) data from USGS
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
year <- as.numeric(substr(args[3],7,10))

# MCD12Q1 & MCD12Q2
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/')

url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/',year,'.01.01/',sep='')
system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD12Q1.A*.',tile,'.006.*.hdf" ',url,sep=''))    

url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q2.006/',year,'.01.01/',sep='')
system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD12Q2.A*.',tile,'.006.*.hdf" ',url,sep=''))    


# MYD21A1D & MYD21A1N
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/')
mod_dates <- seq(1,366)
for(i in 1:365){
  d = mod_dates[i]
  d_str = sprintf("%03d",d)
  dates <- as.Date(d,format='%Y.%m.%d',origin=paste((year-1),'.12.31',sep=''))
  mm <- substr(dates,6,7)
  dd <- substr(dates,9,10)
  
  url <- paste('https://e4ftl01.cr.usgs.gov/MOLA/MYD21A1D.006/',year,'.',mm,'.',dd,'/',sep='')
  system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MYD21A1D.A',year,d_str,'.',tile,'.006.*.hdf" ',url,sep=''))    
  
  url <- paste('https://e4ftl01.cr.usgs.gov/MOLA/MYD21A1N.006/',year,'.',mm,'.',dd,'/',sep='')
  system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MYD21A1N.A',year,d_str,'.',tile,'.006.*.hdf" ',url,sep=''))    
  
  print(i)
}




#####################################################
library(sp)
library(raster)
library(gdalUtils)
library(rgdal)
library(RColorBrewer)


###
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
# tile <- 'h12v04'

###


path_q1 <- '/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006'
path_q2 <- '/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/MOTA/MCD12Q2.006'
path_ad <- '/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/MOLA/MYD21A1D.006'
path_an <- '/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/MOLA/MYD21A1N.006'

sstr <- paste0('*',tile,'*')
files <- list.files(path_q2,pattern=glob2rx(sstr),recursive=T,full.names=T)
sds <- get_subdatasets(files[1])
baseImg <- raster(sds[1])


###
#Create output directory if it doesn't exist
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/stat/'
outDir <- paste0(path,tile,'/')
if (!dir.exists(outDir)) {dir.create(outDir)}

funmat <- function(vv){
  dat <- matrix(NA,(2400*2400),length(files))
  
  for(i in 1:length(files)){
    sds <- get_subdatasets(files[i])
    rast <- raster(sds[vv])
    rast[rast>32760] <- NA
    doy_offset <- as.integer(as.Date(paste((1999+i),'-12-31',sep='')) - as.Date("1970-1-1"))
    dat[,i] <- values(rast) - doy_offset
  }
  return(dat)
}

dat1 <- funmat(3)
rast1 <- setValues(baseImg,apply(dat1[,1:9],1,median,na.rm=T))
rast2 <- setValues(baseImg,apply(dat1[,10:18],1,median,na.rm=T))
Rast <- rast2-rast1
writeRaster(Rast,filename=paste0(outDir,'diff_midonset.tif'),format="GTiff",overwrite=TRUE)  

dat2 <- funmat(7)
rast1 <- setValues(baseImg,apply(dat2[,1:9],1,median,na.rm=T))
rast2 <- setValues(baseImg,apply(dat2[,10:18],1,median,na.rm=T))
Rast <- rast2-rast1
writeRaster(Rast,filename=paste0(outDir,'diff_middown.tif'),format="GTiff",overwrite=TRUE)  

dat <- dat2 - dat1
rast1 <- setValues(baseImg,apply(dat[,1:9],1,median,na.rm=T))
rast2 <- setValues(baseImg,apply(dat[,10:18],1,median,na.rm=T))
Rast <- rast2-rast1
writeRaster(Rast,filename=paste0(outDir,'diff_gsl.tif'),format="GTiff",overwrite=TRUE)  

dat <- funmat(5) - funmat(2)
rast1 <- setValues(baseImg,apply(dat[,1:9],1,median,na.rm=T))
rast2 <- setValues(baseImg,apply(dat[,10:18],1,median,na.rm=T))
Rast <- rast2-rast1
writeRaster(Rast,filename=paste0(outDir,'diff_rateonset.tif'),format="GTiff",overwrite=TRUE)  

dat <- funmat(8) - funmat(6)
rast1 <- setValues(baseImg,apply(dat[,1:9],1,median,na.rm=T))
rast2 <- setValues(baseImg,apply(dat[,10:18],1,median,na.rm=T))
Rast <- rast2-rast1
writeRaster(Rast,filename=paste0(outDir,'diff_ratedown.tif'),format="GTiff",overwrite=TRUE)  





###
sstr <- paste0('*',tile,'*')
files <- list.files(path_ad,pattern=glob2rx(sstr),recursive=T,full.names=T)

funlst <- function(files){
  dat <- vector('list',length(files))
  for(i in 1:length(files)){
    sds <- get_subdatasets(files[i])
    dat[[i]] <- raster(sds[1])
  }
  aa <- stack(dat)
  aa <- calc(aa,median,na.rm=T)
  return(aa)
}

for(i in 1:12){
  mm <- sprintf('%02d',i)
  
  files1 <- files[substr(files,91,92)==mm & as.numeric(substr(files,86,89))<2010]
  rast1 <- funlst(files1)
  writeRaster(rast1,filename=paste0(outDir,'lst_',mm,'_1.tif'),format="GTiff",overwrite=TRUE)  
  
  files2 <- files[substr(files,91,92)==mm & as.numeric(substr(files,86,89))>2009]
  rast2 <- funlst(files2)
  writeRaster(rast2,filename=paste0(outDir,'lst_',mm,'_2.tif'),format="GTiff",overwrite=TRUE)  
  
  print(i)
}





#####################################################
tile <- 'h12v04'

path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/stat/'
outDir <- paste0(path,tile,'/')

fileLst <- list.files(outDir,pattern=glob2rx('*lst*'),full.names=T)

setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
png(filename='lst.png',width=13,height=8,units='in',res=300)
par(mfrow=c(3,4),oma=c(0,0.5,1,0),mar=c(0,0,0,0))
for(i in 1:12){
  rast01 <- raster(fileLst[2*i]) - raster(fileLst[2*(i-1)+1])  
  rast01[rast01 >  5] <- 5
  rast01[rast01 < -5] <- -5
  plot(rast01,zlim=c(-5,5),col=rev(brewer.pal(11,'RdYlBu')),
       bty="n",axes=F,box=F,legend=F)
}
dev.off()


###
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/stat/'
outDir <- paste0(path,tile,'/')

fileLst <- list.files(outDir,pattern=glob2rx('*dif*'),full.names=T)

rast01 <- raster(fileLst[3]) 
rast01 <- focal(rast01,w=matrix(1,3,3),fun=median,na.rm=T)
rast01[rast01 >  5] <- 5
rast01[rast01 < -5] <- -5

rast02 <- raster(fileLst[2]) 
rast02 <- focal(rast02,w=matrix(1,3,3),fun=median,na.rm=T)
rast02[rast02 >  5] <- 5
rast02[rast02 < -5] <- -5

rast03 <- raster(fileLst[1]) 
rast03 <- focal(rast03,w=matrix(1,3,3),fun=median,na.rm=T)
rast03[rast03 >  5] <- 5
rast03[rast03 < -5] <- -5

rast04 <- raster(fileLst[5]) 
rast04 <- focal(rast04,w=matrix(1,3,3),fun=median,na.rm=T)
rast04[rast04 >  5] <- 5
rast04[rast04 < -5] <- -5

rast05 <- raster(fileLst[4]) 
rast05 <- focal(rast05,w=matrix(1,3,3),fun=median,na.rm=T)
rast05[rast05 >  5] <- 5
rast05[rast05 < -5] <- -5


setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
png(filename='metrics.png',width=13,height=8,units='in',res=300)
par(mfrow=c(2,3),oma=c(0,0.5,1,0),mar=c(0,0,0,0))
plot(rast01,zlim=c(-5,5),col=(brewer.pal(11,'Spectral')),
     bty="n",axes=F,box=F,legend=F,colNA='grey85')
plot(rast02,zlim=c(-5,5),col=rev(brewer.pal(11,'Spectral')),
     bty="n",axes=F,box=F,legend=F,colNA='grey85')
plot(rast03,zlim=c(-5,5),col=rev(brewer.pal(11,'Spectral')),
     bty="n",axes=F,box=F,legend=F,colNA='grey85')
plot(rast04,zlim=c(-5,5),col=(brewer.pal(11,'Spectral')),
     bty="n",axes=F,box=F,legend=F,colNA='grey85')
plot(rast05,zlim=c(-5,5),col=(brewer.pal(11,'Spectral')),
     bty="n",axes=F,box=F,legend=F,colNA='grey85')
dev.off()













