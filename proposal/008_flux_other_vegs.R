rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(RColorBrewer)
library(viridis)
library(doMC)
library(tidyverse)


smet <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_AA-Flx_BIF_YY_20200501.csv')

sites <- as.character(unique(smet$SITE_ID))

cord <- matrix(NA,length(sites),2)
igbp <- matrix(NA,length(sites),1)
for(i in 1:length(sites)){
  subs <- smet[smet$SITE_ID==sites[i],]
  
  cord[i,1] <- as.numeric(as.character(subs[subs$VARIABLE=='LOCATION_LAT',5]))[1]
  cord[i,2] <- as.numeric(as.character(subs[subs$VARIABLE=='LOCATION_LONG',5]))[1]
  igbp[i] <- as.character(subs[subs$VARIABLE=='IGBP',5])[1]
}














# ##############################
# tran1 <- matrix(NA,14,6)
# tran2 <- matrix(NA,14,6)
# for(i in 1:6){
#   if(i==1){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-ARM_FLUXNET2015_FULLSET_DD_2003-2012_1-4.csv')
#   }else if(i==2){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-CRT_FLUXNET2015_FULLSET_DD_2011-2013_1-4.csv')
#   }else if(i==3){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne1_FLUXNET2015_FULLSET_DD_2001-2013_1-4.csv')
#   }else if(i==4){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne2_FLUXNET2015_FULLSET_DD_2001-2013_1-4.csv')
#   }else if(i==5){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne3_FLUXNET2015_FULLSET_DD_2001-2013_1-4.csv')
#   }else{
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Twt_FLUXNET2015_FULLSET_DD_2009-2014_1-4.csv')
#   }
# 
#   for(j in 1:14){
#     temp <- flux[as.numeric(substr(flux$TIMESTAMP,1,4))>(1999+j) & as.numeric(substr(flux$TIMESTAMP,1,4))<(2001+j),]
#     nee <- temp$NEE_CUT_MEAN
#     nee[nee== -9999] <- NA
#     if(sum(!is.na(nee))>364&sum(nee<0)>30){
#       neeSm <- smooth.spline(nee,spar=0.6)
#       plot(nee);points(neeSm$y);abline(h=0)
# 
#       midp <- which(neeSm$y==min(neeSm$y))
#       tran1[j,i] <- max(which(neeSm$y[1:midp] > 0))
#       tran2[j,i] <- min(which(neeSm$y[(midp+1):365] > 0))+(midp+1)
#     }
#   }
#   print(i)
# }
# tran1[tran1==Inf] <- NA
# tran1[tran1==-Inf] <- NA
# tran2[tran2==Inf] <- NA
# tran2[tran2==-Inf] <- NA
# 
# 
# carb1 <- matrix(NA,14,6)
# carb2 <- matrix(NA,14,6)
# for(i in 1:6){
#   if(i==1){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-ARM_FLUXNET2015_FULLSET_YY_2003-2012_1-4.csv')
#   }else if(i==2){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-CRT_FLUXNET2015_FULLSET_YY_2011-2013_1-4.csv')
#   }else if(i==3){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne1_FLUXNET2015_FULLSET_YY_2001-2013_1-4.csv')
#   }else if(i==4){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne2_FLUXNET2015_FULLSET_YY_2001-2013_1-4.csv')
#   }else if(i==5){
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ne3_FLUXNET2015_FULLSET_YY_2001-2013_1-4.csv')
#   }else{
#     flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Twt_FLUXNET2015_FULLSET_YY_2009-2014_1-4.csv')
#   }
# 
#   for(j in 1:14){
#     temp <- flux[as.numeric(substr(flux$TIMESTAMP,1,4))>(1999+j) & as.numeric(substr(flux$TIMESTAMP,1,4))<(2001+j),]
#     if(dim(temp)[1]>0){
#       nee <- temp$NEE_CUT_MEAN
#       gpp <- temp$GPP_DT_CUT_MEAN
# 
#       carb1[j,i] <- nee
#       carb2[j,i] <- gpp
#     }
#   }
# 
#   print(i)
# }
# carb1[is.na(tran1)] <- NA
# carb1[is.na(tran2)] <- NA
# carb2[is.na(tran1)] <- NA
# carb2[is.na(tran2)] <- NA

############################################################
# From bash code
args <- commandArgs()
print(args)

ss <- as.numeric(args[3])
# ss <- 1

# Get site coordiante
source('~/R/Codes/latlong2MODIStile.R', echo=TRUE)

cor<- c(36.6058,	-97.4888,
        41.628495,	-83.347086,
        41.1650600,	-96.4766400,
        41.1648700,	-96.4701,
        41.1796700,	-96.4396500,
        38.105533,	-121.6521)

sam <- latlong2tile(cor[2*(ss-1)+1],cor[2*ss])
hh <- sprintf('%02d',sam[1])
vv <- sprintf('%02d',sam[2])
tile <- paste0('h',hh,'v',vv)

print(tile)

mcd12q1_path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01'
sstr <- paste0('*',tile,'*')
file <- list.files(path=mcd12q1_path,pattern=glob2rx(sstr),full.names=T)
sds <- get_subdatasets(file)
lct <- raster(sds[1])

pixNum <- setValues(lct,1:length(lct))

geog_crs = CRS("+proj=longlat +datum=WGS84")
site <- data.frame(cor[2*ss],cor[2*(ss-1)+1])
colnames(site) <- c('lon','lat')
xy   <- site[,c(1,2)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
bb   <- spTransform(bb,crs(pixNum))

pixN <- unlist(extract(pixNum,bb,buffer=700))


##############################
phemet1 <- matrix(NA,14,length(pixN))
phemet2 <- matrix(NA,14,length(pixN))
phemet3 <- matrix(NA,14,length(pixN))
phemet4 <- matrix(NA,14,length(pixN))
phemet5 <- matrix(NA,14,length(pixN))
phemet6 <- matrix(NA,14,length(pixN))
for(i in 1:14){
  yy <- 2000 + i
  doy_offset <- as.integer(as.Date(paste(yy, "-1-1", sep="")) - as.Date("1970-1-1"))

  if(i<19){
    path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/product/',yy,sep='')
    sstr <- paste('*',yy,'*',tile,'*.hdf',sep='')
    file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
    temp <- get_subdatasets(file)

    phemet1[i,] <- values(raster(temp[2]))[pixN] - doy_offset
    phemet2[i,] <- values(raster(temp[3]))[pixN] - doy_offset
    phemet3[i,] <- values(raster(temp[7]))[pixN] - doy_offset
    phemet4[i,] <- values(raster(temp[8]))[pixN] - doy_offset
    phemet5[i,] <- values(raster(temp[10]))[pixN]
    phemet6[i,] <- values(raster(temp[11]))[pixN]

  }else{
    path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/pheno_out/',tile,sep='')
    sstr <- paste('*',tile,'_',yy,sep='')
    file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

    phemet1[i,] <- values(raster(file,band=5))[pixN] - doy_offset
    phemet2[i,] <- values(raster(file,band=6))[pixN] - doy_offset
    phemet3[i,] <- values(raster(file,band=10))[pixN] - doy_offset
    phemet4[i,] <- values(raster(file,band=11))[pixN] - doy_offset
    phemet5[i,] <- values(raster(file,band=3))[pixN]
    phemet6[i,] <- values(raster(file,band=2))[pixN]
  }
  print(i)
}

setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/')
save(phemet1,phemet2,phemet3,phemet4,phemet5,phemet6,
     file=paste('flux_ag_',ss,'.rda',sep=''))

setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
png(filename='flux_ag.png',width=12,height=4.3,unit='in',res=300)

par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(5,5.5,1,1))
col <- brewer.pal(6,'Set1')
for(i in 1:6){
  load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/flux_ag_',i,'.rda'))

  phe2 <- apply(phemet2,1,median,na.rm=T)

  if(i==1){
    plot(phe2,tran1[,i],xlim=c(30,250),ylim=c(30,250),
         xlab='MODIS SOS (DOY)',ylab='NEE Source to Sink (DOY)',
         cex.lab=2,cex.axis=1.7,cex=2,pch=21,bg=col[i])
    abline(0,1,lty=5)
    dat <- cbind(phe2,tran1[,i])
  }else{
    points(phe2,tran1[,i],cex=2,pch=21,bg=col[i])
    temp <- cbind(phe2,tran1[,i])
    dat <- rbind(dat,temp)
  }
}
dat[c(6,11),] <- NA
reg <- summary(lm(dat[,2]~dat[,1]))
abline(lm(dat[,2]~dat[,1]),lwd=2)
legend('topleft',c('US-ARM','US-CRT','US-Ne1','US-Ne2','US-Ne3','US-Twt'),
       pch=21,pt.bg=col,cex=1.7,bty='n',pt.cex=2)
rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
rmse <- round(sqrt(mean((dat[,1]-dat[,2])^2,na.rm=T)),1)
mtext(expression(paste(italic(r)^2,' = 0.67')),side=1,line=-5.5,adj=0.8,cex=1.3)
mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-3.5,adj=0.83,cex=1.3)
mtext(paste0('RMSE = ',rmse),side=1,line=-2,adj=0.91,cex=1.3)

for(i in 1:6){
  load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/flux_ag_',i,'.rda'))

  phe <- apply(phemet3,1,median,na.rm=T)

  if(i==1){
    plot(phe,tran2[,i],xlim=c(80,380),ylim=c(80,380),
         xlab='MODIS EOS (DOY)',ylab='NEE Sink to Source (DOY)',
         cex.lab=2,cex.axis=1.7,cex=2,pch=21,bg=col[i])
    abline(0,1,lty=5)
    dat <- cbind(phe,tran2[,i])
  }else{
    points(phe,tran2[,i],cex=2,pch=21,bg=col[i])
    temp <- cbind(phe,tran2[,i])
    dat <- rbind(dat,temp)
  }
}
dat[c(6,11),] <- NA
reg <- summary(lm(dat[,2]~dat[,1]))
abline(lm(dat[,2]~dat[,1]),lwd=2)
rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
rmse <- round(sqrt(mean((dat[,1]-dat[,2])^2,na.rm=T)),1)
mtext(expression(paste(italic(r)^2,' = 0.55')),side=1,line=-5.5,adj=0.75,cex=1.3)
mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-3.5,adj=0.78,cex=1.3)
mtext(paste0('RMSE = ',rmse),side=1,line=-2,adj=0.91,cex=1.3)


for(i in 1:6){
  load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/flux_ag_',i,'.rda'))

  phe1 <- apply(phemet2,1,median,na.rm=T)
  phe2 <- apply(phemet3,1,median,na.rm=T)
  phe3 <- apply(phemet5,1,median,na.rm=T)
  phe4 <- apply(phemet6,1,median,na.rm=T)

  phe <- (phe2-phe1)

  if(i==1){
    plot(phe,carb2[,i],xlim=c(30,200),ylim=c(500,2500),
         xlab='MODIS GSL (EOS - SOS)',
         ylab=expression(paste('Annual GPP (g',CO[2],' ',m^-2,')')),
         cex.lab=2,cex.axis=1.7,cex=2,pch=21,bg=col[i])
    dat <- cbind(phe,carb2[,i],phe1,phe2,phe3,phe4)
  }else{
    points(phe,carb2[,i],cex=2,pch=21,bg=col[i])
    temp <- cbind(phe,carb2[,i],phe1,phe2,phe3,phe4)
    dat <- rbind(dat,temp)
  }
}
dat[c(6,11),] <- NA
reg <- summary(lm(dat[,2]~dat[,1]*dat[,3]*dat[,4]*dat[,5]*dat[,6]))
reg <- summary(lm(dat[,2]~dat[,1]))
abline(lm(dat[,2]~dat[,1]),lwd=2)
rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
rmse <- round(sqrt(mean((dat[,1]-dat[,2])^2,na.rm=T)),1)
mtext(expression(paste(italic(r)^2,' = 0.079')),side=1,line=-3.5,adj=0.9,cex=1.3)
mtext(expression(paste(italic(p),' = 0.05')),side=1,line=-1.5,adj=0.85,cex=1.3)

dev.off()