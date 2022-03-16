rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(RColorBrewer)
library(viridis)
library(doMC)

# ##############################
# # From bash code
# args <- commandArgs()
# print(args)
# 
# ss <- as.numeric(args[3])
# # ss <- 7

#############################################
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

sites <- sites[cord[,1]>0&igbp=='DBF']
cord <- cord[cord[,1]>0&igbp=='DBF',]

print(length(sites)==dim(cord)[1])

# ####
# # Get site coordiante
# source('~/R/Codes/latlong2MODIStile.R', echo=TRUE)
# 
# sam <- latlong2tile(cord[ss,1],cord[ss,2])
# hh <- sprintf('%02d',sam[1])
# vv <- sprintf('%02d',sam[2])
# tile <- paste0('h',hh,'v',vv)
# 
# print(tile)
# 
# mcd12q1_path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01'
# sstr <- paste0('*',tile,'*')
# file <- list.files(path=mcd12q1_path,pattern=glob2rx(sstr),full.names=T)
# 
# # if(length(file)==0){
# #   for(year in 2001:2019){
# #     # MCD12Q1
# #     setwd('/projectnb/modislc/users/mkmoon/mcd12q1/')
# #     url <- paste('https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/',year,'.01.01/',sep='')
# #     system(paste('wget --user=mkmoon --password=M159k258! -l1 -r --no-parent -A "MCD12Q1.A*.',tile,'.006.*.hdf" ',url,sep=''))
# #   }
# # }
# 
# sds <- get_subdatasets(file)
# lct <- raster(sds[1])
# 
# pixNum <- setValues(lct,1:length(lct))
# 
# 
# geog_crs = CRS("+proj=longlat +datum=WGS84")
# site <- data.frame(cord[ss,2],cord[ss,1])
# colnames(site) <- c('lon','lat')
# xy   <- site[,c(1,2)]
# bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
# bb   <- spTransform(bb,crs(pixNum))
# 
# pixN <- unlist(extract(pixNum,bb,buffer=1500))
# lctN <- values(lct)[pixN]
# 
# # ##############################
# # ncol <- 2400
# # nrow <- 2400
# # nbands <- 365
# # cnt <- ncol*nrow*nbands
# #
# # splined <- matrix(NA,(365*20),length(pixN))
# # for(j in 1:20){
# #   data <- readBin(paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/',
# #                         tile,'/',(2000+j),'/c6_tas.',tile,'.',(2000+j),'.evi2.bip',sep=''),
# #                   what="integer",n=cnt,size=2,endian="little")
# #   data <- array(data,c(nbands, ncol, nrow))
# #   data <- aperm(data, c(3,2,1)) #for transposing
# #   data <- brick(data)
# #   for(i in 1:365){
# #     rast <- data[[i]]
# #     splined[((365*(j-1))+i),] <- values(rast)[pixN]
# #   }
# #   print(j)
# # }
# #
# # save(splined,
# #      file=paste('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/splined_',ss,'.rda',sep=''))
# 
# 
# ##############################
# phemet1 <- matrix(NA,20,length(pixN))
# phemet2 <- matrix(NA,20,length(pixN))
# phemet3 <- matrix(NA,20,length(pixN))
# phemet4 <- matrix(NA,20,length(pixN))
# phemet5 <- matrix(NA,20,length(pixN))
# phemet6 <- matrix(NA,20,length(pixN))
# phemet7 <- matrix(NA,20,length(pixN))
# phemet8 <- matrix(NA,20,length(pixN))
# phemet9 <- matrix(NA,20,length(pixN))
# phemet0 <- matrix(NA,20,length(pixN))
# for(i in 1:15){
#   yy <- 2000 + i
#   doy_offset <- as.integer(as.Date(paste(yy, "-1-1", sep="")) - as.Date("1970-1-1"))
# 
#   if(i<19){
#     path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/product/',yy,sep='')
#     sstr <- paste('*',yy,'*',tile,'*.hdf',sep='')
#     file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
#     temp <- get_subdatasets(file)
# 
#     phemet1[i,] <- ifelse(values(raster(temp[2]))[pixN]>32000,NA, values(raster(temp[2]))[pixN] - doy_offset)
#     phemet2[i,] <- ifelse(values(raster(temp[3]))[pixN]>32000,NA, values(raster(temp[3]))[pixN] - doy_offset)
#     phemet3[i,] <- ifelse(values(raster(temp[4]))[pixN]>32000,NA, values(raster(temp[4]))[pixN] - doy_offset)
#     phemet4[i,] <- ifelse(values(raster(temp[5]))[pixN]>32000,NA, values(raster(temp[5]))[pixN] - doy_offset)
#     phemet5[i,] <- ifelse(values(raster(temp[6]))[pixN]>32000,NA, values(raster(temp[6]))[pixN] - doy_offset)
#     phemet6[i,] <- ifelse(values(raster(temp[7]))[pixN]>32000,NA, values(raster(temp[7]))[pixN] - doy_offset)
#     phemet7[i,] <- ifelse(values(raster(temp[8]))[pixN]>32000,NA, values(raster(temp[8]))[pixN] - doy_offset)
#     phemet8[i,] <- ifelse(values(raster(temp[9]))[pixN]>32000,NA, values(raster(temp[9]))[pixN])
#     phemet9[i,] <- ifelse(values(raster(temp[10]))[pixN]>32000,NA, values(raster(temp[10]))[pixN])
#     phemet0[i,] <- ifelse(values(raster(temp[11]))[pixN]>32000,NA, values(raster(temp[11]))[pixN])
# 
#   }else{
#     path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/pheno_out/',tile,sep='')
#     sstr <- paste('*',tile,'_',yy,sep='')
#     file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
# 
#     phemet1[i,] <- ifelse(values(raster(file,band=5))[pixN]>32000,NA, values(raster(file,band=5))[pixN] - doy_offset)
#     phemet2[i,] <- ifelse(values(raster(file,band=7))[pixN]>32000,NA, values(raster(file,band=7))[pixN] - doy_offset)
#     phemet3[i,] <- ifelse(values(raster(file,band=8))[pixN]>32000,NA, values(raster(file,band=8))[pixN] - doy_offset)
#     phemet4[i,] <- ifelse(values(raster(file,band=11))[pixN]>32000,NA, values(raster(file,band=11))[pixN] - doy_offset)
#   }
#   print(i)
# }
# 
# setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/ano_flux')
# save(phemet1,phemet2,phemet3,phemet4,phemet5,phemet6,phemet7,phemet8,phemet9,phemet0,
#      file=paste('ano_flux_g_',ss,'.rda',sep=''))


##############################

ss <- 1
par(mfrow=c(4,4))
for(ss in c(1:4,6,9:24)){
  # load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/splined_',ss,'.rda'))
  load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/ano_flux/ano_flux_g_',ss,'.rda'))


  # if(ss==1){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Ha1_FLUXNET2015_FULLSET_DD_1991-2012_1-4.csv')
  # }else if(ss==2){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-MMS_FLUXNET2015_FULLSET_DD_1999-2014_1-4.csv')
  # }else if(ss==3){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-Oho_FLUXNET2015_FULLSET_DD_2004-2013_1-4.csv')
  # }else if(ss==4){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-UMB_FLUXNET2015_FULLSET_DD_2000-2014_1-4.csv')
  # }else if(ss==5){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-UMd_FLUXNET2015_FULLSET_DD_2007-2014_1-4.csv')
  # }else if(ss==6){
  #   flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-WCr_FLUXNET2015_FULLSET_DD_1999-2014_1-4.csv')
  # }

  path <- '/projectnb/modislc/users/mkmoon/Fluxnet2015_newer'
  sstr <- paste('FLX_',sites[ss],'_FLUXNET2015_FULLSET_DD_*.csv',sep='')
  file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

  flux <- read.csv(file)

  # years <- as.numeric(unique(substr(flux$TIMESTAMP,1,4)))
  # years <- years[years>2000]

  ##############################
  anoMat <- matrix(NA,20,18)
  # anoMat[,20] <- sites[ss]
  # anoMat[,21] <- 2001:2020
  for(i in 1:20){
    xx <- c(median(phemet1[i,],na.rm=T),median(phemet2[i,],na.rm=T),median(phemet3[i,],na.rm=T),median(phemet4[i,],na.rm=T),median(phemet5[i,],na.rm=T),
            median(phemet6[i,],na.rm=T),median(phemet7[i,],na.rm=T),median(phemet8[i,],na.rm=T),median(phemet9[i,],na.rm=T),median(phemet0[i,],na.rm=T))
    # yy <- apply(splined[(365*(i-1)+1):(365*i),],1,median,na.rm=T)
    gg <- flux[as.numeric(substr(flux$TIMESTAMP,1,4))>(1999+i)&as.numeric(substr(flux$TIMESTAMP,1,4))<(2001+i),255]
    gg[gg < -9995] <- NA

    # anoMat[i,1] <- sum(yy[xx[1]:xx[4]])*0.01
    # anoMat[i,2] <- sum(yy[xx[1]:xx[2]])*0.01
    # anoMat[i,3] <- sum(yy[xx[2]:xx[3]])*0.01
    # anoMat[i,4] <- sum(yy[xx[3]:xx[4]])*0.01

    if(sum(!is.na(gg))>360){
      plot(gg)
      gg <- smooth.spline(gg,spar=0.55)$y
      anoMat[i,1] <- sum(gg[xx[2]:xx[6]])
      anoMat[i,2] <- sum(gg[xx[1]:xx[4]])
      anoMat[i,3] <- sum(gg[xx[4]:xx[5]])
      anoMat[i,4] <- sum(gg[xx[5]:xx[7]])
    }

    anoMat[i,5:14] <- xx

    anoMat[i,15] <- xx[6]-xx[2]
    anoMat[i,16] <- xx[4]-xx[1]
    anoMat[i,17] <- xx[5]-xx[4]
    anoMat[i,18] <- xx[7]-xx[5]


  }
  # anoMat[3,] <- NA
  anoltm <- apply(anoMat,2,mean,na.rm=T)

  # for(i in 1:18){
  #   anoMat[,i] <- anoMat[,i] - anoltm[i]
  # }

  if(ss==1){
    AnoMat <- anoMat
  }else{
    AnoMat <- rbind(AnoMat,anoMat)
  }
  print(ss)
}

# colnames(AnoMat) <- c('EVIann','EVIgre','EVImid','EVIsen',
#                       'GPPann','GPPgre','GPPmid','GPPsen',
#                       'GSLann','GSLgre','GSLmid','GSLsen',
#                       'PHEgup','PHEmat','PHEsen','PHEdor','EVIamp',NA)

AnoMat[is.na(AnoMat[,1])|AnoMat[,1]<50,] <- NA

AnoMat <- as.data.frame(AnoMat)
AnoMat <- cbind(AnoMat,rep(sites[c(1:4,6,9:24)],each=20),rep(2001:2020,21))
AnoMat <- na.omit(AnoMat)


colnames(AnoMat) <- c('GPPann','GPPgre','GPPmid','GPPsen',
                      'PHEgup','PHEmgp','PHEpek','PHEmat','PHEsen',
                      'PHEmsn','PHEdor','EVImin','EVIamp','EVIare',
                      'GSLann','GSLgre','GSLmid','GSLsen',
                      'site','year')



write.csv(AnoMat,file='/projectnb/modislc/users/mkmoon/TAscience/Data/flux_sites_ext.csv')



###################################


reg <- lm(AnoMat[,1]~as.matrix(AnoMat[,5:18]))
pre <- predict(reg)
summary(lm(AnoMat[,1]~as.matrix(AnoMat[,c(5:18)])))

par(oma=c(1,2,1,1),mar=c(4,5,3,3))
plot(AnoMat[,1],pre,pch=21,bg=AnoMat[,19],
     xlim=c(-500,2000),ylim=c(-500,2000),
     xlab='Observed',ylab='Predicted',
     cex=1.3,cex.lab=1.5,cex.axis=1.5)
abline(0,1,lty=5)
abline(v=0)
abline(h=0)

pcp <- prcomp(AnoMat[,5:19],center=T,scale=T)
x <- scale(AnoMat[,5:19], center = T, scale = T)
pcas <- as.matrix(x) %*% pcp$rotation
reg <- lm(AnoMat[,1]~pcas)
pre <- predict(reg)
summary(lm(AnoMat[,1]~pcas))
plot(AnoMat[,1],pre)
abline(0,1)

setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
png(filename='flux_mr.png',width=12,height=4.3,unit='in',res=300)

par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(5,5.5,1,1))
Col <- brewer.pal(6,'Set1')
plot(na.omit(AnoMat[,5]),reg$fitted.values,cex=2,pch=21,
     xlim=c(300,2200),ylim=c(300,2200),
     cex.axis=1.7,cex.lab=2,
     ylab='Estimated Annual GPP',
     xlab='Measured Annual GPP',
     bg=c(rep(Col[1],12),rep(Col[2],14),rep(Col[3],10),rep(Col[4],14),rep(Col[5],8),rep(Col[6],10)))
abline(0,1,lwd=2,lty=5)
mtext(expression(paste(italic(r)^2,' = 0.83')),side=1,line=-4.5,adj=0.8,cex=1.3)
mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-2.5,adj=0.83,cex=1.3)

dev.off()


# ####
# setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
# png(filename=paste0('evi_gpp.png'),width=13,height=4,unit='in',res=300)
# 
# Col <- brewer.pal(6,'Set1')
# par(mfrow=c(1,4),oma=c(1,1,1,0),mar=c(5,5,1,0.5),mgp=c(2.5,1,0))
# plot(AnoMat[,1],AnoMat[,5],cex=1.8,pch=21,
#      xlim=c(6000,14000),ylim=c(500,2500),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Annual EVIarea',
#      ylab=expression(paste('Annual GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,5]~(AnoMat[,1])),lwd=2)
# reg <- summary(lm(AnoMat[,5]~AnoMat[,1]))
# legend('topleft',c('US-Ha1','US-MMS','US-Oho','US-UMB','US-UMD','US-WCr'),
#        pch=21,pt.bg=Col,cex=1.6,bty='n',pt.cex=2)
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.20')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,2],AnoMat[,6],cex=1.8,pch=21,
#      xlim=c(1000,3500),ylim=c(50,500),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Greenup EVIarea',
#      ylab=expression(paste('Geenup GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,6]~(AnoMat[,2])),lwd=2)
# reg <- summary(lm(AnoMat[,6]~AnoMat[,2]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.03')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' = 0.150')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,3],AnoMat[,7],cex=1.8,pch=21,
#      xlim=c(1800,7500),ylim=c(100,1100),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Mid-season EVIarea',
#      ylab=expression(paste('Mid-season GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,7]~(AnoMat[,3])),lwd=2)
# reg <- summary(lm(AnoMat[,7]~AnoMat[,3]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.02')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' = 0.238')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,4],AnoMat[,8],cex=1.8,pch=21,
#      xlim=c(2500,8000),ylim=c(100,1000),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Senescence EVIarea',
#      ylab=expression(paste('Senescence GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,8]~(AnoMat[,4])),lwd=2)
# reg <- summary(lm(AnoMat[,8]~AnoMat[,4]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.41')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# dev.off()
# 
# 
# 
# 
# setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
# png(filename=paste0('evi_gpp_ano.png'),width=13,height=4,unit='in',res=300)
# 
# Col <- brewer.pal(6,'Set1')
# par(mfrow=c(1,4),oma=c(1,1,1,0),mar=c(5,5,1,0.5),mgp=c(2.5,1,0))
# plot(AnoMat[,1],AnoMat[,5],cex=1.8,pch=21,
#      xlim=c(-1600,1600),ylim=c(-500,500),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Annual EVIarea',
#      ylab=expression(paste('Annual GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,5]~(AnoMat[,1])),lwd=2)
# reg <- summary(lm(AnoMat[,5]~AnoMat[,1]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.01')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' = 0.41')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,2],AnoMat[,6],cex=1.8,pch=21,
#      xlim=c(-1600,1600),ylim=c(-400,400),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Greenup EVIarea',
#      ylab=expression(paste('Geenup GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,6]~(AnoMat[,2])),lwd=2)
# reg <- summary(lm(AnoMat[,6]~AnoMat[,2]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.09')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' = 0.013')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,3],AnoMat[,7],cex=1.8,pch=21,
#      xlim=c(-1600,1600),ylim=c(-400,400),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Mid-season EVIarea',
#      ylab=expression(paste('Mid-season GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,7]~(AnoMat[,3])),lwd=2)
# legend('topleft',c('US-Ha1','US-MMS','US-Oho','US-UMB','US-UMD','US-WCr'),
#        pch=21,pt.bg=Col,cex=1.5,bty='n',pt.cex=2)
# reg <- summary(lm(AnoMat[,7]~AnoMat[,3]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.38')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' < 0.001')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# plot(AnoMat[,4],AnoMat[,8],cex=1.8,pch=21,
#      xlim=c(-1600,1600),ylim=c(-400,400),
#      cex.axis=1.5,cex.lab=1.8,
#      xlab='Senescence EVIarea',
#      ylab=expression(paste('Senescence GPP (g',CO[2],' ',m^-2,')')),
#      bg=c(rep(Col[1],20),rep(Col[2],20),rep(Col[3],20),rep(Col[4],20),rep(Col[5],20),rep(Col[6],20)))
# abline(lm(AnoMat[,8]~(AnoMat[,4])),lwd=2)
# reg <- summary(lm(AnoMat[,8]~AnoMat[,4]))
# rseq <- formatC(reg$r.squared,2,format='fg',flag='#')
# mtext(expression(paste(italic(r)^2,' = 0.08')),side=1,line=-4.5,adj=0.8,cex=1.3)
# mtext(expression(paste(italic(p),' < 0.016')),side=1,line=-2.5,adj=0.83,cex=1.3)
# 
# dev.off()
# 
# 
# 
# plot(AnoMat[,9],AnoMat[,5])
# plot(AnoMat[,10],AnoMat[,5])
# plot(AnoMat[,11],AnoMat[,5])
# plot(AnoMat[,12],AnoMat[,5])
# 
# 
# png(filename=paste0('ano_evi_gpp_',ss,'.png'),width=10,height=7,unit='in',res=300)
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(5,5,1,1))
# for(i in 1:4){
#   j <- i+4
#   plot(AnoMat[,c(i,j)],cex=1.2,
#        # xlim=c(-1200,1200),ylim=c(-420,420),
#        cex.axis=1.2,cex.lab=1.5)
#   text(AnoMat[,i],AnoMat[,j],1:12,pos=4)
# 
#   abline(lm(AnoMat[,j]~AnoMat[,i]))
# 
#   reg <- summary(lm(AnoMat[,j]~AnoMat[,i]))
#   legend('bottomright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
# }
# dev.off()
# 
# png(filename='ano_gsl_gpp_1.png',width=10,height=7,unit='in',res=300)
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(5,5,1,1))
# for(i in 1:4){
#   i <- i+8
#   j <- i-4
#   plot(AnoMat[,c(i,j)],cex=1.2,
#        xlim=c(-30,30),ylim=c(-420,420),
#        cex.axis=1.2,cex.lab=1.5)
#   text(AnoMat[,i],AnoMat[,j],1:12,pos=4)
# 
#   abline(lm(AnoMat[,j]~AnoMat[,i]))
# 
#   reg <- summary(lm(AnoMat[,j]~AnoMat[,i]))
#   legend('bottomright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
# }
# dev.off()
# 
# png(filename='ano_gsl_gpptot_1.png',width=10,height=7,unit='in',res=300)
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(5,5,1,1))
# for(i in 1:4){
#   i <- i+8
#   j <- 5
#   plot(AnoMat[,c(i,j)],cex=1.2,
#        xlim=c(-30,30),ylim=c(-420,420),
#        cex.axis=1.2,cex.lab=1.5)
#   text(AnoMat[,i],AnoMat[,j],1:12,pos=4)
# 
#   abline(lm(AnoMat[,j]~AnoMat[,i]))
# 
#   reg <- summary(lm(AnoMat[,j]~AnoMat[,i]))
#   legend('bottomright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
# }
# dev.off()
# 
# png(filename='ano_phe_gpptot_1.png',width=10,height=7,unit='in',res=300)
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(5,5,1,1))
# for(i in 1:4){
#   i <- i+12
#   j <- 5
#   plot(AnoMat[,c(i,j)],cex=1.2,
#        xlim=c(-30,30),ylim=c(-420,420),
#        cex.axis=1.2,cex.lab=1.5)
#   text(AnoMat[,i],AnoMat[,j],1:12,pos=4)
# 
#   abline(lm(AnoMat[,j]~AnoMat[,i]))
# 
#   reg <- summary(lm(AnoMat[,j]~AnoMat[,i]))
#   legend('bottomright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
# }
# dev.off()
# 
# png(filename='ano_evi_gpptot_1.png',width=10,height=7,unit='in',res=300)
# par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(5,5,1,1))
# for(i in 1:4){
#   j <- 5
#   plot(AnoMat[,c(i,j)],cex=1.2,
#        xlim=c(-120000,120000),ylim=c(-420,420),
#        cex.axis=1.2,cex.lab=1.5)
#   text(AnoMat[,i],AnoMat[,j],1:12,pos=4)
# 
#   abline(lm(AnoMat[,j]~AnoMat[,i]))
# 
#   reg <- summary(lm(AnoMat[,j]~AnoMat[,i]))
#   legend('bottomright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
# }
# dev.off()
# 
# 
# 
# #################################
# ss <- 2
# 
# load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/splined_',ss,'.rda'))
# load(paste0('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/ano_flux_',ss,'.rda'))
# flux <- read.csv('/projectnb/modislc/users/mkmoon/Fluxnet2015_newer/FLX_US-MMS_FLUXNET2015_FULLSET_DD_1999-2014_1-4.csv')
# 
# evi <- matrix(NA,365,14)
# gpp <- matrix(NA,365,14)
# tep <- matrix(NA,365,14)
# for(i in 1:14){
#   xx <- c(median(phemet1[i,],na.rm=T),median(phemet2[i,],na.rm=T),median(phemet3[i,],na.rm=T),median(phemet4[i,],na.rm=T))
#   yy <- apply(splined[(365*(i-1)+1):(365*i),],1,median,na.rm=T)
#   gg <- flux[as.numeric(substr(flux$TIMESTAMP,1,4))>(1999+i)&as.numeric(substr(flux$TIMESTAMP,1,4))<(2001+i),255]
#   gg[gg < -9995] <- NA
#   tt <- flux[as.numeric(substr(flux$TIMESTAMP,1,4))>(1999+i)&as.numeric(substr(flux$TIMESTAMP,1,4))<(2001+i),2]
#   tt[tt < -9995] <- NA
#   
#   gg <- smooth.spline(gg,spar=0.7)$y
#   tt <- smooth.spline(tt,spar=0.7)$y
#   
#   evi[,i] <- yy[1:365]
#   gpp[,i] <- gg[1:365]
#   tep[,i] <- tt[1:365]
# }  
# 
# EVI <- apply(evi,1,mean)
# GPP <- apply(gpp,1,mean)
# TEP <- apply(tep,1,mean)
# 
# 
# 
# #######################
# # png(filename='ano_2012_us-mms.png',width=6,height=6,unit='in',res=300)
# # par(mfrow=c(3,1),oma=c(2,1,0,1),mar=c(0.2,5,1,0),mgp=c(2.8,1,0))
# 
# png(filename='ano_2012_us-mms.png',width=11.9,height=3.55,unit='in',res=300)
# par(mfrow=c(1,3),oma=c(2,1,0,0),mar=c(4,5,1,1),mgp=c(2.8,1,0))
# 
# 
# plot(EVI*0.0001,type='l',lwd=3,axe=F,ann=F,ylim=c(0.1,0.95))
# axis(2,seq(0,1,0.2),cex.axis=2)
# axis(1,at=seq(0,365,100),cex.axis=2)
# mtext('EVI2',2,2.7,cex=2)
# box()
# lines(evi[,12]*0.0001,col='red',lwd=5)
# legend('topleft',c('Mean','2012'),lwd=c(3,5),bty='n',cex=2.2,col=c('black','red'))
# 
# plot(GPP,type='l',lwd=3,axe=F,ann=F,ylim=c(0,8))
# axis(2,seq(0,10,2),cex.axis=2)
# axis(1,at=seq(0,365,100),cex.axis=2)
# mtext('GPP',2,2.7,cex=2)
# box()
# lines(gpp[,12],col='red',lwd=5)
# # legend('topright',expression(paste(mu,'mol',CO[2],' ',m^-2,s^-1)),bty='n',cex=2.2)
# 
# 
# plot(TEP,type='l',lwd=3,axe=F,ann=F,ylim=c(-5,35))
# axis(2,seq(-10,50,10),cex.axis=2)
# axis(1,at=seq(0,365,100),cex.axis=2)
# mtext('Air Temp.',2,2.7,cex=2)
# box()
# lines(tep[,12],col='red',lwd=5)
# 
# mtext('Day of year',1,2.7,cex=1.7,line=-0.5)
# 
# dev.off()
# 
# 
# 
# 
# #############################
# flux <- read.csv('/projectnb/modislc/users/mkmoon/Download/AMF_US-KFS_BASE_HH_7-5.csv',skip=2)
# gg <-dat$FC_PI_F_1_1_1

