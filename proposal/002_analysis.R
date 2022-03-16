rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(ppcor)
library(Kendall)
library(RColorBrewer)
library(scales)
library(maptools)
library(zyp)
library(Kendall)
library(scales)
library(agricolae)
library(profvis)
library(viridis)

###########################################################
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- 'met*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

# SOS long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # sos[sos[,21]<50,1:20] <- NA
  temp <- apply(sos[,1:5],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_sos <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # sos[sos[,21]<50,1:20] <- NA
  temp <- apply(sos[,16:20],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_sos <- do.call(mosaic,rast)

diff_sos <- ltm2_sos-ltm1_sos

# EOS long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # eos[eos[,21]<50,1:20] <- NA
  temp <- apply(eos[,1:5],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_eos <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # eos[eos[,21]<50,1:20] <- NA
  temp <- apply(eos[,16:20],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_eos <- do.call(mosaic,rast)

diff_eos <- ltm2_eos-ltm1_eos

# GSL long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # gsl1[gsl1[,21]<50,1:20] <- NA
  temp <- apply(gsl1[,1:5],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_gsl <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # gsl1[gsl1[,21]<50,1:20] <- NA
  temp <- apply(gsl1[,16:20],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_gsl <- do.call(mosaic,rast)

diff_gsl1 <- ltm2_gsl-ltm1_gsl



### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/domain.shp')
nest <-c('10')
erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]

rastList <- vector('list',9)
for(i in 1:9){
  if(i==1){
    rast <- ltm1_sos
  }else if(i==2){
    rast <- ltm2_sos
  }else if(i==3){
    rast <- diff_sos
  }else if(i==4){
    rast <- ltm1_eos
  }else if(i==5){
    rast <- ltm2_eos
  }else if(i==6){
    rast <- diff_eos
  }else if(i==7){
    rast <- ltm1_gsl
  }else if(i==8){
    rast <- ltm2_gsl
  }else{
    rast <- diff_gsl1
  }
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 5000
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  rastList[[i]] <- rast
  
  print(i)
}



#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(rastList[[1]]))

png(filename='ltm_diff_na.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- rastList[[1]]
temp[temp<90] <- 90
temp[temp>170] <- 170
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(90,170),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(100,120,140,160),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'SOS',cex=2.2,pos=4)
text(-5300000,-2400000,'2001-2005',cex=2,pos=4)
##
temp <- rastList[[4]]
temp[temp<220] <- 220
temp[temp>280] <- 280
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(220,280),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(220,240,260,280),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EOS',cex=2.2,pos=4)
text(-5300000,-2400000,'2001-2005',cex=2.2,pos=4)
##
temp <- rastList[[7]]
temp[temp<70] <- 70
temp[temp>170] <- 170
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(70,170),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(80,100,120,140,160),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'GSL',cex=2.2,pos=4)
text(-5300000,-2400000,'2001-2005',cex=2.2,pos=4)


##
temp <- rastList[[2]]
temp[temp<90] <- 90
temp[temp>170] <- 170
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(90,170),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(100,120,140,160),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'SOS',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)
##
temp <- rastList[[5]]
temp[temp<220] <- 220
temp[temp>280] <- 280
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(220,280),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(220,240,260,280),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EOS',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)
##
temp <- rastList[[8]]
temp[temp<70] <- 70
temp[temp>170] <- 170
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'Spectral'),
     zlim=c(70,170),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(80,100,120,140,160),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'GSL',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)


##
temp <- rastList[[3]]
temp[temp< -11] <- -11
temp[temp>  11] <- 11
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-11,11),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-10,-5,0,5,10),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'SOS',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[6]]
temp[temp< -11] <- -11
temp[temp>  11] <- 11
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-11,11),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-10,-5,0,5,10),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'EOS',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[9]]
temp[temp< -11] <- -11
temp[temp>  11] <- 11
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-11,11),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-10,-5,0,5,10),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'GSL',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

dev.off()



##################################
# Trends
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- '*3by3*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # sos[sos[,21]<50,1:20] <- NA
  temp <- matrix(NA,dim(sos)[1],1)
  for(j in 1:dim(sos)[1]){
    if(sum(!is.na(sos[j,1:20]))>9){
      x <- 1:20
      y <- sos[j,1:20]
      z <- zyp.sen(y~x)
      temp[j] <- z$coefficients[2] 
    }
  }
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
trd_sos <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # eos[eos[,21]<50,1:20] <- NA
  temp <- matrix(NA,dim(eos)[1],1)
  for(j in 1:dim(eos)[1]){
    if(sum(!is.na(eos[j,1:20]))>9){
      x <- 1:20
      y <- eos[j,1:20]
      z <- zyp.sen(y~x)
      temp[j] <- z$coefficients[2] 
    }
  }
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
trd_eos <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # gsl1[gsl1[,21]<50,1:20] <- NA
  temp <- matrix(NA,dim(gsl1)[1],1)
  for(j in 1:dim(gsl1)[1]){
    if(sum(!is.na(gsl1[j,1:20]))>9){
      x <- 1:20
      y <- gsl1[j,1:20]
      z <- zyp.sen(y~x)
      temp[j] <- z$coefficients[2] 
    }
  }
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
trd_gsl1 <- do.call(mosaic,rast)


### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na_ne.shp')
# nest <-c('10')
# erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]
# erepa_sinu <- spTransform(erepa,crs(rt))

trdList <- vector('list',3)
for(i in 1:3){
  if(i==1){
    rast <- trd_sos
  }else if(i==2){
    rast <- trd_eos
  }else{
    rast <- trd_gsl1
  }
  
  # rast <- focal(rast,w=matrix(1,3,3),fun=mean,na.rm=T)
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 1500
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  trdList[[i]] <- rast
  
  print(i)
}

setwd('/projectnb/modislc/users/mkmoon/TAscience/maps/')
writeRaster(trdList[[1]],filename="trd_sos.tif", format="GTiff", overwrite=TRUE)
writeRaster(trdList[[2]],filename="trd_eos.tif", format="GTiff", overwrite=TRUE)
writeRaster(trdList[[3]],filename="trd_gsl.tif", format="GTiff", overwrite=TRUE)


#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(trdList[[1]]))

png(filename='ltm_trd_na.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(2,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- trdList[[1]]
temp[temp< -0.8] <- -0.8
temp[temp>  0.8] <- 0.8
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-0.8,0.8),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.6),
     axis.args=list(at=c(-0.8,-0.4,0,+0.4,+0.8),
                    cex.axis=1.8,font=1)) 
text(-4200000,-2200000,'SOS trend',cex=2,pos=4)
text(-4200000,-2700000,'(days per year)',cex=1.8,pos=4)
##
temp <- trdList[[2]]
temp[temp< -0.8] <- -0.8
temp[temp>  0.8] <- 0.8
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-0.8,0.8),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.6),
     axis.args=list(at=c(-0.8,-0.4,0,+0.4,+0.8),
                    cex.axis=1.8,font=1)) 
text(-4200000,-2200000,'EOS trend',cex=2,pos=4)
text(-4200000,-2700000,'(days per year)',cex=1.8,pos=4)
##
temp <- trdList[[3]]
temp[temp< -0.8] <- -0.8
temp[temp>  0.8] <- 0.8
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-0.8,0.8),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.6),
     axis.args=list(at=c(-0.8,-0.4,0,+0.4,+0.8),
                    cex.axis=1.8,font=1)) 
text(-4200000,-2200000,'GSL trend',cex=2,pos=4)
text(-4200000,-2700000,'(days per year)',cex=1.8,pos=4)

par(fig=c(0,1,0,0.5),oma=c(2,2,2,2),mar=c(4,4,3,2),new=T)
k <- 1
for(i in 1:length(hfsp)){
  if(i==1){
    temp <- cbind(as.numeric(dat[dat$species==hfsp[1],1]),
                  as.numeric(dat[dat$species==hfsp[1],4]))
    if(length(unique(temp[,1]))>10){
      plot(temp,ylim=c(85,155),xlim=c(2001,2020),axes=F,ann=F,
           bg=col_vector[i],pch=21,cex=1.5)  
      
      axis(1,c(2001,2005,2010,2015,2020),cex.axis=1.8)
      axis(2,seq(90,200,20),cex.axis=1.5)
      box()
      mtext('Bud Break dates (DOY)',2,line=2.8,cex=1.3)
      x <- temp[,1]; y <- temp[,2]
      abline(lm(y~x),lty=5,col=col_vector[i],lwd=2)  
    }
  }else{
    temp <- cbind(as.numeric(dat[dat$species==hfsp[i],1]),
                  as.numeric(dat[dat$species==hfsp[i],4]))
    if(length(unique(temp[,1]))>10){
      points(temp,ylim=c(80,170),bg=col_vector[i],pch=21,cex=1.5)    
      x <- temp[,1]; y <- temp[,2]
      abline(lm(y~x),lty=5,col=col_vector[i],lwd=2)  
      k <- k+1
    }
  }
  print(k)
}
legend('bottomleft',bty='n','Each color (line) for each sepcies (total 17 species)',cex=1.8)
mtext('Bud break dates at Harvard Forest',1,line=-20,cex=1.5)

dev.off()




# LST
###########################################################
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/met'

# Long-term mean Ann
ltm1_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2003*'  
  }else if(yy==2){
    sstr <- 'lst*_2004*'  
  }else if(yy==3){
    sstr <- 'lst*_2005*'  
  }else if(yy==4){
    sstr <- 'lst*_2006*'  
  }else{
    sstr <- 'lst*_2007*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    # tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    # temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm1_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm1_LST <- stack(ltm1_lst)
ltm1_LST_ann <- calc(ltm1_LST,mean,na.rm=T)

ltm2_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2016*'  
  }else if(yy==2){
    sstr <- 'lst*_2017*'  
  }else if(yy==3){
    sstr <- 'lst*_2018*'  
  }else if(yy==4){
    sstr <- 'lst*_2019*'  
  }else{
    sstr <- 'lst*_2020*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    # tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    # temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm2_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm2_LST <- stack(ltm2_lst)
ltm2_LST_ann <- calc(ltm2_LST,mean,na.rm=T)

diff_LST_ann <- ltm2_LST_ann-ltm1_LST_ann

# Long-term mean Spring
ltm1_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2003*'  
  }else if(yy==2){
    sstr <- 'lst*_2004*'  
  }else if(yy==3){
    sstr <- 'lst*_2005*'  
  }else if(yy==4){
    sstr <- 'lst*_2006*'  
  }else{
    sstr <- 'lst*_2007*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    # temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm1_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm1_LST <- stack(ltm1_lst)
ltm1_LST_spr <- calc(ltm1_LST,mean,na.rm=T)

ltm2_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2016*'  
  }else if(yy==2){
    sstr <- 'lst*_2017*'  
  }else if(yy==3){
    sstr <- 'lst*_2018*'  
  }else if(yy==4){
    sstr <- 'lst*_2019*'  
  }else{
    sstr <- 'lst*_2020*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    # temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm2_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm2_LST <- stack(ltm2_lst)
ltm2_LST_spr <- calc(ltm2_LST,mean,na.rm=T)

diff_LST_spr <- ltm2_LST_spr-ltm1_LST_spr

# Long-term mean Summer
ltm1_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2003*'  
  }else if(yy==2){
    sstr <- 'lst*_2004*'  
  }else if(yy==3){
    sstr <- 'lst*_2005*'  
  }else if(yy==4){
    sstr <- 'lst*_2006*'  
  }else{
    sstr <- 'lst*_2007*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    tag <- which(substr(file,91,92)=='06'|substr(file,91,92)=='07'|substr(file,91,92)=='08')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    # temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm1_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm1_LST <- stack(ltm1_lst)
ltm1_LST_sum <- calc(ltm1_LST,mean,na.rm=T)

ltm2_lst <- vector('list',5)
for(yy in 1:5){
  if(yy==1){
    sstr <- 'lst*_2016*'  
  }else if(yy==2){
    sstr <- 'lst*_2017*'  
  }else if(yy==3){
    sstr <- 'lst*_2018*'  
  }else if(yy==4){
    sstr <- 'lst*_2019*'  
  }else{
    sstr <- 'lst*_2020*'  
  }
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    tag <- which(substr(file,91,92)=='06'|substr(file,91,92)=='07'|substr(file,91,92)=='08')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    # temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm2_lst[[yy]] <- do.call(mosaic,rast)
  print(yy)
}
ltm2_LST <- stack(ltm2_lst)
ltm2_LST_sum <- calc(ltm2_LST,mean,na.rm=T)

diff_LST_sum <- ltm2_LST_sum-ltm1_LST_sum





### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/domain.shp')
nest <-c('10')
erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]

rastList <- vector('list',9)
for(i in 1:9){
  if(i==1){
    rast <- ltm1_LST_ann
  }else if(i==2){
    rast <- ltm2_LST_ann
  }else if(i==3){
    rast <- diff_LST_ann
  }else if(i==4){
    rast <- ltm1_LST_spr
  }else if(i==5){
    rast <- ltm2_LST_spr
  }else if(i==6){
    rast <- diff_LST_spr
  }else if(i==7){
    rast <- ltm1_LST_sum
  }else if(i==8){
    rast <- ltm2_LST_sum
  }else{
    rast <- diff_LST_sum
  }
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 5000
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  rastList[[i]] <- rast
  
  print(i)
}



#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(rastList[[1]]))

png(filename='ltm_diff_na_lst.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- rastList[[1]]
temp[temp<260] <- 260
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Ann.avg',cex=2.2,pos=4)
text(-5300000,-2400000,'2003-2007',cex=2,pos=4)
##
temp <- rastList[[4]]
temp[temp<260] <- 260
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Spring',cex=2.2,pos=4)
text(-5300000,-2400000,'2003-2007',cex=2.2,pos=4)
##
temp <- rastList[[7]]
temp[temp<260] <- 260
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Summer',cex=2.2,pos=4)
text(-5300000,-2400000,'2003-2007',cex=2.2,pos=4)


##
temp <- rastList[[2]]
temp[temp<260] <- 260
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Ann.avg',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2,pos=4)
##
temp <- rastList[[5]]
temp[temp<260] <- 260
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Spring',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)
##
temp <- rastList[[8]]
temp[temp<270] <- 270
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(260,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(270,280,290,300,310),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'Summer',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)


##
temp <- rastList[[3]]
temp[temp< -2.5] <- -2.5
temp[temp>  2.5] <- 2.5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-2.5,2.5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-2,-1,0,1,2),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'Ann.avg',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[6]]
temp[temp< -2.5] <- -2.5
temp[temp>  2.5] <- 2.5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-2.5,2.5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-2,-1,0,1,2),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'Spring',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[9]]
temp[temp< -2.5] <- -2.5
temp[temp>  2.5] <- 2.5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-2.5,2.5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-2,-1,0,1,2),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'Summer',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

dev.off()








###########################################################
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- 'met*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

# EVImax long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  temp <- apply(pek[,1:5],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_pek <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # pek[pek[,21]<50,1:20] <- NA
  temp <- apply(pek[,16:20],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_pek <- do.call(mosaic,rast)

diff_pek <- ltm2_pek-ltm1_pek

# EVIarea long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # ear[ear[,21]<50,1:20] <- NA
  temp <- apply(ear[,1:5],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_ear <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # ear[ear[,21]<50,1:20] <- NA
  temp <- apply(ear[,16:20],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_ear <- do.call(mosaic,rast)

diff_ear <- ltm2_ear-ltm1_ear


## Regression
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/stat/'
sstr <- 'reg*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  rt <- setValues(panel_grid,regCoef[,3])
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm1_reg <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  rt <- setValues(panel_grid,regCoef[,4])
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm2_reg <- do.call(mosaic,rast)

rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  rt <- setValues(panel_grid,regCoefphe)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_regp <- do.call(mosaic,rast)


### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/domain.shp')
nest <-c('10')
erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]

rastList <- vector('list',9)
for(i in 1:9){
  if(i==1){
    rast <- ltm1_pek
  }else if(i==2){
    rast <- ltm2_pek
  }else if(i==3){
    rast <- diff_pek
  }else if(i==4){
    rast <- ltm1_ear
  }else if(i==5){
    rast <- ltm2_ear
  }else if(i==6){
    rast <- diff_ear
  }else if(i==7){
    rast <- ltm1_reg
  }else if(i==8){
    rast <- ltm2_reg
  }else{
    rast <- ltm_regp
  }
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 5000
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  rastList[[i]] <- rast
  
  print(i)
}



#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(rastList[[1]]))

png(filename='ltm_reg_na.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- rastList[[1]]
temp[temp<2000] <- 2000
temp[temp>7000] <- 7000
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(9,'Greens'))
plot(rast,
     legend.only=T,
     col=brewer.pal(9,'Greens'),
     zlim=c(2000,7000),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(3000,4000,5000,6000),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EVImax',cex=2.2,pos=4)
text(-5300000,-2400000,'2001-2005',cex=2,pos=4)
##
temp <- rastList[[4]]
temp[temp<30] <- 30
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(9,'Greens'))
plot(rast,
     legend.only=T,
     col=brewer.pal(9,'Greens'),
     zlim=c(30,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(50,150,250),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EVIarea',cex=2.2,pos=4)
text(-5300000,-2400000,'2001-2005',cex=2.2,pos=4)
##
temp <- rastList[[7]]
temp[temp<0.05] <- 0.05
temp[temp>0.55] <- 0.55
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(0.05,0.55),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(0.1,0.3,0.5),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'R2 Spr.',cex=2.2,pos=4)
text(-5300000,-2400000,'SOS vs. LST',cex=2.2,pos=4)


##
temp <- rastList[[2]]
temp[temp<2000] <- 2000
temp[temp>7000] <- 7000
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(9,'Greens'))
plot(rast,
     legend.only=T,
     col=brewer.pal(9,'Greens'),
     zlim=c(2000,7000),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(3000,4000,5000,6000),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EVImax',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2,pos=4)
##
temp <- rastList[[5]]
temp[temp<30] <- 30
temp[temp>320] <- 320
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(9,'Greens'))
plot(rast,
     legend.only=T,
     col=brewer.pal(9,'Greens'),
     zlim=c(30,320),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(50,150,250),
                    cex.axis=2,font=1))   
text(-5300000,-1700000,'EVIarea',cex=2.2,pos=4)
text(-5300000,-2400000,'2016-2020',cex=2.2,pos=4)
##
temp <- rastList[[8]]
temp[temp<0.01] <- 0.01
temp[temp>0.5] <- 0.5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'Spectral'))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(0.01,0.5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(0.5,0.25,0.01),c(0.01,0.25,0.5),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'p-value Spr.',cex=2.2,pos=4)
text(-5300000,-2400000,'SOS vs. LST',cex=2.2,pos=4)


##
temp <- rastList[[3]]
temp[temp< -800] <- -800
temp[temp>  800] <- 800
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-800,800),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-800,-400,0,400,800),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'EVImax',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[6]]
temp[temp< -50] <- -50
temp[temp>  50] <- 50
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-50,50),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-50,-25,0,25,50),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'EVIarea',cex=2.2,pos=4)
text(-5300000,-2400000,'Change',cex=2.2,pos=4)

temp <- rastList[[9]]
temp[temp< -0.6] <- -0.6
temp[temp>  0.6] <- 0.6
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(-0.6,0.6),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-0.5,0,0.5),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'r',cex=2.2,pos=4)
text(-5300000,-2400000,'SOS vs. EOS',cex=2.2,pos=4)

dev.off()











# EVI and LST Anomaly
###########################################################
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- 'metrics_h*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

# EVImax and EVIarea long-term mean
rast <- vector('list',length(file))
rast1 <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  temp <- apply(pek[,1:20],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
  
  temp <- apply(ear[,1:20],1,median,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast1[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_pek_all <- do.call(mosaic,rast)

rast1$fun <- mean
rast1$na.rm <- T
ltm_ear_all <- do.call(mosaic,rast1)

# EVImax and EVIarea 2012
rast <- vector('list',length(file))
rast1 <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  rt <- setValues(panel_grid,pek[,12])
  rast[[i]] <- rt
  
  rt <- setValues(panel_grid,ear[,12])
  rast1[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_pek_2012 <- do.call(mosaic,rast)

rast1$fun <- mean
rast1$na.rm <- T
ltm_ear_2012 <- do.call(mosaic,rast1)

ano_pek <- ltm_pek_2012 - ltm_pek_all
ano_ear <- ltm_ear_2012 - ltm_ear_all

####
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/met'

# Long-term mean Ann
ltm_lst_all <- vector('list',18)
ltm_lst_spr <- vector('list',18)
ltm_lst_sum <- vector('list',18)
for(yy in 1:18){
  sstr <- paste0('lst*_',(2002+yy),'*')  
  
  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  rast <- vector('list',length(files))
  rast1 <- vector('list',length(files))
  rast2 <- vector('list',length(files))
  for(i in 1:length(files)){
    load(files[i])
    # tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    # temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    temp <- apply(lstMat,1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast[[i]] <- rt
    
    tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast1[[i]] <- rt
    
    tag <- which(substr(file,91,92)=='06'|substr(file,91,92)=='07'|substr(file,91,92)=='08')
    temp <- apply(lstMat[,tag],1,mean,na.rm=T)
    rt <- setValues(panel_grid,temp)
    rast2[[i]] <- rt
  }
  rast$fun <- mean
  rast$na.rm <- T
  ltm_lst_all[[yy]] <- do.call(mosaic,rast)
  
  rast1$fun <- mean
  rast1$na.rm <- T
  ltm_lst_spr[[yy]] <- do.call(mosaic,rast1)
  
  rast2$fun <- mean
  rast2$na.rm <- T
  ltm_lst_sum[[yy]] <- do.call(mosaic,rast2)
  
  print(yy)
}
ltm_LST <- stack(ltm_lst_all)
ltm_LST_ann <- calc(ltm_LST,mean,na.rm=T)

ltm_LST <- stack(ltm_lst_spr)
ltm_LST_spr <- calc(ltm_LST,mean,na.rm=T)

ltm_LST <- stack(ltm_lst_sum)
ltm_LST_sum <- calc(ltm_LST,mean,na.rm=T)

sstr <- 'lst*_2012*'  
files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
rast <- vector('list',length(files))
rast1 <- vector('list',length(files))
rast2 <- vector('list',length(files))
for(i in 1:length(files)){
  load(files[i])
  # tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
  # temp <- apply(lstMat[,tag],1,mean,na.rm=T)
  temp <- apply(lstMat,1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
  
  tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
  temp <- apply(lstMat[,tag],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast1[[i]] <- rt
  
  tag <- which(substr(file,91,92)=='06'|substr(file,91,92)=='07'|substr(file,91,92)=='08')
  temp <- apply(lstMat[,tag],1,mean,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast2[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_LST_2012 <- do.call(mosaic,rast)

rast1$fun <- mean
rast1$na.rm <- T
ltm_LST_2012_spr <- do.call(mosaic,rast1)

rast2$fun <- mean
rast2$na.rm <- T
ltm_LST_2012_sum <- do.call(mosaic,rast2)


ano_LST_2012     <- ltm_LST_2012 - ltm_LST_ann
ano_LST_2012_spr <- ltm_LST_2012_spr - ltm_LST_spr
ano_LST_2012_sum <- ltm_LST_2012_sum - ltm_LST_sum


### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/domain.shp')
nest <-c('10')
erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]

rastList <- vector('list',5)
for(i in 1:5){
  if(i==1){
    rast <- ano_pek
  }else if(i==2){
    rast <- ano_ear
  }else if(i==3){
    rast <- ano_LST_2012
  }else if(i==4){
    rast <- ano_LST_2012_spr
  }else{
    rast <- ano_LST_2012_sum
  }
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 5000
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  rastList[[i]] <- rast
  
  print(i)
}


setwd('/projectnb/modislc/users/mkmoon/TAscience/maps/')
writeRaster(rastList[[2]],filename="ano_2012_ear.tif", format="GTiff", overwrite=TRUE)
writeRaster(rastList[[4]],filename="ano_2012_lst_spr.tif", format="GTiff", overwrite=TRUE)
writeRaster(rastList[[5]],filename="ano_2012_lst_sum.tif", format="GTiff", overwrite=TRUE)



#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(rastList[[1]]))

png(filename='ano_evi_lst.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- rastList[[3]]
temp[temp< -5] <- -5
temp[temp>  5] <-  5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-5,5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-4,-2,0,2,4),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'LST Ann.',cex=2,pos=4)
text(-5300000,-2400000,'2012 Anomaly',cex=1.8,pos=4)

##
temp <- rastList[[4]]
temp[temp< -5] <- -5
temp[temp>  5] <-  5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-5,5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-4,-2,0,2,4),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'LST MAM',cex=2,pos=4)
text(-5300000,-2400000,'2012 Anomaly',cex=1.8,pos=4)
##
temp <- rastList[[5]]
temp[temp< -5] <- -5
temp[temp>  5] <-  5
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'PiYG')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'PiYG')),
     zlim=c(-5,5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-4,-2,0,2,4),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'LST JJA',cex=2,pos=4)
text(-5300000,-2400000,'2012 Anomaly',cex=1.8,pos=4)

##
temp <- rastList[[1]]
temp[temp< -500] <- -500
temp[temp>  500] <- 500
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-0.05,0.05),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-0.04,-0.02,0,0.02,0.04),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EVImax',cex=2,pos=4)
text(-5300000,-2400000,'2012 Anomaly',cex=1.8,pos=4)

##
temp <- rastList[[2]]
temp[temp< -50] <- -50
temp[temp>  50] <- 50
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=brewer.pal(11,'PiYG'))
plot(rast,
     legend.only=T,
     col=brewer.pal(11,'PiYG'),
     zlim=c(-5,5),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(-4,-2,0,2,4),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'EVIarea',cex=2,pos=4)
text(-5300000,-2400000,'2012 Anomaly',cex=1.8,pos=4)

dev.off()









###########################################################
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- 'metrics_h*'
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

# SOS SD
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # sos[sos[,21]<50,1:20] <- NA
  temp <- apply(sos[,1:20],1,sd,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_sd_sos <- do.call(mosaic,rast)


# EOS SD
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # eos[eos[,21]<50,1:20] <- NA
  temp <- apply(eos[,1:20],1,sd,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_sd_eos <- do.call(mosaic,rast)


# GSL long-term mean
rast <- vector('list',length(file))
for(i in 1:length(file)){
  load(file[i])
  # gsl1[gsl1[,21]<50,1:20] <- NA
  temp <- apply(gsl1[,1:20],1,sd,na.rm=T)
  rt <- setValues(panel_grid,temp)
  rast[[i]] <- rt
}
rast$fun <- mean
rast$na.rm <- T
ltm_sd_gsl <- do.call(mosaic,rast)


### ecoregion by EPA)
erepa <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/domain.shp')
nest <-c('10')
erepa <- erepa[as.character(erepa$NA_L1CODE) %in% nest,]

rastList <- vector('list',3)
for(i in 1:3){
  if(i==1){
    rast <- ltm_sd_sos
  }else if(i==2){
    rast <- ltm_sd_eos
  }else if(i==3){
    rast <- ltm_sd_gsl
  }
  
  pr3 <- projectExtent(erepa,crs(erepa))
  res(pr3) <- 5000
  rast <- projectRaster(rast,pr3)
  rast <- mask(rast,erepa)  
  
  rastList[[i]] <- rast
  
  print(i)
}

setwd('/projectnb/modislc/users/mkmoon/TAscience/maps/')
writeRaster(rastList[[1]],filename="sd_sos.tif", format="GTiff", overwrite=TRUE)
writeRaster(rastList[[2]],filename="sd_eos.tif", format="GTiff", overwrite=TRUE)
writeRaster(rastList[[3]],filename="sd_gsl.tif", format="GTiff", overwrite=TRUE)


#############
setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/')
bound <- shapefile('/projectnb/modislc/users/mkmoon/TAscience/shp/na.shp')
bound <- spTransform(bound,crs(rastList[[1]]))

png(filename='ltm_sd.png',width=12,height=8,unit='in',res=300)

par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(0,0,0.5,0))
##
temp <- rastList[[1]]
temp[temp< 4] <- 4
temp[temp>  16] <-  16
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(4,16),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(5,10,15),
                    cex.axis=2,font=1))    
text(-5300000,-1700000,'SOS',cex=2,pos=4)
text(-5300000,-2400000,'Std Dev (days)',cex=1.8,pos=4)

##
temp <- rastList[[2]]
temp[temp< 4] <- 4
temp[temp>  16] <-  16
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(4,16),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(5,10,15),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'EOS',cex=2,pos=4)
text(-5300000,-2400000,'Std Dev (days)',cex=1.8,pos=4)
##
temp <- rastList[[3]]
temp[temp< 4] <- 4
temp[temp>  16] <-  16
plot(temp,axes=F,box=F,cex.main=1.2,legend=F)
plot(bound,add=T,col='grey95',border='grey95')
plot(erepa,add=T,col='grey75',border='grey75')
plot(temp,axes=F,box=F,cex.main=1.2,legend=F,add=T,
     col=rev(brewer.pal(11,'Spectral')))
plot(rast,
     legend.only=T,
     col=rev(brewer.pal(11,'Spectral')),
     zlim=c(4,16),
     legend.width=1.5,
     legend.shrink=0.6,
     horiz=F,
     smallplot=c(0.06,0.09,0.25,0.68),
     axis.args=list(at=c(5,10,15),
                    cex.axis=2,font=1)) 
text(-5300000,-1700000,'GSL',cex=2,pos=4)
text(-5300000,-2400000,'Std Dev (days)',cex=1.8,pos=4)

dev.off()
