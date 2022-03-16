rm(list = ls())

library(gdalUtils)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(RColorBrewer)
library(viridis)
library(doMC)
library(lmodel2)


##############################
args <- commandArgs()
print(args)

tile <- substr(args[3],1,6)
# tile <- 'h12v04'

path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/metrics'
sstr <- paste0('metrics_',tile,'*')
file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)

load(file)


##
path <- '/projectnb/modislc/users/mkmoon/TAscience/Data/met'

lstAnn <- matrix(NA,(240*240),18)
lstSpr <- matrix(NA,(240*240),18)
lstSum <- matrix(NA,(240*240),18)
for(i in 1:18){
  sstr <- paste0('lst_',tile,'_',(i+2002),'*')

  files <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
  load(files)

  lstAnn[,i] <- apply(lstMat,1,mean,na.rm=T)

  tag <- which(substr(file,91,92)=='03'|substr(file,91,92)=='04'|substr(file,91,92)=='05')
  lstSpr[,i] <- apply(lstMat[,tag],1,mean,na.rm=T)

  tag <- which(substr(file,91,92)=='06'|substr(file,91,92)=='07'|substr(file,91,92)=='08')
  lstSum[,i] <- apply(lstMat[,tag],1,mean,na.rm=T)

  print(i)
}


# ##############################
# setwd('/projectnb/modislc/users/mkmoon/TAscience/figure/ts')
# 
# for(i in 1:30){
#   pp <- sample((240*240),1)
#   
#   if(sum(!is.na(sos[pp,1:20]))==20){
#     PP <- sprintf('%05d',pp)
#     png(filename=paste0('ts_',tile,'_',PP,'.png'),width=12,height=8,unit='in',res=300)
#     
#     par(mfrow=c(3,1),oma=c(2,1,1,1),mar=c(4,5,1,5),mgp=c(2.5,1,0))
#     
#     plot(2000:2020,c(NA,sos[pp,1:20]),pch=19,col='forestgreen',axe=F,ann=F,type='o',cex=1.5)
#     axis(2,seq(0,400,10),cex.axis=1.5)
#     mtext('SOS',2,2.7,cex=1.5)
#     box(lty=1)
#     axis(1,at=seq(2000,2020,5),seq(2000,2020,5),cex.axis=1.5)
#     
#     
#     plot(2000:2020,c(NA,NA,NA,lstAnn[pp,]),pch=19,col='red',axe=F,ann=F,type='o',cex=1.5,
#          ylim=c(min(lstSpr[pp,],na.rm=T),max(lstSum[pp,],na.rm=T)))
#     legend('topleft',c('Ann.avg','Spring','Summer'),pch=19,lty=1,bty='n',cex=2,
#            col=c('red','blue','brown'))
#     axis(2,seq(0,500,5),cex.axis=1.5)
#     mtext('LST',2,2.7,cex=1.5)
#     box(lty=1)
#     axis(1,at=seq(2000,2020,5),seq(2000,2020,5),cex.axis=1.5)
#     points(2000:2020,c(NA,NA,NA,lstSpr[pp,]),pch=19,col='blue',type='o',cex=1.5)
#     points(2000:2020,c(NA,NA,NA,lstSum[pp,]),pch=19,col='brown',type='o',cex=1.5)
#     
#     par(fig=c(0,0.5,0,0.33),oma=c(2,1,0,1),mar=c(4,5,0,5),mgp=c(2.5,1,0),new=T)
#     plot(c(NA,NA,lstAnn[pp,]),sos[pp,1:20],pch=19,cex=1.5,col='red',
#          xlab='Ann.avg LST',ylab='SOS',cex.lab=2,cex.axis=1.5)
#     abline(lm(sos[pp,1:20]~c(NA,NA,lstAnn[pp,])))
#     reg <- summary(lm(sos[pp,1:20]~c(NA,NA,lstAnn[pp,])))
#     legend('topright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
#     par(fig=c(0.5,1,0,0.33),oma=c(2,1,0,1),mar=c(4,5,0,5),mgp=c(2.5,1,0),new=T)
#     plot(c(NA,NA,lstSpr[pp,]),sos[pp,1:20],pch=19,cex=1.5,col='blue',
#          xlab='Spring LST',ylab='SOS',cex.lab=2,cex.axis=1.5)
#     abline(lm(sos[pp,1:20]~c(NA,NA,lstSpr[pp,])))
#     reg <- summary(lm(sos[pp,1:20]~c(NA,NA,lstSpr[pp,])))
#     legend('topright',as.character(round(reg$r.squared,3)),bty='n',cex=3)
#     
#     dev.off() 
#   }
# }


################
mcd12q1_path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01'
sstr <- paste0('*',tile,'*')
file <- list.files(path=mcd12q1_path,pattern=glob2rx(sstr),full.names=T)
sds <- get_subdatasets(file)
lct <- raster(sds[1])

sinu_crs <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

## panel grid
cxx <- 10
cyy <- 10

panel_grid <- matrix(1:length(lct),dim(lct)[1]/cxx,dim(lct)[2]/cyy)
panel_grid <- raster(panel_grid,crs=sinu_crs,
                     xmn=lct@extent@xmin,xmx=lct@extent@xmax,
                     ymn=lct@extent@ymin,ymx=lct@extent@ymax)


###
registerDoMC()

################
regCoef <- foreach(i=1:(240*240),.combine='rbind') %dopar% {
  
  if(sum(!is.na(sos[i,3:20]))>14 & sum(!is.na(lstAnn[i,]))>14 & sum(!is.na(lstSpr[i,]))>14 & sum(!is.na(lstSum[i,]))>14){
    # reg <- summary(lm(sos[i,3:20]~lstAnn[i,]))
    # r1 <- reg$r.squared;  p1 <- reg$coefficients[8]  
    # 
    # reg <- summary(lm(sos[i,3:20]~lstSpr[i,]))
    # r2 <- reg$r.squared;  p2 <- reg$coefficients[8]  
    # 
    # reg <- summary(lm(sos[i,3:20]~lstSum[i,]))
    # r3 <- reg$r.squared;  p3 <- reg$coefficients[8]  
    
    reg <- (lmodel2(sos[i,3:20]~lstAnn[i,]))
    r1 <- reg$r;  s1 <- reg$regression.results$Slope[3]; p1 <- reg$P.param

    reg <- (lmodel2(sos[i,3:20]~lstSpr[i,]))
    r2 <- reg$r;  s2 <- reg$regression.results$Slope[3]; p2 <- reg$P.param

    reg <- (lmodel2(sos[i,3:20]~lstSum[i,]))
    r3 <- reg$r;  s3 <- reg$regression.results$Slope[3]; p3 <- reg$P.param
    
    
    c(r1,s1,p1,r2,s2,p2,r3,s3,p3)  
  }else{
    rep(NA,9)
  }
  
}
print(dim(regCoef))

colnames(regCoef) <- rep(c('r','slope','pval'),3)


################
regCoefphe <- foreach(i=1:(240*240),.combine='rbind') %dopar% {
  
  if(sum(!is.na(sos[i,1:20]))>14 & sum(!is.na(eos[i,1:20]))>14){
    temp <- cbind(eos[i,1:20],sos[i,1:20])
    # temp <- na.omit(temp)
    # cor(temp)[2] 
    reg <- lmodel2(temp[,1]~temp[,2])
    c(reg$r,reg$regression.results$Slope[3],reg$P.param)
  }else{
    rep(NA,3)
  }
}
print(dim(regCoefphe))

colnames(regCoefphe) <- c('r','slope','pval')


#########
if(dim(regCoef)[1]==(240*240)){
  setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/reg/')
  save(regCoef,regCoefphe,panel_grid,
       file=paste0('1_reg_',tile,'.rda'))
}
