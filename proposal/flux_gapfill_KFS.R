library(REddyProc)
# +++ Input data from csv (example needs to be downloaded)
daf <- read.table('/projectnb/modislc/users/mkmoon/autumn/sample_data/Example_DETha98.txt',header=T)

ss=9
sites <- c('US-WCr','US-UMd','US-UMB','CA-TPD','US-Bar','US-MMS','US-Ha2','US-xHa','US-KFS')
lat <- c('-','-','-','-',44.0646,39.3232,42.5393,42.5369,39.0561)
lon <- c('-','-','-','-',-71.28808,-86.4131,-72.1779,-72.17266,-95.1907)

# Read data
Dat <- read.csv('/projectnb/modislc/users/mkmoon/Download/AMF_US-KFS_BASE_HH_7-5.csv',skip=2)
year <- as.numeric(substr(Dat[,1],1,4))

GPP <- matrix(NA,365,12)
for(YY in 1:12){
  options(warn=-1)
  
dat <- Dat[which(year==(2007+YY)),]
yy <- as.numeric(substr(dat[,1],1,4))
mm <- as.numeric(substr(dat[,1],5,6))
dd <- as.numeric(substr(dat[,1],7,8))
doy <- as.numeric(strftime(paste(yy,'-',mm,'-',dd,sep=''),format='%j'))

dates <- as.Date(paste(yy,'-',mm,'-',dd,sep=''))
hh <- as.numeric(substr(dat[,1],9,10))
mi <- as.numeric(substr(dat[,1],11,12))
mi[mi==30] <- 0.5
hm <- as.numeric(hh)+as.numeric(mi)

# Rearrange input data
ndat <- data.frame(yy,doy,hm,dat$FC_1_1_1,dat$LE_1_1_1,dat$H_1_1_1,dat$SW_IN_1_1_1,dat$TA_1_1_1,
                   dat$TS_2_1_1,dat$RH_1_1_1,dat$VPD_PI_1_1_1,dat$USTAR_1_1_1)

colnames(ndat) <- colnames(daf)
ndat[which(ndat[,7]<0),7] <- 0

ndat[ndat==-9999] <- NA
EddyData <- ndat

#+++ If not provided, calculate VPD from Tair and rH
EddyData$VPD <- fCalcVPDfromRHandTair(EddyData$rH, EddyData$Tair)
#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- EddyData %>% 
  filterLongRuns("NEE") %>% 
  fConvertTimeToPosix('YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')
#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EProc <- sEddyProc$new(sites[ss], EddyDataWithPosix, c('NEE','Rg','Tair','VPD', 'Ustar'))
#Location of DE-Tharandt
EProc$sSetLocationInfo(LatDeg = as.numeric(lat[ss]),
                       LongDeg = as.numeric(lon[ss]), TimeZoneHour = -5)  
#
#++ Fill NEE gaps with MDS gap filling algorithm (without prior ustar filtering)
EProc$sMDSGapFill('NEE', FillAll=F)
EProc$sMDSGapFill('Tair', FillAll=F)
EProc$sMRFluxPartition(FluxVar='NEE_f',
                       QFFluxVar='NEE_fqc',
                       QFFluxValue=0L,
                       TempVar='Tair_f',
                       QFTempVar='Tair_fqc',
                       QFTempValue=0,
                       RadVar='Rg',
                       TRef=273.15+15)
#++ Export gap filled and partitioned data to standard data frame
FilledEddyData <- EProc$sExportResults()

#++ Example plots of filled data to screen or to directory \plots
# EProc$sPlotFingerprintY('Reco', Year = 2008)
# plot(EddyDataWithPosix[,1],FilledEddyData$GPP_f)
# 
# gpp <- data.frame(EddyDataWithPosix[,1],FilledEddyData$GPP_f)
# 
# nee <- FilledEddyData$NEE_f
# res <- FilledEddyData$Reco
gpp <- FilledEddyData$GPP_f

# aa <- -nee+res
# plot(aa,gpp)

gpp[gpp<0] <- 0
dgpp <- matrix(gpp[1:(48*365)],48,365)
dgpp <- apply(dgpp,2,sum)
GPP[,YY] <- dgpp
print(YY)
}

par(mfrow=c(2,2))
for(i in 1:12){
  plot(GPP[,i])
}
meGPP <- apply(GPP,1,mean)
plot(meGPP)
points(GPP[,5],col='red')





##########################################
# Reflectance
red <- read.csv('/projectnb/modislc/users/mkmoon/TAscience/Data/kfs/statistics_sur_refl_b01_a.csv')
nir <- read.csv('/projectnb/modislc/users/mkmoon/TAscience/Data/kfs/statistics_sur_refl_b02_a.csv')
RED <- red[which(substr(red$dt,1,4)>2007&substr(red$dt,1,4)<2020),]
NIR <- nir[which(substr(nir$dt,1,4)>2007&substr(nir$dt,1,4)<2020),]

EVI <- matrix(NA,365,12)
for(i in 1:12){
  red <- RED[which(substr(RED$dt,1,4)==(2007+i)),]
  nir <- NIR[which(substr(NIR$dt,1,4)==(2007+i)),]
  
  dd  <- as.numeric(strftime(red$dt, format = "%j"))
  evi <- 2.5*(nir$value_mean-red$value_mean)/(nir$value_mean+2.4*red$value_mean+1.0)

  spl <- smooth.spline(dd,evi,spar=0.55)
  smoothed <- predict(spl,1:365)
    
  EVI[,i] <- smoothed$y
}
meEVI <- apply(EVI,1,mean)
plot(meEVI)
points(EVI[,5],col='red')

# LST
lstd <- read.csv('/projectnb/modislc/users/mkmoon/TAscience/Data/kfs/statistics_LST_Day_1KM.csv')
lstn <- read.csv('/projectnb/modislc/users/mkmoon/TAscience/Data/kfs/statistics_LST_Night_1KM.csv')
LSTd <- lstd[which(substr(lstd$dt,1,4)>2007&substr(lstd$dt,1,4)<2020),]
LSTn <- lstn[which(substr(lstn$dt,1,4)>2007&substr(lstn$dt,1,4)<2020),]

LST <- matrix(NA,365,12)
for(i in 1:12){
  lstd <- LSTd[which(substr(LSTd$dt,1,4)==(2007+i)),]
  lstn <- LSTn[which(substr(LSTn$dt,1,4)==(2007+i)),]
  
  dd  <- as.numeric(strftime(lstd$dt, format = "%j"))
  # lst <- (lstd$value_mean+lstn$value_mean)/2
  lst <- lstd$value_mean
  
  lst[lst<150] <-NA
  dd[is.na(lst)] <-NA
  
  lst <- na.omit(lst)
  dd <- na.omit(dd)
  
  spl <- smooth.spline(dd,lst,spar=0.55)
  smoothed <- predict(spl,1:365)
  
  LST[,i] <- smoothed$y
}
meLST <- apply(LST,1,mean)
plot(meLST)
points(LST[,5],col='red')



#######################
# png(filename='ano_2012_us-mms.png',width=6,height=6,unit='in',res=300)
# par(mfrow=c(3,1),oma=c(2,1,0,1),mar=c(0.2,5,1,0),mgp=c(2.8,1,0))

png(filename='ano_2012_us-kfs.png',width=11.9,height=3.6,unit='in',res=300)
par(mfrow=c(1,3),oma=c(2,1,0,0),mar=c(4,5,1,1),mgp=c(2.8,1,0))


plot(meEVI,type='l',lwd=3,axe=F,ann=F,ylim=c(0.1,0.7))
axis(2,seq(0,1,0.2),cex.axis=2)
axis(1,at=seq(0,365,100),cex.axis=2)
mtext('EVI2',2,2.7,cex=2)
box()
lines(EVI[,5],col='red',lwd=5)
legend('topleft',c('Mean (2008-2019)','2012'),lwd=c(3,5),bty='n',cex=2.2,col=c('black','red'))


spl <- smooth.spline(1:364,meGPP[1:364],spar=0.55)
smoothed <- predict(spl,1:365)

plot(smoothed$y*0.020833333,type='l',lwd=3,axe=F,ann=F,ylim=c(0,9))
axis(2,seq(0,10,2),cex.axis=2)
axis(1,at=seq(0,365,100),cex.axis=2)
mtext('GPP',2,2.7,cex=2)
box()

spl <- smooth.spline(1:364,GPP[1:364,5],spar=0.55)
smoothed <- predict(spl,1:365)
lines(smoothed$y*0.020833333,col='red',lwd=5)
# legend('topright',expression(paste(mu,'mol',CO[2],' ',m^-2,s^-1)),bty='n',cex=2.2)


plot(meLST,type='l',lwd=3,axe=F,ann=F,ylim=c(270,330))
axis(2,seq(200,400,20),cex.axis=2)
axis(1,at=seq(0,365,100),cex.axis=2)
mtext('LST',2,2.7,cex=2)
box()
lines(LST[,5],col='red',lwd=5)

mtext('Day of year',1,2.7,cex=1.7,line=-0.5)

dev.off()