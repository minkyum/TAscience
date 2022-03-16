rm(list = ls())

library(gdalUtils)
library(raster)
library(rgdal)
library(RColorBrewer)
library(zyp)
library(Kendall)
library(scales)
library(agricolae)
library(profvis)

############################################################
# From bash code
args <- commandArgs()
print(args)

tt <- as.numeric(args[3])
# tt <- 15

sinu_crs <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

## Get MCD12Q1
mat_lct <- matrix(NA,(2400*2400),19)
for(i in 1:19){
  mcd12q1_path <- paste('/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/',(i+2000),'.01.01',sep='')
  file <- list.files(path=mcd12q1_path,full.names=T)
  tile_list <- substr(file,106,111)
  sds <- get_subdatasets(file[tt])
  lct <- raster(sds[1])
  mat_lct[,i] <- values(lct)
  
  print(i)
}

lct_ch <- matrix(NA,(2400*2400),1)
for(i in 1:(2400*2400)){
  temp <- mat_lct[i,]
  if(length(unique(temp))==1){
    lct_ch[i] <- 1
  }else{
    lct_ch[i] <- 0
  }
}
lct_ch <- setValues(lct,lct_ch)
rm(mat_lct)

mcd12q1_path <- '/projectnb/modislc/users/mkmoon/mcd12q1/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2019.01.01'
file <- list.files(path=mcd12q1_path,full.names=T)
tile_list <- substr(file,106,111)
sds <- get_subdatasets(file[tt])
lct <- raster(sds[1])
lct_mat <- as.matrix(lct)

## panel grid
cxx <- 10
cyy <- 10

panel_grid <- matrix(1:length(lct),dim(lct)[1]/cxx,dim(lct)[2]/cyy)
panel_grid <- raster(panel_grid,crs=sinu_crs,
                     xmn=lct@extent@xmin,xmx=lct@extent@xmax,
                     ymn=lct@extent@ymin,ymx=lct@extent@ymax)

## portion of deciduous pixels
pofdp <- matrix(1,dim(panel_grid)[1],dim(panel_grid)[2])
poflc <- matrix(1,dim(panel_grid)[1],dim(panel_grid)[2])
for(i in 1:dim(panel_grid)[1]){
  for(j in 1:dim(panel_grid)[2]){
    temp <- lct_mat[(cxx*(i-1)+1):(cxx*i),(cyy*(j-1)+1):(cyy*j)]
    temp[temp==0] <- NA
    temp <- na.omit(temp)
    
    pofdp[i,j] <- sum(temp==4|temp==5)/length(temp)*100
    
    aa <- as.matrix(table(temp))
    poflc[i,j] <- as.numeric(row.names(aa)[which(aa==max(aa))])[1]
  }
  if(i%%10==0) print(i)
}
pdp_raster <- raster(pofdp,crs=sinu_crs,
                     xmn=lct@extent@xmin,xmx=lct@extent@xmax,
                     ymn=lct@extent@ymin,ymx=lct@extent@ymax)
plc_raster <- raster(poflc,crs=sinu_crs,
                     xmn=lct@extent@xmin,xmx=lct@extent@xmax,
                     ymn=lct@extent@ymin,ymx=lct@extent@ymax)


### load SOS
sos <- matrix(NA,length(panel_grid),22)
eos <- matrix(NA,length(panel_grid),22)
gsl1 <- matrix(NA,length(panel_grid),22)
gsl2 <- matrix(NA,length(panel_grid),22)
rup <- matrix(NA,length(panel_grid),22)
rdn <- matrix(NA,length(panel_grid),22)
rsn <- matrix(NA,length(panel_grid),22)
pek <- matrix(NA,length(panel_grid),22)
ear <- matrix(NA,length(panel_grid),22)
for(ii in 1:20){
  yy <- 2000 + ii
  if(ii<19){
    path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/product/',yy,sep='')
    sstr <- paste('*',yy,'*',tile_list[tt],'*.hdf',sep='')
    file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
    temp <- get_subdatasets(file)
    rast1 <- raster(temp[3])  
    rast2 <- raster(temp[7])  
    rast3 <- raster(temp[2])  
    rast4 <- raster(temp[5])  
    rast5 <- raster(temp[6])  
    rast6 <- raster(temp[8])  
    rast71 <- raster(temp[9])
    rast72 <- raster(temp[10])      
    rast8 <- raster(temp[11])
  }else{
    path <- paste('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/pheno_out/',tile_list[tt],sep='')
    sstr <- paste('*',tile_list[tt],'_',yy,sep='')
    file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
    rast1 <- raster(file,band=6)  
    rast2 <- raster(file,band=10)  
    rast3 <- raster(file,band=5)  
    rast4 <- raster(file,band=7)  
    rast5 <- raster(file,band=9)  
    rast6 <- raster(file,band=11)  
    rast71 <- raster(file,band=3)
    rast72 <- raster(file,band=4)    
    rast8 <- raster(file,band=2)
  }
  
  # exclude non-forest pixels
  doy_offset <- as.integer(as.Date(paste(yy, "-1-1", sep="")) - as.Date("1970-1-1"))
  
  values(rast1)[values(rast1)>32000] <- NA
  rast1 <- rast1 - doy_offset
  # values(rast1)[values(rast1)<0] <- NA
  # values(rast1)[values(rast1)>365] <- NA

  values(rast2)[values(rast2)>32000] <- NA
  rast2 <- rast2 - doy_offset
  # values(rast2)[values(rast2)<0] <- NA
  # values(rast2)[values(rast2)>365] <- NA
  
  values(rast3)[values(rast3)>32000] <- NA
  rast3 <- rast3 - doy_offset
  # values(rast3)[values(rast3)<0] <- NA
  # values(rast3)[values(rast3)>365] <- NA
  
  values(rast4)[values(rast4)>32000] <- NA
  rast4 <- rast4 - doy_offset
  # values(rast4)[values(rast4)<0] <- NA
  # values(rast4)[values(rast4)>365] <- NA
  
  values(rast5)[values(rast5)>32000] <- NA
  rast5 <- rast5 - doy_offset
  # values(rast5)[values(rast5)<0] <- NA
  # values(rast5)[values(rast5)>365] <- NA
  
  values(rast6)[values(rast6)>32000] <- NA
  rast6 <- rast6 - doy_offset
  # values(rast6)[values(rast6)<0] <- NA
  # values(rast6)[values(rast6)>365] <- NA
  
  values(rast71)[values(rast71)>32000] <- NA
  values(rast72)[values(rast72)>32000] <- NA
  rast7 <- rast71 + rast72
  # values(rast7)[values(rast7)<0] <- NA
  # values(rast7)[values(rast7)>365] <- NA
  
  values(rast8)[values(rast8)>32760] <- NA
  # values(rast8)[values(rast8)<0] <- NA
  # values(rast8)[values(rast8)>365] <- NA
  
  # aggregate as panel
  for(jj in 1:9){
    if(jj==1){
      dd <- as.matrix(rast1) # mid-gup
    }else if(jj==2){
      dd <- as.matrix(rast2) # mid-gdn
    }else if(jj==3){
      dd <- as.matrix(rast2-rast1) # mid-gsl
    }else if(jj==4){
      dd <- as.matrix(rast6-rast3) # gsl
    }else if(jj==5){
      dd <- as.matrix(rast4-rast3) # rate gup
    }else if(jj==6){
      dd <- as.matrix(rast5-rast4) # rate gdn
    }else if(jj==7){
      dd <- as.matrix(rast6-rast5) # rate sen
    }else if(jj==8){
      dd <- as.matrix(rast7) # peak
    }else{
      dd <- as.matrix(rast8) # evi area
    }
    
    dat <- matrix(NA,dim(panel_grid)[1],dim(panel_grid)[2])
    ll <- as.matrix(lct)
    lc <- as.matrix(lct_ch)
    for(i in 1:dim(panel_grid)[1]){
      for(j in 1:dim(panel_grid)[2]){
        temp1 <- dd[(cxx*(i-1)+1):(cxx*i),(cyy*(j-1)+1):(cyy*j)] 
        temp2 <- ll[(cxx*(i-1)+1):(cxx*i),(cyy*(j-1)+1):(cyy*j)] 
        temp3 <- lc[(cxx*(i-1)+1):(cxx*i),(cyy*(j-1)+1):(cyy*j)] 
        # temp1[temp3==0] <- NA 
        # temp <- temp1[temp2==1|temp2==2|temp2==3|temp2==4|temp2==5|temp2==6|temp2==7|temp2==8|temp2==9|temp2==10|temp2==11]
        temp <- temp1
        if(sum(!is.na(temp))>9){
          dat[i,j] <- median(temp,na.rm=T)        
        }else{
          dat[i,j] <- NA
        }
      }
    }
    rast <- raster(dat,crs=sinu_crs,
                   xmn=panel_grid@extent@xmin,xmx=panel_grid@extent@xmax,
                   ymn=panel_grid@extent@ymin,ymx=panel_grid@extent@ymax)
    if(jj==1){
      sos[,ii] <- values(rast)    
    }else if(jj==2){
      eos[,ii] <- values(rast)    
    }else if(jj==3){
      gsl1[,ii] <- values(rast)    
    }else if(jj==4){
      gsl2[,ii] <- values(rast)    
    }else if(jj==5){
      rup[,ii] <- values(rast)    
    }else if(jj==6){
      rdn[,ii] <- values(rast)    
    }else if(jj==7){
      rsn[,ii] <- values(rast)    
    }else if(jj==8){
      pek[,ii] <- values(rast)    
    }else{
      ear[,ii] <- values(rast)    
    }
    print(paste(ii,';',jj))  
  }
}

sos[,21] <- values(pdp_raster)
eos[,21] <- values(pdp_raster)
gsl1[,21] <- values(pdp_raster)
gsl2[,21] <- values(pdp_raster)
rup[,21] <- values(pdp_raster)
rdn[,21] <- values(pdp_raster)
rsn[,21] <- values(pdp_raster)
pek[,21] <- values(pdp_raster)
ear[,21] <- values(pdp_raster)

sos[,22] <- values(plc_raster)
eos[,22] <- values(plc_raster)
gsl1[,22] <- values(plc_raster)
gsl2[,22] <- values(plc_raster)
rup[,22] <- values(plc_raster)
rdn[,22] <- values(plc_raster)
rsn[,22] <- values(plc_raster)
pek[,22] <- values(plc_raster)
ear[,22] <- values(plc_raster)


setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/')
save(sos,eos,gsl1,gsl2,rup,rdn,rsn,panel_grid,pek,ear,
     file=paste('metrics_',tile_list[tt],'.rda',sep=''))  



# ################################################################################################
# setwd('/projectnb/modislc/users/mkmoon/NEphenology/data/c6/')
# path <- '/projectnb/modislc/users/mkmoon/NEphenology/data/c6'
# sstr <- '*pa*'
# file <- list.files(path=path,pattern=glob2rx(sstr),full.names=T)
# 
# stat <- matrix(NA,18,2)
# for(j in 1:18){
#   for(i in 1:15){
#     load(file[i])
#     if(i==1){
#       temp <- sos
#     }else{
#       temp <- rbind(temp,sos)
#     }
#   }
#   sos <- temp[which(temp[,17]>=(j*5)),1:16]
#   
#   par(mfrow=c(2,1),oma=c(2,2,2,2),mar=c(4,4,2,1))
#   # SOS trend
#   x <- 1:16
#   y <- apply(sos,2,mean,na.rm=T)
#   plot(x,y,axe=F,ann=F)
#   z <- zyp.sen(y~x)
#   abline(z$coefficients[1],z$coefficients[2],col=2)
#   MannKendall(y)
#   title(paste('Portion of forest cover within panel: ',j*5,'%',sep=''),cex.main=1.5,line=0.5)
#   text(13,max(y),paste('Trend = ',round(z$coefficients[2],3),' days/year',sep=''),cex=1.5,pos=1)
#   axis(1,at=seq(1,16,2),seq(2001,2016,2),cex.axis=1.5)
#   axis(2,at=seq(0,200,5),cex.axis=1.5)
#   box(lty=1)
#   mtext('SOS (day of year)',2,line=2.5,cex=1.5)
#   stat[j,1] <- z$coefficients[2]
#   
#   
#   # Geographical variation
#   x <- 1:16
#   y <- apply(sos,2,sd,na.rm=T)
#   plot(x,y,axe=F,ann=F)
#   z <- zyp.sen(y~x)
#   abline(z$coefficients[1],z$coefficients[2],col=2)
#   MannKendall(y)
#   text(13,max(y),paste('Trend = ',round(z$coefficients[2],3),' days/year',sep=''),cex=1.5,pos=1)
#   axis(1,at=seq(1,16,2),seq(2001,2016,2),cex.axis=1.5)
#   axis(2,at=seq(0,200,2),cex.axis=1.5)
#   box(lty=1)
#   mtext('SD (days)',2,line=2.5,cex=1.5)
#   stat[j,2] <- z$coefficients[2]
# 
# print(j)
# }
# 
# par(mfrow=c(2,1),oma=c(2,2,2,2),mar=c(4,4,1,1))
# plot(stat[,1],axe=F,ann=F)
# axis(1,at=seq(2,18,2),seq(10,90,10),cex.axis=1.5)
# axis(2,at=seq(-5,0,0.02),cex.axis=1.3)
# box(lty=1)
# mtext('Portion of forest cover within panel (%)',1,line=2.5,cex=1.5)
# mtext('SOS trend (days/year)',2,line=2.5,cex=1.5)
# 
# plot(stat[,2],axe=F,ann=F)
# axis(1,at=seq(2,18,2),seq(10,90,10),cex.axis=1.5)
# axis(2,at=seq(-5,0,0.01),cex.axis=1.3)
# box(lty=1)
# mtext('Portion of forest cover within panel (%)',1,line=2.5,cex=1.5)
# mtext('SD (days)',2,line=2.5,cex=1.5)


# num_cycles # 1
# evi_area_cycle1 # 2
# evi_amp_cycle1 # 3
# evi_min_cycle1 # 4
# ogi_cycle1 # 5
# midgup_cycle1 # 6
# mat_cycle1 # 7
# peak_cycle1 # 8
# sen_cycle1 # 9
# midgdown_cycle1 # 10
# dor_cycle1 # 11
# overall_qa_cycle1 # 12
# detailed_qa_cycle1 # 13
# evi_area_cycle2 # 14
# evi_amp_cycle2 # 15
# evi_min_cycle2 # 16
# ogi_cycle2 # 17
# midgup_cycle2 # 18
# mat_cycle2 # 19
# peak_cycle2 # 20
# sen_cycle2 # 21
# midgdown_cycle2 # 22
# dor_cycle2 # 23
# overall_qa_cycle2 # 24
# detailed_qa_cycle2=default_value # 25
