##########################################################################
# 
# Author: Mirco Migliavacca (version 1 2012/5/26)
# SCRIPT FOR THE EXTRACTION OF SPRING PHENOPHASES FROM HF DATASET
#
# Adapted by Emery R. Boose 2015/5/4 to use current file formats
# Minor updates ERB 2019/5/7
#
# Input file:
#     hf003-03-spring.csv
#
# Output files:
#     hf003-05-spring-mean-ind.csv (dates by year and tree)
#     hf003-06-spring-mean-spp.csv (dates by year and species)
# 
# Variables:
#     bb = budbreak
#     l75 = 75% leaf development
#     fbb = flowering budbreak
#     fopn = flowers open
#     sd.bb,sd.l75,sd.fbb,sd.fopn: standard deviations of dates
#
##########################################################################


spring <- function(dat, treecode) {	

  subdat <- subset(dat, tree.id==treecode)
  niterp <- 200000

  # start yearly loop	
  years <- unique(subdat$year)
  nyears <- length(years)
  i <- 1

  for (n in years) {
    #extract year		
    ydat <- subset(subdat, year==n)
		
    ifelse(i==1, tree.id <- treecode, tree.id <- c(tree.id, treecode))
		
    ifelse(i==1, species <- unlist(strsplit(treecode, "-"))[1], 
      species <- c(species, unlist(strsplit(treecode, "-"))[1]))
		
    if (length(ydat$bbrk[is.na(ydat$bbrk)==FALSE]) > 2) { 
      interp <- approx(ydat$doy, ydat$bbrk, n=niterp)
	
      ifelse(i==1, bb<-interp$x[interp$y >= 50][1], 
        bb <- c(bb, interp$x[interp$y >= 50][1]))	# 50% development
		}
		
    if (length(ydat$bbrk[is.na(ydat$bbrk)==FALSE]) <= 2) { 
      ifelse(i==1, bb <- NA, bb <- c(bb, NA))	# NA if values cannot be calculated
    }
		
    if (length(ydat$l75[is.na(ydat$l75)==FALSE]) > 2) { 
      interp <- approx(ydat$doy, ydat$l75, n=niterp)
      ifelse(i==1, l75 <- interp$x[interp$y >= 50][1],
        l75 <- c(l75, interp$x[interp$y >= 50][1]))	# 50% development (max 80)
    }
		
    if (length(ydat$l75[is.na(ydat$l75)==FALSE]) <= 2) { 
      ifelse(i==1, l75 <- NA, l75 <- c(l75, NA))	# NA if values cannot be calculated
    }
		
    if (length(ydat$fbud[is.na(ydat$fbud)==FALSE]) > 0) { 	
      day1 <- which(ydat$fbud=="BB")[1]	
      if (is.na(day1)==TRUE) { 
        ifelse(i==1, fbb <- NA, fbb <- c(fbb, NA))	# NA if values cannot be calculated
      }
			
      if (is.na(day1)==FALSE) { 
        ifelse(day1==1, tmpday <- ydat$doy[day1], tmpday <- (ydat$doy[day1]+ydat$doy[day1-1])/2)						
        ifelse(i==1, fbb <- tmpday, fbb <- c(fbb,fbb <- tmpday))	# flower budburst
      }
    }
		
    if (length(ydat$fbud[is.na(ydat$fbud)==FALSE]) == 0) { 
      ifelse(i==1, fbb <- NA, fbb <- c(fbb, NA))	# NA if values cannot be calculated
    }
		
    if (length(ydat$fopn[is.na(ydat$fopn)==FALSE]) > 2) { 
      interp <- approx(ydat$doy, ydat$fopn, n=niterp)
      ifelse(i==1, fopn <- interp$x[interp$y>=50][1],
        fopn <- c(fopn,interp$x[interp$y>=50][1]))	# 50% development (max 80)
    }
		
    if (length(ydat$fopn[is.na(ydat$fopn)==FALSE]) <= 2) { 
      ifelse(i==1, fopn <- NA, fopn <- c(fopn, NA))	# NA if values cannot be calculated
    }
		
    i = i + 1
  }	
	
  return(data.frame(year=years, tree.id=tree.id, species=species, bb.doy=bb, l75.doy=l75, fbb.doy=fbb, fopn.doy=fopn, stringsAsFactors = FALSE))	
}

file.in <- "/projectnb/modislc/users/mkmoon/TAscience/Data/hf/hf003-03-spring.csv"

file.out.tree <- "/projectnb/modislc/users/mkmoon/TAscience/Data/hf/hf003-05-spring-mean-ind.csv"
file.out.species <- "/projectnb/modislc/users/mkmoon/TAscience/Data/hf/hf003-06-spring-mean-spp.csv"

dat <- read.csv(file(file.in), header=TRUE, stringsAsFactors=FALSE)

tmp <- as.POSIXct(dat$date)
dat$year <- format(tmp, format="%Y")

# start tree ID loop

treecodes <- unique(dat$tree)
treecodes <- sort(treecodes)

for (loop in treecodes) {
  #extract year	
  out <- spring(dat, loop)
  ifelse(loop==treecodes[1], out.frame <- out, out.frame <- rbind(out.frame, out))
  cat("..Analyzing ", loop, "\n")
}

# convert julian day values to integer
out.frame.int <- out.frame
for (i in 4:7) {
  out.frame.int[ , i] <- round(out.frame.int[ , i])
}
tmp.frame <- out.frame


##########################################################################
hfsp <- unique(tmp.frame$species)
dat <- tmp.frame[tmp.frame$year>2000,]
library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

par(oma=c(2,2,2,2),mar=c(4,4,1,1))
k <- 1
for(i in 1:length(hfsp)){
  if(i==1){
    temp <- cbind(as.numeric(dat[dat$species==hfsp[1],1]),
                  as.numeric(dat[dat$species==hfsp[1],4]))
    if(length(unique(temp[,1]))>10){
      plot(temp,ylim=c(90,155),xlim=c(2001,2020),axes=F,ann=F,
           bg=col_vector[i],pch=21)  
      axis(1,c(2001,2005,2010,2015,2020),cex=1.8)
      axis(2,seq(80,200,10),cex=1.8)
      box()
      mtext('Bud Break dates (DOY)',2,line=2.5,cex=1.5)
      x <- temp[,1]; y <- temp[,2]
      abline(lm(y~x),lty=5,col=col_vector[i],lwd=2)  
    }
  }else{
    temp <- cbind(as.numeric(dat[dat$species==hfsp[i],1]),
                  as.numeric(dat[dat$species==hfsp[i],4]))
    if(length(unique(temp[,1]))>10){
      points(temp,ylim=c(80,170),bg=col_vector[i],pch=21)    
      x <- temp[,1]; y <- temp[,2]
      abline(lm(y~x),lty=5,col=col_vector[i],lwd=2)  
      k <- k+1
    }
  }
  print(k)
}

