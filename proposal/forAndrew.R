

tsos <- matrix(NA,3,dim(phemet2)[2])

for(i in 1:dim(phemet2)[2]){
  x <- 1:20
  y1 <- phemet2[,i]
  y2 <- phemet3[,i]
  y3 <- phemet3[,i] - phemet2[,i]
  z1 <- zyp.sen(y1~x)
  z2 <- zyp.sen(y2~x)
  z3 <- zyp.sen(y3~x)
  tsos[1,i] <- z1$coefficients[2]   
  tsos[2,i] <- z2$coefficients[2]   
  tsos[3,i] <- z3$coefficients[2]   
}

dat <- rbind(round(phemet2),round(phemet3),round(aa),tsos)
colnames(dat) <- c('pix1','pix2','pix3','pix4','pix5',
                   'pix6','pix7','pix8','pix9','pix10',
                   'pix11','pix12','pix13','pix14','pix15',
                   'pix16','pix17','pix18','pix19','pix20',
                   'pix21','pix22','pix23','pix24','pix25',
                   'pix26','pix27','pix28','pix29','pix30',
                   'pix31','pix32')

rownames(dat) <- c('sos_2001','sos_2002','sos_2003','sos_2004','sos_2005','sos_2006','sos_2007','sos_2008','sos_2009','sos_2010',
                   'sos_2011','sos_2012','sos_2013','sos_2014','sos_2015','sos_2016','sos_2017','sos_2018','sos_2019','sos_2020',
                   'eos_2001','eos_2002','eos_2003','eos_2004','eos_2005','eos_2006','eos_2007','eos_2008','eos_2009','eos_2010',
                   'eos_2011','eos_2012','eos_2013','eos_2014','eos_2015','eos_2016','eos_2017','eos_2018','eos_2019','eos_2020',
                   'gsl_2001','gsl_2002','gsl_2003','gsl_2004','gsl_2005','gsl_2006','gsl_2007','gsl_2008','gsl_2009','gsl_2010',
                   'gsl_2011','gsl_2012','gsl_2013','gsl_2014','gsl_2015','gsl_2016','gsl_2017','gsl_2018','gsl_2019','gsl_2020',
                   'TheilSen_sos','TheilSen_eos','TheilSen_gsl')
dat <- as.data.frame(dat)

write.csv(dat,file='/projectnb/modislc/users/mkmoon/TAscience/Data/stat/trend_hf.csv')

