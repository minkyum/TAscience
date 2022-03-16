tiles <- c('h09v04','h09v05','h09v06',
           'h10v02','h10v03','h10v04','h10v05','h10v06',
           'h11v02','h11v03','h11v04','h11v05',
           'h12v02','h12v03','h12v04','h12v05',
           'h13v02','h13v03','h13v04',
           'h14v02','h14v03','h14v04',
           'h08v04','h08v05','h08v06','h09v02','h09v03')

#################
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/')
for(i in 20:24){
  system(paste('qsub -V -j y -pe omp 2 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_001.sh ',i,sep=''))          
}


#################
tiles <- c('h09v04','h09v05','h09v06','h10v02','h10v03',
           'h10v04','h10v05','h10v06','h11v02','h11v03',
           'h11v04','h11v05','h12v02','h12v03','h12v04',
           'h12v05','h13v02','h13v03','h13v04','h14v02',
           'h14v03','h14v04','h08v04','h08v05','h08v06','h09v02','h09v03')

setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/')
for(i in 25:27){ # 5; 13; 2; 7; 15; --- 1; 2; 3; 4; 5;    6; 8-10; 11,12,14,16; 17; 18; 19-21; 22; 23; 24
  for(year in c(2008:2015)){
    system(paste('qsub -N DN_',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_004.sh ',tiles[i],year,sep=''))            
  }
} 

# Get omitted
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/e4ftl01.cr.usgs.gov/')
for(i in c(25:27)){
  for(year in c(2008:2015)){
    system(paste('qsub -N om_',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_omit.sh ', tiles[i],year,sep=''))
  }
}


#################
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/met/')
for(i in c(1:27)){
  for(year in 2008:2015){
    system(paste('qsub -V -j y -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_005.sh ',tiles[i],year,sep=''))              
  }
}

setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/reg/')
for(i in c(1:27)){
  system(paste('qsub -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_006.sh ',tiles[i],sep=''))              
}

setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/stat/ano_flux')
for(ss in 1:24){
  system(paste('qsub -V -j y -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_007.sh ',ss,sep=''))              
}



setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/metrics/')
for(i in c(1:6)){
  system(paste('qsub -V -j y -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_008.sh ',i,sep=''))              
}


#############################
setwd('/projectnb/modislc/users/mkmoon/TAscience/Data/forJ/')
for(j in 2:8){ # metrics 2-8
  for(i in 1:166){ # sites
    ii <- sprintf('%03d',i)
    jj <- sprintf('%01d',j)
    system(paste('qsub -V -j y -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/TAscience/run_script_009.sh ',ii,jj,sep=''))                
  }
}

