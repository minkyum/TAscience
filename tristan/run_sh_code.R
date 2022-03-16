# Get NBAR data
# tiles <- c('h12v03','h12v04','h20v07','h20v08','h23v02','h27v05')
# tiles <- c('h11v09','h12v04','h21v02')
tiles <- c('h09v04','h09v05','h09v06','h10v03','h10v04',
           'h10v05','h10v06','h11v03','h11v04','h11v05',
           'h12v03','h12v04','h12v05','h13v03','h13v04',
           'h14v03','h14v04')
# tiles <- c('h10v02','h11v02','h12v02','h13v02','h14v02',
#            'h08v04','h08v05','h08v06','h09v02','h09v03')
tiles <- c('h11v09')

setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/e4ftl01.cr.usgs.gov/')
for(i in 1){
  for(year in 2014:2020){
    system(paste('qsub -N NB2',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_nbara2.sh ', tiles[i],year,sep=''))
    system(paste('qsub -N NB4',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_nbara4.sh ', tiles[i],year,sep=''))
  }  
}


# Get NBAR omitted
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/e4ftl01.cr.usgs.gov/')
for(i in 1){
  for(year in 2014:2020){
    system(paste('qsub -N om2_',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_nbara2_omit.sh ', tiles[i],year,sep=''))
    system(paste('qsub -N om4_',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_nbara4_omit.sh ', tiles[i],year,sep=''))
  }
}

# Run "change_folder_name.R" code
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/')
for(i in 1){
  for(year in 2014:2020){
    system(paste('qsub -N cfld',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_change_folder.sh ', tiles[i], year,sep=''))  
  }
}

# Run spline code 
tiles <- c('h12v04')

setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/h12v04')
for(i in 1){
  #Create output directory if it doesn't exist
  for(start_yr in 2017){
    # outDir <- paste0('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/spline_out/',tiles[i],'/',(start_yr+2))
    # if (!dir.exists(outDir)) {dir.create(outDir)}
    system(paste('qsub /projectnb/modislc/users/mkmoon/LCD_C6/C6_1/run_spline_fromJosh.sh ', tiles[i],' ',start_yr,sep=''))    
  }
}


# To see splined results
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/spline/ts/log/')
for(i in 1){
  for(year in 2019){
    for(j in 1:150){
      system(paste('qsub -N ck_sp -V -pe omp 4 -l h_rt=03:00:00 -l mem_total=32G /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_ck_spline.sh ',tiles[i],year,sep=''))  
    }
  }  
}


# Sub AnnualPhenology
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/des_sp/pheno/')
for(i in 1){
  #Create output directory if it doesn't exist
  # outDir <- paste0('/projectnb/modislc/users/mkmoon/LCD_C6/v7/pheno/bu/',tiles[i],'/')
  # if (!dir.exists(outDir)) {dir.create(outDir)}
  for(year in 2017){
    system(paste('qsub -V -pe omp 16 -l h_rt=06:00:00 -l mem_total=98G /usr3/graduate/mkmoon/GitHub/MCD12Q2/MCD12Q2C6_subAnnualPhenology_fromJosh.sh ',tiles[i],' ',year,sep=''))  
  }
}


###################
# Check NBAR
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/nbar/')
for(i in 1:3){
  system(paste('qsub -N NbCk',tiles[i],' -V -pe omp 2 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_check_nbar.sh ',tiles[i],sep=''))  
}


###################
# Get MCD12Q1 2020
setwd('/projectnb/modislc/projects/sat/data/mcd12q/q1')
for(year in 2002:2020){
  system(paste('qsub -N Q1',year,' -V -pe omp 8 -l h_rt=12:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_12q1.sh ',year,sep=''))  
}


###################
# Compare with C5 - quality check
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/C5_C6_1/')
for(i in 1:6){
    for(j in 1:10){
      system(paste('qsub -V -l h_rt=01:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare.sh ',i,j,sep=''))    
    }
}

# Compare with C5 - hist
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/hist/')
for(j in 1:10){
  system(paste('qsub -V -l h_rt=02:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare_hist.sh ',j,sep=''))    
}

# Compare with C5 - 1to1
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/1_to_1/')
for(j in 1:10){
  system(paste('qsub -V -l h_rt=02:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_compare_1to1.sh ',j,sep=''))    
}


#### Load C6 for phenocam comparison
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/load_c6/')
for(i in 1:2){
  system(paste('qsub -V -pe omp 4 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_load_c6.sh ',i,sep=''))    
}


#### MSLSP comparison
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/C6_1/figures/hls/')
for(i in 1:6){
  for(j in 1:4){
    system(paste('qsub -V -pe omp 2 -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_c6_hls.sh ',i,j,sep=''))      
  }
}




####################################################
# V7
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v7/')
system('qsub -V -pe omp 2 -j y -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/download_wget.sh')    


####################################################
# V8
setwd('/projectnb/modislc/users/mkmoon/LCD_C6/v8/figure/')
for(tt in 1:3){
  for(yy in 2003:2019){
    system(paste('qsub -V -pe omp 2 -j y -l h_rt=06:00:00 /usr3/graduate/mkmoon/GitHub/MCD12Q2/run_v8.sh ',tt,yy,sep=''))        
  }
}

