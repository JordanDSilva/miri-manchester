library(celestial)
library(ProSpect)
library(Rfits)
library(data.table)
library(magicaxis)
library(foreach)
library(Highlander)
library(doParallel)

# calculate Vmax

main_dir = "~/Documents/Epochs/"
fname = "smiles"

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))
phot_cat = Rfits_read(
  paste0(main_dir, "/data/catalogues/", fname, "_profound_multiband_cat.fits")
)
prospect_cat = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
))
prospect_photom = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_fit_photom.csv"
))

Vmax_mag = function(parm, Data){
  
  z=Data$z
  wave=Data$wave
  lum=Data$lum
  ApMag_lim=Data$ApMag_lim
  filter=Data$filter
  
  photom_mag_values = photom_lum(
    wave = wave, 
    lum = lum, 
    z = z+parm, 
    outtype = 'magAB', 
    filters = filter,
    ref = "737"
  ) 
  
  photom_flux_values = 10^( -0.4 * (photom_mag_values - 8.9) ) * Inoue14_IGM(wave = c(56076.384, 76128.766), z = z+parm, Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
  
  mean_flux = mean(photom_flux_values, na.rm = TRUE)
  mean_mag = -2.5*log10(mean_flux) + 8.9
  
  return(mean_mag)
}

dzmax = vector(length = dim(prospect_cat)[1] )
for(i in 1:dim(prospect_cat)[1]){
# for(i in 1:10){
  if(i %% 100 == 0){
    message(i)
  }
  
  SKIP = FALSE
  if(prospect_cat$FitType[i] == "StellarAGN"){
    sed_stub = paste0(
      "~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
    file_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
  } else if(prospect_cat$FitType[i] == "Stellar"){
    sed_stub = paste0(
      "~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
    file_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
  } else {
    SKIP = TRUE
  }
  
  if(file.exists(sed_stub) & SKIP == FALSE){
    fit_sed = readRDS(
      sed_stub
    )
    
    fit = readRDS(
      file_stub
    )
    
    det_limit = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9
    Data = list(
      wave = fit_sed$restWave,
      lum = fit_sed$lumQ50,
      z = prospect_cat$zQ50[i], 
      ApMag_lim = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9,
      filter = list(filter_transforms$JWST_MIRI_F560W, filter_transforms$JWST_MIRI_F770W)
    )
    
    parm_zgrid = seq(
      0.0, 40.0, 0.01
    )
    
    registerDoParallel(cores = 8)
    redshift_mags = foreach(z = parm_zgrid, .combine = c) %dopar% {
      Vmax_mag(z, Data)
    }
    stopImplicitCluster()
    
    if(sum(redshift_mags >= det_limit) == 0){
      dzmax_test = 100
    }else{
      dzmax_test = parm_zgrid[redshift_mags >= det_limit][1]
    }
  }else{
    dzmax_test = 0
  }
  dzmax[i] = dzmax_test
}

# dzmax[is.na(dzmax)] = 100
df = data.frame(
  "zmax" = prospect_cat$zQ50 + dzmax,
  "dzmax" = dzmax
)

fwrite(
  df, 
  paste0(main_dir, "/data/depths/", "zmax_5sigma_data.csv")
)
