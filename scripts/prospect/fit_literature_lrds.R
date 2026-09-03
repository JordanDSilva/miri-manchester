# LRD forks
library('ProSpectData', lib.loc = "~/Documents/ProSpectStuff/ProSpectLRDLib")
library('ProSpect', lib.loc = "~/Documents/ProSpectStuff/ProSpectLRDLib")
library(celestial)
library(ParmOff)
library(Highlander)
library(foreach)
library(matrixStats)
library(doParallel)
library(magicaxis)
library(data.table)
# library(ProSpectUtils)
devtools::load_all("~/Documents/ProSpectStuff/ProSpectUtils/")


set.seed(666)

# main_dir = "/scratch/"

CMANiters = 2000
LDNiters = 2000

min_pc_err = 0.1

main_dir = "~/Documents/Epochs/"
# main_dir = "/foo/" # from docker/singularity container

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)

virgil_photometry = data.frame(
  Instrument = c("JWST/NIRCam", "HST/ACS_WFC", "HST/WFC3_IR", "JWST/NIRCam",
                 "HST/WFC3_IR", "HST/WFC3_IR", "JWST/NIRCam", "HST/WFC3_IR",
                 "JWST/NIRCam", "JWST/NIRCam", "JWST/NIRCam", "JWST/NIRCam",
                 "JWST/NIRCam", "JWST/NIRCam", "JWST/NIRCam", "JWST/NIRCam",
                 "JWST/NIRCam", "JWST/NIRCam", "JWST/NIRCam", "JWST/MIRI",
                 "JWST/MIRI",   "JWST/MIRI",   "JWST/MIRI",  "JWST/MIRI"),
  Filter = c("F090W", "F850LP", "F105W", "F115W", "F125W", "F140W", "F150W",
             "F160W", "F182M", "F200W", "F210M", "F277W", "F335M", "F356W",
             "F410M", "F430M", "F444W", "F460M", "F480M", "F560W", "F770W",
             "F1000W", "F1280W", "F1500W"),
  Flux_Jy = c(12.86,  3.38, 22.07, 30.54, 27.06, 22.90, 28.92, 28.52,
              24.63, 30.32, 36.70, 40.11, 38.88, 93.84, 37.01, 37.79,
              49.33, 46.99, 67.58, 107.70, 248.20, 399.10, 382.70, 379.50)*1e-9,
  FluxErr_Jy = c(2.25, 2.70, 1.95, 2.00, 2.51, 2.59, 1.85, 2.82,
                 2.65, 2.04, 3.10, 1.85, 2.26, 4.32, 1.89, 4.80,
                 2.27, 6.25, 7.91, 3.03, 30.13, 23.79, 161.10, 296.20)*1e-9
)
virgil_photometry = virgil_photometry[virgil_photometry$Instrument != "HST/WFC3_IR", ]
virgil_photometry$cenwave = c(9024.004, 9031.483, 11540.782, 15008.986, 18452.162, 19876.006, 20955.949, 27762.299, 33622.534, 35654.890, 40829.883, 42812.955, 44017.433, 46296.169, 48165.464, 56076.384, 76128.766, 99568.560, 128160.125, 267365.990)
virgil_Data = list(
  "fluxJy" = virgil_photometry$Flux_Jy,
  "fluxerrJy" = virgil_photometry$FluxErr_Jy,
  "cenwave" = virgil_photometry$cenwave,
  "z" = 6.6312,
  "LumDist_Mpc" = cosdistLumDist(z = 6.6312, ref = "737"),
  "filters" = filter_transforms[c("JWST_NIRCAM_F090W", "HST_ACSWFC_F850LP", "JWST_NIRCAM_F115W", "JWST_NIRCAM_F150W", 
                                  "JWST_NIRCAM_F182M", "JWST_NIRCAM_F200W", "JWST_NIRCAM_F210M", "JWST_NIRCAM_F277W", 
                                  "JWST_NIRCAM_F335M", "JWST_NIRCAM_F356W", "JWST_NIRCAM_F410M", "JWST_NIRCAM_F430M",
                                  "JWST_NIRCAM_F444W", "JWST_NIRCAM_F460M", "JWST_NIRCAM_F480M", "JWST_MIRI_F560W",
                                  "JWST_MIRI_F770W", "JWST_MIRI_F1000W", "JWST_MIRI_F1280W", "JWST_MIRI_F1500W")]
)
virgil_Data$name = "Virgil"

cliff_photometry <- data.frame(
  filter = c("NIRCam/F090W", "NIRCam/F115W", "NIRCam/F150W", "NIRCam/F200W",
             "NIRCam/F277W", "NIRCam/F356W", "NIRCam/F410M", "NIRCam/F444W",
             "MIRI/F770W", "MIRI/F1800W"),
  f_nu_uJy = c(0.0148, 0.028, 0.0269, 0.322, 1.090, 1.334, 1.597, 1.725, 1.54, 1.96),
  err_upper = c(0.0024, 0.004, 0.0022, 0.029, 0.004, 0.005, 0.011, 0.012, 0.03, 0.5),
  err_lower = c(0.0019, 0.004, 0.0017, 0.028, 0.004, 0.006, 0.013, 0.013, 0.03, 0.5), 
  cenwave =   c(9024.004,11540.782, 15008.986, 19876.006, 27762.299, 35654.890, 40829.883, 44017.433, 76128.766, 180396.821)
)
cliff_Data = list(
  "fluxJy" = cliff_photometry$f_nu_uJy*1e-6,
  "fluxerrJy" = 0.5*(cliff_photometry$err_lower + cliff_photometry$err_upper)*1e-6,
  "cenwave" = cliff_photometry$cenwave,
  "z" = 3.55, 
  "LumDist_Mpc" = cosdistLumDist(z = 3.55, ref = "737"),
  "filters" = filter_transforms[c("JWST_NIRCAM_F090W", "JWST_NIRCAM_F115W", "JWST_NIRCAM_F150W", "JWST_NIRCAM_F200W",
                                  "JWST_NIRCAM_F277W", "JWST_NIRCAM_F356W", "JWST_NIRCAM_F410M", "JWST_NIRCAM_F444W",
                                  "JWST_MIRI_F770W", "JWST_MIRI_F1800W")]
)
cliff_Data$name = "Cliff"

mom_photometry <- data.frame(
  Filter = c("F090W", "F115W", "F150W", "F200W", "F277W",
             "F356W", "F410M", "F444W", "F770W", "F1800W"),
  f_nu_uJy    = c(11, 15,  4,  5,  8,  93, 206, 243,  465, 1269),
  err = c( 8,  8,  6,  5,  4,   5,  10,  12,   52,  615),
  cenwave =   c(9024.004,11540.782, 15008.986, 19876.006, 27762.299, 35654.890, 40829.883, 44017.433, 76128.766, 180396.821)
)
mom_Data = list(
  "fluxJy" = mom_photometry$f_nu_uJy*1e-9,
  "fluxerrJy" = mom_photometry$err*1e-9,
  "cenwave" = mom_photometry$cenwave,
  "z" = 7.7569,
  "LumDist_Mpc" = cosdistLumDist(z = 7.7569, ref = "737"),
  "filters" = filter_transforms[c("JWST_NIRCAM_F090W", "JWST_NIRCAM_F115W", "JWST_NIRCAM_F150W", "JWST_NIRCAM_F200W",
                                  "JWST_NIRCAM_F277W", "JWST_NIRCAM_F356W", "JWST_NIRCAM_F410M", "JWST_NIRCAM_F444W",
                                  "JWST_MIRI_F770W", "JWST_MIRI_F1800W")]
)
mom_Data$name = "MoMBH1"

pfunc_no_tau_birth=function(parm){
  prior = sum(
    c(
      100*erf(-5 * (parm['tau_screen'] - parm['tau_birth'])) - 100, ## Column density of dust in the ISM should be less than that of the birth clouds
      dnorm(parm['alpha_SF_birth'], mean=2, sd=1, log=TRUE),
      dnorm(parm['alpha_SF_screen'], mean=2, sd=1, log=TRUE),
      0
    ),
    na.rm = TRUE
  )
  return(prior)
}
constraintfunc = function(parm, Data){
  if(parm[grep('tau_screen', Data$parm.names)] > parm[grep('tau_birth', Data$parm.names)]){ ## V-band optical depth in the ISM must be less than the birth clouds
    parm[grep('tau_screen', Data$parm.names)] = parm[grep('tau_birth', Data$parm.names)]
  }
  return(parm)
}

do_photoz_bc03_snorm_fritz_norm_eazy = function(galData){
  
  save_stub = paste0(main_dir, "/data/literature_lrds/", galData$name, "_fritz_fit.rds")
  print(save_stub)
  
  flux_data = data.frame(
    "filter" = names(galData$filters),
    "cenwave" = galData$cenwave,
    "flux" = as.numeric(galData$fluxJy)
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(galData$fluxerrJy)^2 + (min_pc_err * flux_data$flux)^2 )
  
  filtout = galData$filters
  
  Redshift = galData$z
  agemax = cosdistUniAgeAtz(z = Redshift, ref = "737") * 1e9
  magemax = agemax/1e9
  Zagemax = agemax/1e9
  LumDist_Mpc = cosdistLumDist(z = Redshift, ref = "737")
  
  Data = list(
    flux = flux_data,
    mode = "photom",
    arglist = list(
      z = Redshift,
      agemax = agemax,
      magemax = magemax, 
      Zagemax = Zagemax,
      LumDist_Mpc = LumDist_Mpc,
      
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      
      emission = TRUE,
      photoz = FALSE,
      IGMabsorb = "Inoue14",
      
      mtrunc = 2,
      
      AGNrm = 60,
      AGNbe = -0.5,
      AGNal = 4.0,
      AGNct = 100,
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7 
    ),
    
    speclib = BC03hr,
    Dale = Dale_NormTot,
    filtout = filtout,
    Inoue14_LAFcoef = Inoue14_LAFcoef,
    Inoue14_DLAcoef = Inoue14_DLAcoef,
    
    SFH=SFHfunc,
    AGN = Fritz,
    Dale_M2L_func = Dale_M2L_func,
    
    parm.names = c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta'), 
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst"),
    logged = c(
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE,
      TRUE, FALSE, TRUE
    ),
    intervals = list(
      lo = c(
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0, 
        -4, 
        40, 0.001, -1
      ), 
      hi = c(
        4, magemax, 2, 1, 
        1.5, 1, 4, 4, 
        -1.3, 
        49, 89.990, 1
      )
    ), 
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_no_tau_birth,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  
  start = Sys.time()
  highout = Highlander(
    startpoint, 
    Data, 
    ProSpectSEDlike,
    Niters = c(CMANiters, LDNiters), 
    NfinalMCMC = LDNiters,
    optim_iters = 2, 
    seed = 666,
    lower=Data$intervals$lo, 
    upper=Data$intervals$hi
  )
  print(Sys.time() - start)
  
  Data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=Data)
  
  out = list(
    save_stub = save_stub,
    "bestfit" = bestfit,
    "highout" = highout,
    "galInfo" = cbind( "galName" = galData$name, "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}
do_photoz_bc03_snorm_LRDMBB_norm_eazy_fit = function(galData){
  save_stub = paste0(main_dir, "/data/literature_lrds/", galData$name, "_lrdmbb_fit.rds")
  print(save_stub)
  
  flux_data = data.frame(
    "filter" = names(galData$filters),
    "cenwave" = galData$cenwave,
    "flux" = as.numeric(galData$fluxJy)
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(galData$fluxerrJy)^2 + (min_pc_err * flux_data$flux)^2 )
  
  filtout = galData$filters
  
  Redshift = galData$z
  agemax = cosdistUniAgeAtz(z = Redshift, ref = "737") * 1e9
  magemax = agemax/1e9
  Zagemax = agemax/1e9
  LumDist_Mpc = cosdistLumDist(z = Redshift, ref = "737")
  
  Data = list(
    flux = flux_data,
    mode = "photom",
    arglist = list(
      z = Redshift,
      agemax = agemax,
      magemax = magemax, 
      Zagemax = Zagemax,
      LumDist_Mpc = LumDist_Mpc,
      
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      
      emission = TRUE,
      photoz = FALSE,
      IGMabsorb = "Inoue14",
      
      mtrunc = 2,
      
      AGNrm = 60,
      AGNct = 160,
      AGNq = 0,
      AGNp = 1,
      AGNta = 3,
      AGNan = 45,
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7
    ),
    
    speclib = BC03hr,
    Dale = Dale_NormTot,
    filtout = filtout,
    
    SFH=SFHfunc,
    AGN = LRDBBSKIRTOR,
    Dale_M2L_func = Dale_M2L_func,
    parm.names = c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal', 'AGNlum', 'AGNteff', 'AGNbeta', 'tau_AGN'), 
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst"),
    logged = c(
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE,
      TRUE, TRUE, FALSE, TRUE
    ),
    intervals = list(
      lo = c(
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0, 
        -4, 
        40, log10(2000), 0, -5.0
      ), 
      hi = c(
        4, magemax, 2, 1, 
        1.5, 1, 4, 4, 
        -1.3, 
        49, log10(10000), 3, 1.5
      )
    ), 
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_no_tau_birth,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  
  highout = Highlander(
    startpoint, 
    Data, 
    ProSpectSEDlike,
    Niters = c(CMANiters, LDNiters), 
    NfinalMCMC = LDNiters,
    optim_iters = 2, 
    seed = 666,
    lower=Data$intervals$lo, 
    upper=Data$intervals$hi
  )
  
  Data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=Data)
  
  out = list(
    save_stub = save_stub,
    "bestfit" = bestfit,
    "highout" = highout,
    "galInfo" = cbind( "galName" = galData$name, "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}
do_photoz_bc03_snorm_LRDLIU_norm_eazy_fit = function(galData){
  save_stub = paste0(main_dir, "/data/literature_lrds/", galData$name, "_lrdliu_fit.rds")
  print(save_stub)
  
  flux_data = data.frame(
    "filter" = names(galData$filters),
    "cenwave" = galData$cenwave,
    "flux" = as.numeric(galData$fluxJy)
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(galData$fluxerrJy)^2 + (min_pc_err * flux_data$flux)^2 )
  
  filtout = galData$filters
  
  Redshift = galData$z
  agemax = cosdistUniAgeAtz(z = Redshift, ref = "737") * 1e9
  magemax = agemax/1e9
  Zagemax = agemax/1e9
  LumDist_Mpc = cosdistLumDist(z = Redshift, ref = "737")
  
  Data = list(
    flux = flux_data,
    mode = "photom",
    arglist = list(
      z = Redshift,
      agemax = agemax,
      magemax = magemax, 
      Zagemax = Zagemax,
      LumDist_Mpc = LumDist_Mpc,
      
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      
      emission = TRUE,
      photoz = FALSE,
      IGMabsorb = "Inoue14",
      
      mtrunc = 2,

      AGNrm = 60,
      AGNct = 160,
      AGNq = 0,
      AGNp = 1,
      AGNta = 3,
      AGNan = 45,

      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7
    ),

    speclib = BC03hr,
    Dale = Dale_NormTot,
    filtout = filtout,

    SFH=SFHfunc,
    AGN = LRDLIU,
    Dale_M2L_func = Dale_M2L_func,
    parm.names = c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal', 'AGNlum', 'AGNteff', 'AGNlogg', 'tau_AGN'),
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst"),
    logged = c(
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE,
      TRUE, TRUE, FALSE, TRUE
    ),
    intervals = list(
      lo = c(
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0,
        -4,
        40, log10(2000), -4.0, -5.0
      ),
      hi = c(
        4, magemax, 2, 1,
        1.5, 1, 4, 4,
        -1.3,
        49, log10(10000), 1.5, 1.5
      )
    ),
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_no_tau_birth,
    constraints = NULL,
    verbose = FALSE
  )

  Data$constraints = function(parm){constraintfunc(parm, Data)}

  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44

  highout = Highlander(
    startpoint,
    Data,
    ProSpectSEDlike,
    Niters = c(CMANiters, LDNiters),
    NfinalMCMC = LDNiters,
    optim_iters = 2,
    seed = 666,
    lower=Data$intervals$lo,
    upper=Data$intervals$hi
  )

  Data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=Data)

  out = list(
    save_stub = save_stub,
    "bestfit" = bestfit,
    "highout" = highout,
    "galInfo" = cbind( "galName" = galData$name, "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}
do_photoz_bc03_snorm_sirocco_norm_eazy_fit = function(galData){
  save_stub = paste0(main_dir, "/data/literature_lrds/", galData$name, "_sirocco_fit.rds")
  print(save_stub)
  
  flux_data = data.frame(
    "filter" = names(galData$filters),
    "cenwave" = galData$cenwave,
    "flux" = as.numeric(galData$fluxJy)
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(galData$fluxerrJy)^2 + (min_pc_err * flux_data$flux)^2 )
  
  filtout = galData$filters
  
  Redshift = galData$z
  agemax = cosdistUniAgeAtz(z = Redshift, ref = "737") * 1e9
  magemax = agemax/1e9
  Zagemax = agemax/1e9
  LumDist_Mpc = cosdistLumDist(z = Redshift, ref = "737")
  
  Data = list(
    flux = flux_data,
    mode = "photom",
    arglist = list(
      z = Redshift,
      agemax = agemax,
      magemax = magemax, 
      Zagemax = Zagemax,
      LumDist_Mpc = LumDist_Mpc,
      
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      
      emission = TRUE,
      photoz = FALSE,
      IGMabsorb = "Inoue14",
      
      mtrunc = 2,
      
      AGNrm = 60,
      AGNct = 160,
      AGNq = 0,
      AGNp = 1,
      AGNta = 3,
      AGNan = 45,
      AGNZ = -1,
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7
    ),
    
    speclib = BC03hr,
    Dale = Dale_NormTot,
    filtout = filtout,
    
    SFH=SFHfunc,
    AGN = SIROCCO,
    Dale_M2L_func = Dale_M2L_func,
    parm.names = c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal', 'AGNlum', 'AGNmdot', 'tau_AGN'),
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst"),
    logged = c(
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE,
      TRUE, FALSE, TRUE
    ),
    intervals = list(
      lo = c(
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0,
        -4,
        40, 0.5, -5.0
      ),
      hi = c(
        4, magemax, 2, 1,
        1.5, 1, 4, 4,
        -1.3,
        49, 8, 1.5
      )
    ),
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_no_tau_birth,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  
  highout = Highlander(
    startpoint,
    Data,
    ProSpectSEDlike,
    Niters = c(CMANiters, LDNiters),
    NfinalMCMC = LDNiters,
    optim_iters = 2,
    seed = 666,
    lower=Data$intervals$lo,
    upper=Data$intervals$hi
  )
  
  Data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=Data)
  
  out = list(
    save_stub = save_stub,
    "bestfit" = bestfit,
    "highout" = highout,
    "galInfo" = cbind( "galName" = galData$name, "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}

loop_through_fit = function(galData){
  for(do_fit in list(
    # do_photoz_bc03_snorm_fritz_norm_eazy,
    # do_photoz_bc03_snorm_LRDMBB_norm_eazy_fit,
    # do_photoz_bc03_snorm_LRDLIU_norm_eazy_fit,
    do_photoz_bc03_snorm_sirocco_norm_eazy_fit
  )){
  # for(do_fit in list(do_photoz_bc03_snorm_fritz_norm_eazy)){
    fit = do_fit(galData)
    
    summ = suppressMessages(
      ProSpectUtils::calc_astro(
        fit$bestfit,
        fit$highout,
        newDust = TRUE,
        specz = TRUE,
        cores = ncoresPost,
      )
    )
    
    fit$astro = summ$summary
    fit$post = summ$Posterior
    
    fit$bestfit$Data$prior = NULL
    fit$bestfit$Data$constraints = NULL
    fit$bestfit$Data$speclib = NULL
    fit$bestfit$Data$AGN = NULL
    fit$bestfit$Data$Dale = NULL
    
    saveRDS(
      object = fit,
      file = fit$save_stub
    )
  }
}

# loop_through_fit(virgil_Data)
loop_through_fit(cliff_Data)
loop_through_fit(mom_Data)


plot_spectrum_ribbon <- function(spectrum,
                                 key,
                                 xvar   = "restWave",
                                 col    = "black",
                                 alpha  = 0.3,
                                 add    = FALSE,
                                 ylim   = NULL,
                                 xlim   = NULL,
                                 xlab   = "Wavelength",
                                 ylab   = "Luminosity",
                                 ...) {
  
  # build the Q16/Q50/Q84/Q95 column names from the key
  q16_name <- paste0(key, "Q16")
  q50_name <- paste0(key, "Q50")
  q84_name <- paste0(key, "Q84")
  q95_name <- paste0(key, "Q95")   # used only for auto ylim
  
  x   <- spectrum[[xvar]]
  q16 <- spectrum[[q16_name]]
  q50 <- spectrum[[q50_name]]
  q84 <- spectrum[[q84_name]]
  
  if (is.null(x) || is.null(q16) || is.null(q50) || is.null(q84)) {
    stop(sprintf("Missing one or more required columns: %s, %s, %s, %s",
                 xvar, q16_name, q50_name, q84_name))
  }
  
  # keep only rows that are finite (drops NA/NaN/Inf) AND strictly positive
  # (needed both to filter bad data and because log-scale can't plot 0/negative anyway)
  keep <- is.finite(x)   & x   > 0 &
    is.finite(q16) & q16 > 0 & 
    is.finite(q50) & q50 > 0 & 
    is.finite(q84) & q84 > 0 
  
  x   <- x[keep]
  q16 <- q16[keep]
  q50 <- q50[keep]
  q84 <- q84[keep]
  
  if (length(x) == 0) {
    stop(sprintf("No valid (finite, positive) data remaining for key '%s' after filtering", key))
  }
  
  # default axis ranges
  if (is.null(xlim)) xlim <- range(x, na.rm = TRUE)
  if (is.null(ylim)) {
    q95 <- spectrum[[q95_name]]
    upper_ref <- if (!is.null(q95)) q95 else q84
    upper_ref <- upper_ref[is.finite(upper_ref) & upper_ref > 0]
    ylim <- pmax(1e-1, range(upper_ref, na.rm = TRUE))
  }
  
  # convert named colour to RGBA with transparency for the fill
  poly_col <- adjustcolor(col, alpha.f = alpha)
  
  if (!add) {
    plot(NA, xlim = xlim, ylim = ylim, log = "xy",
         xlab = xlab, ylab = ylab, ...)
  }
  
  # shaded Q16-Q84 ribbon as a polygon
  polygon(
    c(x, rev(x)),
    c(q16, rev(q84)),
    col    = poly_col,
    border = NA
  )
  
  # median line on top
  lines(x, q50, col = adjustcolor(col, 0.5), lwd = 2)
  
  invisible(list(x = x, q16 = q16, q50 = q50, q84 = q84))
}

plot_sfh_zh <- function(spectrum,
                        key,
                        xvar,
                        col     = "black",
                        alpha   = 0.3,
                        ylab    = NULL,
                        ylab2   = NULL,
                        xlab    = "Age (Gyr)",
                        add     = FALSE,
                        lwd     = 2,
                        ...) {
  
  # build column names from the key
  q16_name <- paste0(key, "Q16")
  q50_name <- paste0(key, "Q50")
  q84_name <- paste0(key, "Q84")
  
  x   <- spectrum[[xvar]]
  q16 <- spectrum[[q16_name]]
  q50 <- spectrum[[q50_name]]
  q84 <- spectrum[[q84_name]]
  
  if (is.null(x) || is.null(q16) || is.null(q50) || is.null(q84)) {
    stop(sprintf("Missing required column(s): %s, %s, %s, %s",
                 xvar, q16_name, q50_name, q84_name))
  }
  
  # drop NA/Inf rows across all series, keep everything aligned
  keep <- is.finite(x) & is.finite(q16) & is.finite(q50) & is.finite(q84)
  x   <- x[keep]
  q16 <- q16[keep]
  q50 <- q50[keep]
  q84 <- q84[keep]
  
  if (length(x) == 0) {
    stop(sprintf("No valid data remaining for key '%s' after filtering", key))
  }
  
  poly_col <- adjustcolor(col, alpha.f = alpha)
  
  if (!add) {
    # first call: draw primary series + ribbon on the standard left axis
    ylim <- range(c(q16, q84), na.rm = TRUE)
    
    plot(
      x, q50,
      type = 'l',
      lwd  = lwd,
      col  = col,
      xlab = xlab,
      ylab = ylab,
      ylim = ylim,
      ...
    )
    
    polygon(
      c(x, rev(x)),
      c(q16, rev(q84)),
      col    = poly_col,
      border = NA
    )
    lines(x, q50, col = col, lwd = lwd)   # redraw median on top of ribbon
    
  } else {
    # second call: overlay this series + ribbon on a new right-hand axis
    y_max <- max(q84, na.rm = TRUE)
    y_min <- min(0, min(q16, na.rm = TRUE))
    
    par(usr = c(
      par()$usr[1:2],
      y_min - y_max * 0.04,
      y_max * 1.04
    ))
    
    polygon(
      c(x, rev(x)),
      c(q16, rev(q84)),
      col    = poly_col,
      border = NA
    )
    lines(x, q50, col = col, lwd = lwd)
    axis(4, col = col, col.axis = col)
    
    if (!is.null(ylab2)) {
      mtext(ylab2, side = 4, line = 3, col = col)
    }
  }
  
  invisible(list(x = x, q16 = q16, q50 = q50, q84 = q84))
}

compute_fit_sed = function(file_stub, AGN_model, fittype){
  out = readRDS(file_stub)
  
  save_stub = gsub("fit.rds", "sed.rds", file_stub)
  
  post = out$post
  
  Data = out$bestfit$Data
  Data$speclib = BC03hr
  
  Data$AGN = AGN_model
  Data$Dale = Dale_NormTot
  
  lum_samples = foreach(j = 1:dim(post)[1]) %do% {
    
    post_sample = out$highout$LD_last$Posterior1[j, names(out$bestfit$parm)]
    out_sample = ProSpectSEDlike(
      post_sample,
      Data
    )
    
    if(is.null(Data$AGN)){
      AGN = rep(0,10)
    }else{
      AGN = out_sample$SEDout$AGN$lum
    }
    
    if(!is.null(Data$AGN)){
      photom_AGN_source = photom_lum(
        wave = out_sample$SEDout$AGN$wave,
        lum = out_sample$SEDout$AGN$lum,
        filters = Data$filtout, 
        outtype = "Jy", 
        z = Data$arglist$z, 
        ref = "737"
      )
      AGN_Dale_frac = sum( out_sample$SEDout$AGN$lum[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4] * c(0, diff(out_sample$SEDout$AGN$wave[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4])) ) / 
        sum( out_sample$SEDout$FinalLum$lum[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4] * c(0, diff(out_sample$SEDout$FinalLum$wave[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4])) )
    }else{
      photom_AGN_source = rep(0, length(Data$filtout))
      AGN_Dale_frac = 0
    }
    
    photom_all_source = photom_lum(
      wave = out_sample$SEDout$FinalLum$wave,
      lum = out_sample$SEDout$FinalLum$lum,
      filters = Data$filtout, 
      outtype = "Jy", 
      z = Data$arglist$z, 
      ref = "737"
    )
    
    photom_all = photom_all_source * Inoue14_IGM(wave = Data$flux$cenwave, z = Data$arglist$z, Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
    photom_AGN = photom_AGN_source *  Inoue14_IGM(wave = Data$flux$cenwave, z = Data$arglist$z, Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
    
    names(photom_AGN) = paste0(Data$flux$filter, "_fit_AGN_flux")
    names(photom_all) = paste0(Data$flux$filter, "_fit_flux")
    names(photom_AGN_source) = paste0(Data$flux$filter, "_source_fit_AGN_flux")
    names(photom_all_source) = paste0(Data$flux$filter, "_source_fit_flux")
    
    list(
      "rest" = out_sample$SEDout$FinalLum$lum,
      "StarsUnAtten" = out_sample$SEDout$StarsUnAtten$lum,
      "StarsAtten" = out_sample$SEDout$StarsAtten$lum,
      "DustEmit" = out_sample$SEDout$DustEmit$lum,
      "AGN" = AGN,
      "obs" = out_sample$SEDout$FinalFlux$flux, 
      "agevec" =  out_sample$SEDout$Stars$agevec/1e9,
      "sfhvec" = out_sample$SEDout$Stars$SFR,
      "Zvec" = out_sample$SEDout$Stars$Zvec,
      "photom" = c(
        photom_all,
        photom_AGN,
        photom_all_source,
        photom_AGN_source,
        "AGN_Dale_frac" = AGN_Dale_frac
      )
    )
  }
  
  rest_samples         = do.call(rbind, lapply(lum_samples, function(x) x$rest))
  StarsUnAtten_samples = do.call(rbind, lapply(lum_samples, function(x) x$StarsUnAtten))
  StarsAtten_samples   = do.call(rbind, lapply(lum_samples, function(x) x$StarsAtten))
  DustEmit_samples     = do.call(rbind, lapply(lum_samples, function(x) x$DustEmit))
  AGN_samples          = do.call(rbind, lapply(lum_samples, function(x) x$AGN))
  obs_samples          = do.call(rbind, lapply(lum_samples, function(x) x$obs))
  phot_samples         = do.call(rbind, lapply(lum_samples, function(x) x$photom))
  
  med_phot = colMedians(as.matrix(phot_samples))
  Q16_phot = colQuantiles(as.matrix(phot_samples), probs = 0.16)
  Q84_phot = colQuantiles(as.matrix(phot_samples), probs = 0.84)
  Q05_phot = colQuantiles(as.matrix(phot_samples), probs = 0.05)
  Q95_phot = colQuantiles(as.matrix(phot_samples), probs = 0.95)
  
  names(med_phot) = paste0(names(med_phot), "Q50")
  names(Q16_phot) = paste0(names(Q16_phot), "Q16")
  names(Q84_phot) = paste0(names(Q84_phot), "Q84")
  names(Q05_phot) = paste0(names(Q05_phot), "Q05")
  names(Q95_phot) = paste0(names(Q95_phot), "Q95")
  
  row_phot <- c(med_phot, Q16_phot, Q84_phot, Q05_phot, Q95_phot)
  
  max_age_len <- max(sapply(lum_samples, function(x) length(x$agevec)))
  agevec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
    v <- x$agevec
    length(v) <- max_age_len   # pads with NA automatically if v is shorter
    v
  }))
  sfhvec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
    v <- x$sfhvec
    length(v) <- max_age_len   # pads with NA automatically if v is shorter
    v
  }))
  Zvec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
    v <- x$Zvec
    length(v) <- max_age_len   # pads with NA automatically if v is shorter
    v
  }))
  
  spectrum = list()
  
  # fitted photom
  spectrum$phot = row_phot
  spectrum$phot_post = phot_samples
  
  # rest frame total
  spectrum$restWave     = out$bestfit$SEDout$FinalLum$wave
  spectrum$lumQ50       = colMedians(as.matrix(rest_samples), na.rm = TRUE)
  spectrum$lumQ16       = colQuantiles(as.matrix(rest_samples), probs = 0.16, na.rm = TRUE)
  spectrum$lumQ84       = colQuantiles(as.matrix(rest_samples), probs = 0.84, na.rm = TRUE)
  spectrum$lumQ05       = colQuantiles(as.matrix(rest_samples), probs = 0.05, na.rm = TRUE)
  spectrum$lumQ95       = colQuantiles(as.matrix(rest_samples), probs = 0.95, na.rm = TRUE)
  
  # unattenuated stars
  spectrum$StarsUnAttenWave = out$bestfit$SEDout$StarsUnAtten$wave
  spectrum$StarsUnAttenQ50 = colMedians(as.matrix(StarsUnAtten_samples), na.rm = TRUE)
  spectrum$StarsUnAttenQ16 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.16, na.rm = TRUE)
  spectrum$StarsUnAttenQ84 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.84, na.rm = TRUE)
  spectrum$StarsUnAttenQ05 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.05, na.rm = TRUE)
  spectrum$StarsUnAttenQ95 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.95, na.rm = TRUE)
  
  # attenuated stars
  spectrum$StarsAttenWave = out$bestfit$SEDout$StarsAtten$wave
  spectrum$StarsAttenQ50   = colMedians(as.matrix(StarsAtten_samples), na.rm = TRUE)
  spectrum$StarsAttenQ16   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.16, na.rm = TRUE)
  spectrum$StarsAttenQ84   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.84, na.rm = TRUE)
  spectrum$StarsAttenQ05   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.05, na.rm = TRUE)
  spectrum$StarsAttenQ95   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.95, na.rm = TRUE)
  
  # dust emission
  spectrum$DustWave = out$bestfit$SEDout$DustEmit$wave
  spectrum$DustEmitQ50     = colMedians(as.matrix(DustEmit_samples), na.rm = TRUE)
  spectrum$DustEmitQ16     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.16, na.rm = TRUE)
  spectrum$DustEmitQ84     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.84, na.rm = TRUE)
  spectrum$DustEmitQ05     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.05, na.rm = TRUE)
  spectrum$DustEmitQ95     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.95, na.rm = TRUE)
  
  # AGN
  if(!is.null(out$bestfit$SEDout$AGN)){
    spectrum$AGNWave = out$bestfit$SEDout$AGN$wave
    spectrum$AGNQ50          = colMedians(as.matrix(AGN_samples), na.rm = TRUE)
    spectrum$AGNQ16          = colQuantiles(as.matrix(AGN_samples), probs = 0.16, na.rm = TRUE)
    spectrum$AGNQ84          = colQuantiles(as.matrix(AGN_samples), probs = 0.84, na.rm = TRUE)
    spectrum$AGNQ05          = colQuantiles(as.matrix(AGN_samples), probs = 0.05, na.rm = TRUE)
    spectrum$AGNQ95          = colQuantiles(as.matrix(AGN_samples), probs = 0.95, na.rm = TRUE)
  }else{
    spectrum$AGNWave = 1:10
    spectrum$AGNQ50 = rep(1e-99, 10)
    spectrum$AGNQ16 = rep(1e-99, 10)
    spectrum$AGNQ84 = rep(1e-99, 10)
    spectrum$AGNQ05 = rep(1e-99, 10)
    spectrum$AGNQ95 = rep(1e-99, 10)
  }
  
  # observed frame flux
  spectrum$obswave         = out$bestfit$SEDout$FinalFlux$wave
  spectrum$fluxQ50         = colMedians(as.matrix(obs_samples), na.rm = TRUE)
  spectrum$fluxQ16         = colQuantiles(as.matrix(obs_samples), probs = 0.16, na.rm = TRUE)
  spectrum$fluxQ84         = colQuantiles(as.matrix(obs_samples), probs = 0.84, na.rm = TRUE)
  spectrum$fluxQ05         = colQuantiles(as.matrix(obs_samples), probs = 0.05, na.rm = TRUE)
  spectrum$fluxQ95         = colQuantiles(as.matrix(obs_samples), probs = 0.95, na.rm = TRUE)
  
  # observed frame flux
  spectrum$obswave         = out$bestfit$SEDout$FinalFlux$wave
  spectrum$fluxQ50         = colMedians(as.matrix(obs_samples), na.rm = TRUE)
  spectrum$fluxQ16         = colQuantiles(as.matrix(obs_samples), probs = 0.16, na.rm = TRUE)
  spectrum$fluxQ84         = colQuantiles(as.matrix(obs_samples), probs = 0.84, na.rm = TRUE)
  spectrum$fluxQ05         = colQuantiles(as.matrix(obs_samples), probs = 0.05, na.rm = TRUE)
  spectrum$fluxQ95         = colQuantiles(as.matrix(obs_samples), probs = 0.95, na.rm = TRUE)
  
  # observed frame flux
  spectrum$ageQ50         = colMedians(as.matrix(agevec_samples), na.rm = TRUE)
  spectrum$ageQ16         = colQuantiles(as.matrix(agevec_samples), probs = 0.16, na.rm = TRUE)
  spectrum$ageQ84         = colQuantiles(as.matrix(agevec_samples), probs = 0.84, na.rm = TRUE)
  spectrum$ageQ05         = colQuantiles(as.matrix(agevec_samples), probs = 0.05, na.rm = TRUE)
  spectrum$ageQ95         = colQuantiles(as.matrix(agevec_samples), probs = 0.95, na.rm = TRUE)
  
  spectrum$SFHQ50         = colMedians(as.matrix(sfhvec_samples), na.rm = TRUE)
  spectrum$SFHQ16         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.16, na.rm = TRUE)
  spectrum$SFHQ84         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.84, na.rm = TRUE)
  spectrum$SFHQ05         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.05, na.rm = TRUE)
  spectrum$SFHQ95         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.95, na.rm = TRUE)
  
  spectrum$ZHQ50         = colMedians(as.matrix(Zvec_samples), na.rm = TRUE)
  spectrum$ZHQ16         = colQuantiles(as.matrix(Zvec_samples), probs = 0.16, na.rm = TRUE)
  spectrum$ZHQ84         = colQuantiles(as.matrix(Zvec_samples), probs = 0.84, na.rm = TRUE)
  spectrum$ZHQ05         = colQuantiles(as.matrix(Zvec_samples), probs = 0.05, na.rm = TRUE)
  spectrum$ZHQ95         = colQuantiles(as.matrix(Zvec_samples), probs = 0.95, na.rm = TRUE)
  
  pdf(
    file = paste0(dirname(save_stub), "/", out$galInfo[1,"galName"], "_", fittype, "_plot.pdf"), 
    width = 10, 
    height = 12,
    pointsize = 16
  )
  par(
    mfrow = c(3,1), 
    mar = c(3.5, 3.5, 1.0, 3.5),   # bottom, left, top, right — trimmed top/right, kept left for y-labels
    oma = c(2.5, 2.5, 1.0, 2.5),   # shared bottom margin for a common x-axis label across all 3 panels
    cex = 1.1,                      # dialed back from 2.0 -- try this first
    cex.axis = 1.1,
    cex.lab = 1.1,
    lwd = 2,
    mgp = c(2.2, 0.7, 0)            # pulls axis title/labels/ticks closer to the plot
  )
  obs_xlo = min(Data$flux$cenwave) - 100
  obs_xhi = max(Data$flux$cenwave) + 1000
  obs_ylo = min(pmax(1e-10, Data$flux$flux - 3*Data$flux$fluxerr, na.rm = TRUE)) * 0.5
  obs_yup = max(Data$flux$flux + 3*Data$flux$fluxerr, na.rm = TRUE) * 2.5
  
  plot_spectrum_ribbon(
    spectrum = spectrum, key = "flux", xvar = "obswave", xlim = c(obs_xlo, obs_xhi), ylim = c(obs_ylo, obs_yup), xlab = "Wavelength [micron]", ylab = "Flux [Jy]"
  )
  for(index in 1:length(Data$flux$cenwave)){
    x = Data$flux$cenwave[index]
    
    q50 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ50")])
    q16 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ16")])
    q84 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ84")])
    
    polygon(
      c(x-1e3, x-1e3, x+1e3, x+1e3),
      c(q16, q84, q84, q16),
      col    = adjustcolor(col = "forestgreen", 0.4),
      border = NA
    )
    points(
      x, q50, col = "forestgreen", pch = 16
    )
    
    q50 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ50")])
    q16 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ16")])
    q84 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ84")])
    
    polygon(
      c(x-1e3, x-1e3, x+1e3, x+1e3),
      c(q16, q84, q84, q16),
      col    = adjustcolor(col = "magenta", 0.4),
      border = NA
    )
    points(
      x, q50, col = "magenta", pch = 16
    )
  }
  points(
    Data$flux$cenwave,
    Data$flux$flux,
    col = "black",
    pch = 1
  )
  magerr(
    Data$flux$cenwave,
    Data$flux$flux,
    ylo = Data$flux$fluxerr, 
    col = "black"
  )
  legend(
    x = "bottomright",
    c(
      NA,
      NA,
      NA,
      "Data",
      "Fit",
      "Fit photom",
      "Fit AGN photom"
    ),
    pch = c(NA, NA, NA, 1, NA, 16, 16),
    lty = c(NA, NA, NA, NA, 1, NA, NA),
    lwd = c(NA, NA, NA, NA, 5, NA, NA),
    col = c(NA, NA, NA, "black", "grey", "forestgreen", "magenta"),
    cex = 0.5
  )
  
  plot_spectrum_ribbon(spectrum, key = "StarsUnAtten", xvar = "StarsUnAttenWave", col = "blue", xlim = c(10, 1e7))
  plot_spectrum_ribbon(spectrum, key = "StarsAtten", xvar = "StarsUnAttenWave", col = "cornflowerblue", add = TRUE)
  plot_spectrum_ribbon(spectrum, key = "AGN", xvar = "AGNWave", col = "purple", add = TRUE)
  plot_spectrum_ribbon(spectrum, key = "DustEmit", xvar = "DustWave", col = "brown", add = TRUE)
  plot_spectrum_ribbon(spectrum, key = "lum", xvar = "restWave", col = "black", add = TRUE, alpha = 0.2)
  
  plot_sfh_zh(spectrum, key = "SFH", xvar = "ageQ50", col = "blue", ylab = "SFR [Msun yr^-1]", xlab = "Lookback time since observation [Gyr]", log = "x")
  plot_sfh_zh(spectrum, key = "ZH", xvar = "ageQ50", col = "red", ylab2 = "Z [mass fraction]", add = TRUE, log = "x")
  dev.off()
  
  dir.create(
    dirname(save_stub), 
    recursive = TRUE
  )
  saveRDS(
    spectrum, 
    save_stub
  )
}

compute_fit_sed(paste0(main_dir, "/data/literature_lrds/Cliff_lrdmbb_fit.rds"), LRDBBSKIRTOR, "lrdmbb")
compute_fit_sed(paste0(main_dir, "/data/literature_lrds/Cliff_fritz_fit.rds"), Fritz, "fritz")
compute_fit_sed(paste0(main_dir, "/data/literature_lrds/Cliff_lrdliu_fit.rds"), LRDLIU, "lrdliu")
compute_fit_sed(paste0(main_dir, "/data/literature_lrds/Cliff_sirocco_fit.rds"), SIROCCO, "sirocco")

# compute_fit_sed(paste0(main_dir, "/data/literature_lrds/MoMBH1_lrdmbb_fit.rds"), LRDBBSKIRTOR, "lrdmbb")
# compute_fit_sed(paste0(main_dir, "/data/literature_lrds/MoMBH1_fritz_fit.rds"), Fritz, "fritz")
# compute_fit_sed(paste0(main_dir, "/data/literature_lrds/MoMBH1_lrdliu_fit.rds"), LRDLIU, "lrdliu")
# compute_fit_sed(paste0(main_dir, "/data/literature_lrds/MoMBH1_sirocco_fit.rds"), SIROCCO, "sirocco")
# traceback()