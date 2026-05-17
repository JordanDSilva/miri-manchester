library(ProSpect)
library(ProSpectUtils)
library(Highlander)
library(celestial)
library(magicaxis)
library(data.table)
library(pracma)
library(foreach)
library(doParallel)
library(matrixStats)
## Fit the mock photometry with the eazyPDFs

CMANiters = 100
LDNiters = 500

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles"

load(
  paste0(
    main_dir, "/prospect_data_for_fitting.rda"
  )
)
zPDF = data.frame(
  fread(
    paste0(main_dir, "/", fname, "_mock_photometry_eazy_zPDF.csv"))
)
photometry = data.frame(
  fread(
    paste0(main_dir, "/", fname, "_mock_photometry.csv")
  )
)
astro = data.frame(
  fread(
    paste0(main_dir, "/", fname, "_mock_astro_cat.csv")
  )
)

df_flux = photometry[, grepl("flux$", names(photometry))]
df_flux_err = photometry[, grepl("flux_err$", names(photometry))]

args <- commandArgs(trailingOnly = TRUE)
if (length(args[1]) == 0 | is.na(args[1])) {
  idx <- 1
} else if (grepl(":", args[1])) {
  bounds <- as.integer(strsplit(args[1], ":")[[1]])
  idx <- bounds[1]:bounds[2]
} else {
  idx <- as.integer(args[1])
}

if (length(args[2]) == 0 | is.na(args[2])) {
  main_dir = "~/Documents/Epochs/source_inject/"
} else {
  main_dir <- args[2]
}
# if (length(args[3]) == 0 | is.na(args[3])) {
#   fittype = "photoZ_bc03_snorm_fritz"
# } else {
#   fittype <- args[3]
# }

do_photoz_bc03_snorm_fritz_fit = function(kk){
  message("Running gal: ", kk)
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) 
  filtout = filter_transforms
  
  zPDF_load = zPDF[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = zPDF$zgrid,
    y = log(zPDF_load),
    rule = 2
  )
  
  zeazy = zPDF$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc(parm) + zPDF_func(10^parm['z']) ## check z logged
    return(prior)
  }
  
  Data = list(
    flux = flux_data,
    mode = "photom",
    arglist = list(
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      
      emission = TRUE,
      photoz = TRUE,
      IGMfunc = "Inoue14",
      
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
    
    SFH=SFHfunc,
    AGN = Fritz,
    Dale_M2L_func = Dale_M2L_func,
    parm.names = c('z', 'mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta'), 
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst"),
    logged = c(
      TRUE,
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE,
      TRUE, FALSE, TRUE
    ),
    intervals = list(
      lo = c(
        -2,
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0, 
        -4, 
        35, 0.001, -1
      ), 
      hi = c(
        log10(40), 
        4, 13.46728, 2, 1, 
        1.5, 1, 4, 4, 
        -1.3, 
        49, 89.990, 1
      )
    ), 
    fit = 'optim',
    N = dim(flux_data)[1],
    prior = pfunc_zPDF,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  startpoint["z"] = log10(zeazy) # log10(1.5) ## start point in redshift that most galaxies are probably close to z=1 ? 
  
  Data$fit = "LD"
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
  
  astro = calc_astro(
    bestfit,
    highout
  )
  gc()
  return(astro)
}

# if(fittype == "photoZ_bc03_snorm_fritz"){
#   do_fit = do_photoz_bc03_snorm_fritz_fit
# }else if (fittype == "photoZ_bc03_snorm"){
#   do_fit = do_photoz_bc03_snorm_fit
# }

registerDoParallel(cores = ncores)
foreach(kk = idx) %dopar% {
  
  prospect_stub = paste0(
    main_dir, "/prospect/", "/galID_", kk, ".rds"
  )
  
  fit = do_photoz_bc03_snorm_fritz_fit(kk)
 
  saveRDS(
    fit,
    prospect_stub
  )
  return(NULL)
}
stopImplicitCluster()

