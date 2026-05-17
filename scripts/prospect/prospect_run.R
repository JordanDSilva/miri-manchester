library(ProSpect)
library(ProSpectUtils)
library(celestial)
library(ParmOff)
library(Highlander)
library(foreach)
library(matrixStats)
library(doParallel)
library(magicaxis)
library(data.table)

set.seed(666)

# main_dir = "/scratch/"

CMANiters = 500
LDNiters = 1000
ncores = 1 
max_slurm_jobs = 1 # should be the maximum of your job array
slurm_paralel_job_array = FALSE

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
  main_dir = "~/Documents/Epochs/"
} else {
  main_dir <- args[2]
}
if (length(args[3]) == 0 | is.na(args[3])) {
  fittype = "photoZ_bc03_snorm_fritz"
} else {
  fittype <- args[3]
}
dir.create(
  paste0(main_dir, "/prospect/", fittype), 
  recursive = TRUE
)

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)
## Frame info, filter info, throughputs loaded 
## prior and constraints loaded
## eazy zPDFs loaded

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

# match_nirspec_idx = coordmatch(
#   coordref = smiles_miri_nirspec[, c("RA", "DEC")],
#   coordcompare = phot_cat$ProFound[, c("RAmax", "Decmax")],
#   rad = 0.3
# )
# nirspec_z = smiles_miri_nirspec$Z_SPEC[match_nirspec_idx$bestmatch$refID]
# idx = match_nirspec_idx$bestmatch$compareID[idx]

## standard BC03 with Fritz AGN model 
do_photoz_bc03_snorm_fritz_fit = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/prospect/", fittype, "/galID_",kk, ".rds" 
  )

  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error

  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc(parm) + 
      zPDF_func(10^parm['z'])
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
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_zPDF,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  startpoint["z"] = log10(zeazy) ## start point in redshift that most galaxies are probably close to z=1 ? 
  
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
    "galInfo" = cbind( phot_cat$ProFound[kk, c("RAmax", "Decmax", "ID")], "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}
## standard BC03 with NO Fritz/AGN model
do_photoz_bc03_snorm_fit = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/prospect/", fittype, "/galID_",kk, ".rds" 
  )
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc(parm) + 
      zPDF_func(10^parm['z'])
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
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7
    ),
    
    speclib = BC03hr,
    Dale = Dale_NormTot,
    filtout = filter_transforms,
    
    SFH=SFHfunc,
    Dale_M2L_func = Dale_M2L_func,
    parm.names = c('z', 'mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal'), 
    mon.names = c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst", paste("flux.",filters_info$filters,sep='')),
    logged = c(
      TRUE,
      TRUE, FALSE, TRUE, FALSE,
      TRUE, TRUE, FALSE, FALSE,
      TRUE
    ),
    intervals = list(
      lo = c(
        -2,
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0, 
        -4
      ), 
      hi = c(
        log10(40), 
        4, 13.46728, 2, 1, 
        1.5, 1, 4, 4, 
        -1.3
      )
    ), 
    fit = 'LD',
    N = dim(flux_data)[1],
    prior = pfunc_zPDF,
    constraints = NULL,
    verbose = FALSE
  )
  
  Data$constraints = function(parm){constraintfunc(parm, Data)}
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["z"] = log10(zeazy) ## start point in redshift that most galaxies are probably close to z=1 ? 

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
    "galInfo" = cbind( phot_cat$ProFound[kk, c("RAmax", "Decmax", "ID")], "LP" = highout$LP, "DIC" = highout$LD_last$DIC1[3])
  )
  return(out)
}

if(fittype == "photoZ_bc03_snorm_fritz"){
  do_fit = do_photoz_bc03_snorm_fritz_fit
}else if (fittype == "photoZ_bc03_snorm"){
  do_fit = do_photoz_bc03_snorm_fit
}

index_parallel_job_array = function(INDEX, ncores, N, start_Nindex, end_index){
  # will find optimal lots of ncores to fit N galaxies starting at start_Nindex:... inside of a 1 to max_index's slurm JOB array
  # Example N = 200000, ncores = 7, max_index = 800
  # At the minimum I will need floor(20000/(800*7)) = 35 lots of 7 indices 
  # JOB = 1 indices will be 1:7, then 5601:5608, then 11201:11207 ... 35 lots of these < 200000
  # 800*7*35 = 196000 so there will be 4000 remainders 
  # floor(4000/7) = 571 so tack on 7 lots of remainders for JOB = 1-571, and the last 3 stragglers on to JOB=572
  # I could fit 1:100000 galaxies N=100000, then fit 100001:200000 but need to set start_Nindex=100001
  
  maxjobs = end_index * ncores
  how_many_groups = floor(N/maxjobs)
  N_remaining = N - maxjobs*how_many_groups
  
  a = INDEX*ncores - ncores 
  b = INDEX*ncores - 1
  
  if(INDEX == 1){ 
    a = 1 
    b = ncores
    c = 0
  }else{
    c = 1
  }
  job_vec_init = (a:b) + c
  
  job_vec = c()
  for(k in 1:how_many_groups){
    job_vec = c(
      job_vec, job_vec_init + (k-1)*maxjobs
    )
  }
  
  if(N_remaining > 0 ){
    remaining_indices = (maxjobs*how_many_groups + 1):N
    job_vec = as.numeric(na.omit( c(job_vec, remaining_indices[job_vec_init]) )) + (start_Nindex-1)
  }else{
    job_vec = as.numeric(na.omit( c(job_vec) )) + (start_Nindex-1)
  }
  
  return( job_vec ) 
}

if(slurm_paralel_job_array){
  idx_idx = idx 
  
  idx = index_parallel_job_array(
    INDEX = idx_idx, 
    ncores = ncores, 
    N = dim(df_flux)[1], 
    start_Nindex = 1,
    end_index = max_slurm_jobs
  )
}

registerDoParallel(cores = ncores)
foreach(kk = idx) %dopar% {
  
  fit = do_fit(kk)

  summ = suppressMessages(
    ProSpectUtils::calc_astro(
      fit$bestfit,
      fit$highout,
      newDust = TRUE,
      specz = FALSE,
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
  return(NULL)
}
stopImplicitCluster()

if(tail(idx,1) == dim(df_flux)[1]){
  df_combine = foreach(ii = 1:tail(idx,1), .combine = rbind) %do% {
    message("Running galID: ", ii)
    read_stub = paste0(
      main_dir, "/prospect/", fittype, "/galID_",ii, ".rds"
    )
    
    out = readRDS(read_stub)
    astro = out$astro
    astro['DIC'] = out$galInfo$DIC
    
    rm(out)
    gc()
    return(astro)
  }
  fwrite(
    df_combine, paste0(main_dir, "/prospect/", fittype, "_astro_cat.csv")
  )
}

