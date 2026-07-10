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
max_slurm_jobs = 84 # should be the maximum of your job array
slurm_paralel_job_array = TRUE
overwrite = FALSE

do_spec = FALSE
do_robust = TRUE
redo_fits = FALSE

min_pc_err = 0.1

args <- commandArgs(trailingOnly = TRUE)
if (length(args[1]) == 0 | is.na(args[1])) {
  idx <- 1
} else if (grepl(":", args[1])) {
  bounds <- as.integer(strsplit(args[1], ":")[[1]])
  idx <- bounds[1]:bounds[2]
} else if (grepl(",", args[1])) {
  idx <- as.integer(strsplit(args[1], ",")[[1]])
}else {
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

if(do_spec){
  save_dir = paste0(
    main_dir, "/data/prospect/", fittype, "_spec/" 
  )
}else{
  save_dir = paste0(
    main_dir, "/data/prospect/", fittype, "/" 
  )
}

dir.create(
  save_dir, recursive = TRUE
)

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)
## Frame info, filter info, throughputs loaded 
## prior and constraints loaded
## eazy zPDFs loaded

robust_ids = which(
  phot_cat$ProFound$flag == 0
)
redo_ids = unlist(
 fread(
  paste0(main_dir, "/data/prospect/poor_fits_IDs.txt")
))

spec_ids = c(
  2193, 725, 781, 2797, 3870, 3123, 5507, 2269, 1903, 3618, 2294, 2362, 3327, 1342, 1869, 1245, 1733, 1810, 196, 4459, 743, 279, 2074, 2591, 2247, 578, 3356, 3012, 3294, 2535,
  491, 1445, 971, 681, 3211, 2020, 4503, 2384, 1779, 3679, 1481, 684, 1456, 2214, 1655, 1715, 4488, 2838, 122, 4542, 945, 2342, 2761, 2758, 826, 261, 3281, 912, 1489, 895, 594,
  1102, 1739, 3380, 372, 5491, 658, 1762, 3014, 746, 600, 1941, 1846, 1095, 629, 1203, 690, 1307, 928, 1186, 30, 1161, 1513, 4950, 919, 1727, 706, 2652, 2357, 1272, 1013, 1116,
  779, 1309, 2864, 1870, 1276, 291, 4500, 383, 1343, 3742, 716, 1735, 1094, 3246, 3377, 790, 3483, 621, 482, 3612, 2052, 108, 1625, 3198, 23, 2531, 1469, 1601, 2536, 5631, 3825,
  1544, 4610, 1806, 2932, 1136, 2827, 1145, 251, 191, 1334, 142, 8, 1225, 2199, 1129, 2058, 915, 4561, 541, 1586, 461, 2712, 4844, 2347, 4479, 1809, 2644, 1007, 585, 3401, 2675,
  2437, 3010, 190, 2683, 1453, 1745, 5219, 2257, 2962, 2843, 2235, 298
)

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

## standard BC03 with Fritz AGN model, normal prior on eazy zbest, no tau birth prior
do_photoz_bc03_snorm_fritz_norm_eazy = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    save_dir, "/galID_",kk, ".rds" 
  )
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (min_pc_err * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs_normals[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs_normals$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc_no_tau_birth(parm) + 
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
## Standard Bc03 with No Fritz, normal prior on eazy zbest, no tau birth prior
do_photoz_bc03_snorm_norm_eazy = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    save_dir, "/galID_",kk, ".rds" 
  )
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (min_pc_err * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs_normals[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs_normals$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc_no_tau_birth(parm) + 
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
## Standard Bc03 with SKIRTOR, normal prior on eazy zbest, no tau birth prior
do_photoz_bc03_snorm_skirtor_norm_eazy = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    save_dir, "/galID_",kk, ".rds" 
  )
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (min_pc_err * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs_normals[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs_normals$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
  pfunc_zPDF = function(parm){
    prior = pfunc_no_tau_birth(parm) + 
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
      AGNq = 1.0,
      AGNp = 1.0,
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
    AGN = SKIRTOR,
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

if(fittype == "photoZ_bc03_snorm_fritz_norm_eazy"){
  do_fit = do_photoz_bc03_snorm_fritz_fit
}else if (fittype == "photoZ_bc03_snorm_norm_eazy"){
  do_fit = do_photoz_bc03_snorm_fit
}else if (fittype == "photoZ_bc03_snorm_skirtor_norm_eazy"){
  do_fit = do_photoz_bc03_snorm_skirtor_norm_eazy
}else{
  stop("Incorrect fit function")
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
  if(do_robust){
    job_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = length(robust_ids), 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
    fit_idx = robust_ids[job_idx]
  }else if(do_spec){
    job_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = length(spec_ids), 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
    fit_idx = spec_ids[job_idx]
  }else if(redo_fits){
    job_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = length(redo_ids), 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
    fit_idx = redo_ids[job_idx]
  }else{
    fit_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = dim(df_flux)[1], 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
  }
}else{
  if(do_robust){
    fit_idx = robust_ids[idx]
  }else if(do_spec){
    fit_idx = spec_ids[idx]
  }else if(redo_fits){
    fit_idx = redo_ids[idx]
  }else{
    fit_idx = idx
  }
}

registerDoParallel(cores = ncores)
foreach(kk = fit_idx) %dopar% {
  
  save_stub = paste0(
    save_dir, "/galID_",kk, ".rds" 
  )
  
  if(file.exists(save_stub) & !overwrite){
    message("file exists")
  }else{
    fit = do_fit(kk)
    
    # fit$bestfit$Data$arglist$escape_frac = 1
    # new = ProSpectSEDlike(
    #   fit$bestfit$parm,
    #   fit$bestfit$Data
    # )
    
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
  }
  return(NULL)
}
stopImplicitCluster()
