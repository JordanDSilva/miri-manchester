library(ProSpect)
library(ProSpectUtils)
library(celestial)
library(ParmOff)
library(Highlander)
library(nloptr)
library(foreach)
library(matrixStats)
library(doParallel)
library(magicaxis)
library(data.table)

set.seed(666)

CMANiters = 500
LDNiters = 1000
ncores = 1 ## will set to 1 for SLURM job array 
ncoresPost = 1 ## how many cores for posterior sampling
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

index_parallel_job_array = function(INDEX, ncores){
  a = INDEX*ncores - ncores
  b = INDEX*ncores - 1 
  
  if(INDEX == 1){ 
    a = 1 
    b = ncores
    c = 0
  }else{
    c = 1
  }
  return( (a:b) + c ) 
}

if(slurm_paralel_job_array){
  idx_idx = idx 
  idx = index_parallel_job_array(idx_idx, ncores)
}

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)
## Frame info, filter info, throughputs loaded 
## prior and constraints loaded
## eazy zPDFs loaded

dir.create(
  paste0(main_dir, "/prospect/", fittype), 
  recursive = TRUE
)

# match_nirspec_idx = coordmatch(
#   coordref = smiles_miri_nirspec[, c("RA", "DEC")],
#   coordcompare = phot_cat$ProFound[, c("RAmax", "Decmax")],
#   rad = 0.3
# )
# nirspec_z = smiles_miri_nirspec$Z_SPEC[match_nirspec_idx$bestmatch$refID]
# idx = match_nirspec_idx$bestmatch$compareID

start = Sys.time()
registerDoParallel(cores = ncores)
df_astro = foreach(ii = idx, .combine = rbind) %dopar% {
  message("Running galID: ", ii)
  read_stub = paste0(
    main_dir, "/prospect/", fittype, "/galID_",ii, ".rds"
  )
  save_stub = paste0(
    main_dir, "/prospect/", fittype, "/astro_galID_",ii, ".rds"
  )
  
  out = readRDS(read_stub)
  bestfit = out$bestfit
  
  ## put all of the heavy stuff back to make the catalogue
  bestfit$Data$speclib = BC03hr
  if(fittype == "photoZ_bc03_snorm_fritz"){
    bestfit$Data$AGN = Fritz
  }
  bestfit$Data$Dale = Dale_NormTot
  bestfit$Data$constraints = function(parm){constraintfunc(parm, bestfit$Data)}
  bestfit$Data$prior = pfunc
  
  summ = suppressMessages(
    ProSpectUtils::calc_astro(
      bestfit,
      out$highout,
      newDust = TRUE,
      specz = FALSE,
      cores = ncoresPost,
    )
  )
  
  out$astro = summ$summary
  out$post = summ$Posterior
  
  saveRDS(
    out, 
    save_stub
  )
  
  rm(out)
  rm(bestfit)
  gc()
  return(NULL)
}
stopImplicitCluster()

if(tail(idx,1) == dim(phot_cat$ProFound)[1]){
  registerDoParallel(cores = ncores)
  df_combine = foreach(ii = 1:tail(idx,1), .combine = rbind) %dopar% {
    message("Running galID: ", ii)
    read_stub = paste0(
      main_dir, "/prospect/", fittype, "/astro_galID_",ii, ".rds"
    )
    
    out = readRDS(read_stub)
    astro = out$astro
    astro$DIC = out$galInfo$DIC
    
    rm(out)
    gc()
    return(astro)
  }
  stopImplicitCluster()
  fwrite(
    df_combine, paste0(main_dir, "/prospect/", fittype, "_astro_cat.csv")
  )
}
print(Sys.time() - start)
