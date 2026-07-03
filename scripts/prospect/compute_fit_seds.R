library(ProSpect)
library(matrixStats)
library(foreach)
library(doParallel)

ncores = 1 
max_slurm_jobs = 84 # should be the maximum of your job array
slurm_paralel_job_array = TRUE

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
  fittype = "photoZ_bc03_snorm_norm_eazy_no_tau_birth"
} else {
  fittype <- args[3]
}

prospect_dir = paste0(
  main_dir, "/data/prospect/", fittype, "/" 
)
save_dir = paste0(
  main_dir, "/data/prospect_seds/", fittype, "/" 
)
dir.create(
  save_dir, recursive = TRUE
)

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)

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
  fit_idx = index_parallel_job_array(
    INDEX = idx, 
    ncores = ncores, 
    N = dim(phot_cat$Phot32)[1], 
    start_Nindex = 1,
    end_index = max_slurm_jobs
  )
}else{
  fit_idx = idx
}

registerDoParallel(cores = ncores)
for(i in fit_idx){
  message(i)
  
  file_stub = paste0(
    prospect_dir, "galID_", i, ".rds"
  )
  
  if(file.exists(file_stub)){
    save_stub = paste0(
      save_dir, "galID_", i, ".rds"
    )
    
    out = readRDS(file_stub)
    
    post = out$post
    
    Data = out$bestfit$Data
    Data$speclib = BC03hr
    
    if(!is.null(out$bestfit$SEDout$AGN)){
      Data$AGN = Fritz
    }else{
      Data$AGN = NULL
    }
    Data$Dale = Dale_NormTot
    
    lum_samples = foreach(j = 1:dim(post)[1]) %dopar% {
      
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
      list(
        "rest" = out_sample$SEDout$FinalLum$lum,
        "StarsUnAtten" = out_sample$SEDout$StarsUnAtten$lum,
        "StarsAtten" = out_sample$SEDout$StarsAtten$lum,
        "DustEmit" = out_sample$SEDout$DustEmit$lum,
        "AGN" = AGN,
        "obs" = out_sample$SEDout$FinalFlux$flux
      )
    }
    
    rest_samples         = do.call(rbind, lapply(lum_samples, function(x) x$rest))
    StarsUnAtten_samples = do.call(rbind, lapply(lum_samples, function(x) x$StarsUnAtten))
    StarsAtten_samples   = do.call(rbind, lapply(lum_samples, function(x) x$StarsAtten))
    DustEmit_samples     = do.call(rbind, lapply(lum_samples, function(x) x$DustEmit))
    AGN_samples          = do.call(rbind, lapply(lum_samples, function(x) x$AGN))
    obs_samples          = do.call(rbind, lapply(lum_samples, function(x) x$obs))
    
    spectrum = list()
    
    # rest frame total
    spectrum$restWave     = out$bestfit$SEDout$FinalLum$wave
    spectrum$lumQ50       = colMedians(as.matrix(rest_samples))
    spectrum$lumQ16       = colQuantiles(as.matrix(rest_samples), probs = 0.16)
    spectrum$lumQ84       = colQuantiles(as.matrix(rest_samples), probs = 0.84)
    spectrum$lumQ05       = colQuantiles(as.matrix(rest_samples), probs = 0.05)
    spectrum$lumQ95       = colQuantiles(as.matrix(rest_samples), probs = 0.95)
    
    # unattenuated stars
    spectrum$StarsUnAttenWave = out$bestfit$SEDout$StarsUnAtten$wave
    spectrum$StarsUnAttenQ50 = colMedians(as.matrix(StarsUnAtten_samples))
    spectrum$StarsUnAttenQ16 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.16)
    spectrum$StarsUnAttenQ84 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.84)
    spectrum$StarsUnAttenQ05 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.05)
    spectrum$StarsUnAttenQ95 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.95)
    
    # attenuated stars
    spectrum$StarsAttenWave = out$bestfit$SEDout$StarsAtten$wave
    spectrum$StarsAttenQ50   = colMedians(as.matrix(StarsAtten_samples))
    spectrum$StarsAttenQ16   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.16)
    spectrum$StarsAttenQ84   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.84)
    spectrum$StarsAttenQ05   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.05)
    spectrum$StarsAttenQ95   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.95)
    
    # dust emission
    spectrum$DustWave = out$bestfit$SEDout$DustEmit$wave
    spectrum$DustEmitQ50     = colMedians(as.matrix(DustEmit_samples))
    spectrum$DustEmitQ16     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.16)
    spectrum$DustEmitQ84     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.84)
    spectrum$DustEmitQ05     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.05)
    spectrum$DustEmitQ95     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.95)
    
    # AGN
    if(!is.null(out$bestfit$SEDout$AGN)){
      spectrum$AGNWave = out$bestfit$SEDout$AGN$wave
      spectrum$AGNQ50          = colMedians(as.matrix(AGN_samples))
      spectrum$AGNQ16          = colQuantiles(as.matrix(AGN_samples), probs = 0.16)
      spectrum$AGNQ84          = colQuantiles(as.matrix(AGN_samples), probs = 0.84)
      spectrum$AGNQ05          = colQuantiles(as.matrix(AGN_samples), probs = 0.05)
      spectrum$AGNQ95          = colQuantiles(as.matrix(AGN_samples), probs = 0.95)
    }

    # observed frame flux
    spectrum$obswave         = out$bestfit$SEDout$FinalFlux$wave
    spectrum$fluxQ50         = colMedians(as.matrix(obs_samples))
    spectrum$fluxQ16         = colQuantiles(as.matrix(obs_samples), probs = 0.16)
    spectrum$fluxQ84         = colQuantiles(as.matrix(obs_samples), probs = 0.84)
    spectrum$fluxQ05         = colQuantiles(as.matrix(obs_samples), probs = 0.05)
    spectrum$fluxQ95         = colQuantiles(as.matrix(obs_samples), probs = 0.95)
    
    dir.create(
      dirname(save_stub), 
      recursive = TRUE
    )
    saveRDS(
      spectrum, 
      save_stub
    )
  }
}
stopImplicitCluster()