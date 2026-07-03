library(ProSpect)
library(matrixStats)
library(foreach)
library(doParallel)

## Compute AGNfractions in filters?

ncores = 1

args <- commandArgs(trailingOnly = TRUE)

if (length(args[1]) == 0 | is.na(args[1])) {
  main_dir = "~/Documents/Epochs/"
} else {
  main_dir <- args[1]
}

if (length(args[2]) == 0 | is.na(args[2])) {
  fittype = "photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth"
} else {
  fittype <- args[2]
}

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)

Ndim = dim(phot_cat$ProFound)[1]

results_list <- vector("list", Ndim)
for(i in 1:Ndim){
  message(i)
  
  save_stub = paste0(
    main_dir, "/data/prospect/", fittype, "/galID_", i, ".rds"
  )
  if(file.exists(save_stub)){
    out = readRDS(save_stub)
    
    post = out$post
    
    Data = out$bestfit$Data
    Data$speclib = BC03hr
    
    if(grepl("fritz", fittype)){
      Data$AGN = Fritz
    }
    
    Data$Dale = Dale_NormTot
    
    registerDoParallel(cores = ncores)
    phot_samples = foreach(j = 1:dim(post)[1], .combine = rbind) %dopar% {
      
      post_sample = out$highout$LD_last$Posterior1[j, names(out$bestfit$parm)]
      out_sample = ProSpectSEDlike(
        post_sample,
        Data
      )
      
      if(!is.null(Data$AGN)){
        photom_AGN_source = photom_lum(
          wave = out_sample$SEDout$AGN$wave,
          lum = out_sample$SEDout$AGN$lum,
          filters = filter_transforms, 
          outtype = "Jy", 
          z = 10^out_sample$parm[1], 
          ref = "737"
        )
        AGN_Dale_frac = sum( out_sample$SEDout$AGN$lum[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4] * c(0, diff(out_sample$SEDout$AGN$wave[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4])) ) / sum( out_sample$SEDout$FinalLum$lum[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4] * c(0, diff(out_sample$SEDout$FinalLum$wave[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4])) )
      }else{
        photom_AGN_source = rep(0, length(filter_transforms))
        AGN_Dale_frac = 0
      }
      
      photom_all_source = photom_lum(
        wave = out_sample$SEDout$FinalLum$wave,
        lum = out_sample$SEDout$FinalLum$lum,
        filters = filter_transforms, 
        outtype = "Jy", 
        z = 10^out_sample$parm[1], 
        ref = "737"
      )
      
      photom_all = photom_all_source * Inoue14_IGM(wave = filters_info$cenwave, z = 10^out_sample$parm[1], Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
      photom_AGN = photom_AGN_source *  Inoue14_IGM(wave = filters_info$cenwave, z = 10^out_sample$parm[1], Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
      
      names(photom_AGN) = paste0(filters_info$filter_names, "_fit_AGN")
      names(photom_all) = paste0(filters_info$filter_names, "_fit")
      names(photom_AGN_source) = paste0(filters_info$filter_names, "_source_fit_AGN")
      names(photom_all_source) = paste0(filters_info$filter_names, "_source_fit")
      
      ff = c(
        photom_all,
        photom_AGN,
        photom_all_source,
        photom_AGN_source,
        "AGN_Dale_frac" = AGN_Dale_frac
      )
    }
    stopImplicitCluster()
    
    med_ = colMedians(as.matrix(phot_samples))
    Q16_ = colQuantiles(as.matrix(phot_samples), probs = 0.16)
    Q84_ = colQuantiles(as.matrix(phot_samples), probs = 0.84)
    Q05_ = colQuantiles(as.matrix(phot_samples), probs = 0.05)
    Q95_ = colQuantiles(as.matrix(phot_samples), probs = 0.95)
    
    names(med_) = paste0(names(med_), "_Q50")
    names(Q16_) = paste0(names(Q16_), "_Q16")
    names(Q84_) = paste0(names(Q84_), "_Q84")
    names(Q05_) = paste0(names(Q05_), "_Q05")
    names(Q95_) = paste0(names(Q95_), "_Q95")
    
    row_ <- c(ID = as.integer(i), med_, Q16_, Q84_, Q05_, Q95_)
    
    results_list[[i]] <- row_
  }else{
    results_list[[i]] <- c(
      "ID" = as.integer(i),
      rep(NA, length(filter_transforms) * 2 * 2 * 5 + 5)
    )
  }
}
results_df <- as.data.frame(do.call(rbind, results_list))

write.csv(
  results_df, 
  paste0(main_dir, "/data/prospect/", fittype, "_fit_photom.csv"), 
  row.names = FALSE
)