library(foreach)
library(doParallel)
library(data.table)

set.seed(666)

ncores = 8 ## will set to 1 for SLURM job array 
ncoresPost = 1 ## how many cores for posterior sampling
do_spec = FALSE
do_lrd = FALSE

args <- commandArgs(trailingOnly = TRUE)
if (length(args[1]) == 0 | is.na(args[1])) {
  idx <- 1
} else if (grepl(":", args[1])) {
  bounds <- as.integer(strsplit(args[1], ":")[[1]])
  idx <- bounds[1]:bounds[2]
} else if (grepl(",", args[1])) {
  idx <- as.integer(strsplit(args[1], ",")[[1]])
} else {
  idx <- as.integer(args[1])
}

if (length(args[2]) == 0 | is.na(args[2])) {
  main_dir = "~/Documents/JWST-MIRI-AGN/"
} else {
  main_dir <- args[2]
}
if (length(args[3]) == 0 | is.na(args[3])) {
  fittype = "photoZ_bc03_snorm_fritz_norm_eazy"
} else {
  fittype <- args[3]
}

# if(do_lrd){
#   save_dir = paste0(
#     main_dir, "data/prospect/", fittype, "_lrd/"
#   )
#   sed_dir = paste0(
#     main_dir, "data/prospect_seds/", fittype, "_lrd/"
#   )
#   cat_stub = paste0(main_dir, "data/prospect/", fittype, "_lrd_astro_cat.csv")
# }else{
#   save_dir = paste0(
#     main_dir, "data/prospect/", fittype, "/"
#   )
#   sed_dir = paste0(
#     main_dir, "data/prospect_seds/", fittype, "/"
#   )
#   cat_stub = paste0(main_dir, "data/prospect/", fittype, "_astro_cat.csv")
# }
# 
# load(
#   paste0(main_dir, "/data/prospect_data_for_fitting.rda")
# )
# 
# if(do_lrd){
#   idx = unlist(
#      fread(
#        paste0(main_dir, "/data/prospect/LRD_colour_selection_IDs.txt")
#      )
#    )
# }

template_cols = names(
  readRDS(
    list.files(save_dir, full.names = TRUE)[1]
  )$astro
)
template_cols_sed = names(
  readRDS(
    list.files(sed_dir, full.names = TRUE)[1]
  )$phot
)

print(idx)

registerDoParallel(cores = ncores)
df_combine = data.frame(foreach(ii = idx, .combine = rbind) %dopar% {
  message("Running galID: ", ii)
  read_stub = paste0(
    save_dir, "/galID_",ii, ".rds"
  )
  read_sed_stub = paste0(
    sed_dir, "/galID_",ii, ".rds"
  )
  
  if(file.exists(read_stub)){
    out = readRDS(read_stub)
    sed = readRDS(read_sed_stub)
    
    astro = out$astro
    astro = unlist( c( astro, sed$phot[!grepl("ID", names(sed$phot))] ) )
    
    
    # astro$AGN_20m_fracQ50 = sum( c(0, diff(sed$AGNWave[sed$AGNWave >= 5e4 & sed$AGNWave <= 20e4])) * sed$AGNQ50[sed$AGNWave >= 5e4 & sed$AGNWave <= 20e4] ) / sum( c(0, diff(sed$restWave[sed$restWave >= 5e4 & sed$restWave <= 20e4])) * sed$lumQ50[sed$restWave >= 5e4 & sed$restWave <= 20e4] )
    # # astro$AGN_20m_fracQ50 = sum( c(0, diff(sed$AGNWave))[sed$AGNWave <= 20e4] * sed$AGNQ50[sed$AGNWave <= 20e4] ) / sum( c(0, diff(sed$restWave))[sed$restWave <= 20e4] * sed$lumQ50[sed$restWave <= 20e4] )
    # # astro$AGN_20m_fracQ50 = sum( c(0, diff(sed$AGNWave))[sed$AGNWave <= 20e4] * sed$AGNQ50[sed$AGNWave <= 20e4] ) / sum( c(0, diff(sed$restWave))[sed$restWave <= 20e4] * sed$lumQ50[sed$restWave <= 20e4] )
    
    UV_ratio_post = 10^(-0.4 * (out$post[, 'UV1500AGN'] - out$post[, 'UV1500'] ) )
    astro$UV1500_AGN_ratioQ50 = median(UV_ratio_post, na.rm = TRUE)
    astro$UV1500_AGN_ratioQ16 = quantile(UV_ratio_post, probs = 0.16, na.rm = TRUE)
    astro$UV1500_AGN_ratioQ84 = quantile(UV_ratio_post, probs = 0.84, na.rm = TRUE)
    
    astro$ID = out$galInfo$ID
    astro$RA = out$galInfo$RAmax
    astro$Dec = out$galInfo$Decmax
    astro$DIC = out$galInfo$DIC
    astro$NfinalMCMC = dim(out$highout$LD_last$Posterior1)[1]
    astro$selected = 1L
    
    rm(out)
    gc()
  }
  else{
    astro = setNames(as.data.frame(matrix(NA, nrow = 1, ncol = length(c(template_cols, template_cols_sed)))), c(template_cols, template_cols_sed))
    astro$ID = ii
    astro$RA = phot_cat$ProFound$RAmax[ii]
    astro$Dec = phot_cat$ProFound$Decmax[ii]
    astro$DIC = NA
    astro$NfinalMCMC = NA
    astro$selected = 0L
    
    astro$UV1500_AGN_ratioQ50 = NA
    astro$UV1500_AGN_ratioQ16 = NA
    astro$UV1500_AGN_ratioQ84 = NA
  }
  return(astro)
})
stopImplicitCluster()

df_combine$CATID = 1:dim(df_combine)[1]
fwrite(
  df_combine, 
  cat_stub
)
