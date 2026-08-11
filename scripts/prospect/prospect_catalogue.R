library(foreach)
library(doParallel)
library(data.table)

set.seed(666)

ncores = 8 ## will set to 1 for SLURM job array 
ncoresPost = 1 ## how many cores for posterior sampling
do_spec = FALSE
do_lrd = TRUE

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
  fittype = "photoZ_bc03_snorm_fritz_norm_eazy"
} else {
  fittype <- args[3]
}

if(do_spec){
  save_dir = paste0(
    main_dir, "data/prospect/", fittype, "_spec/"
  )
  sed_dir = paste0(
    main_dir, "data/prospect_seds/", fittype, "_spect/"
  )
  cat_stub = paste0(main_dir, "data/prospect/", fittype, "_spec_astro_cat.csv")
}else{
  save_dir = paste0(
    main_dir, "data/prospect/", fittype, "/"
  )
  sed_dir = paste0(
    main_dir, "data/prospect_seds/", fittype, "/"
  )
  cat_stub = paste0(main_dir, "data/prospect/", fittype, "_astro_cat.csv")
}

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)

spec_ids = c(
  2193, 725, 781, 2797, 3870, 3123, 5507, 2269, 1903, 3618, 2294, 2362, 3327, 1342, 1869, 1245, 1733, 1810, 196, 4459, 743, 279, 2074, 2591, 2247, 578, 3356, 3012, 3294, 2535,
  491, 1445, 971, 681, 3211, 2020, 4503, 2384, 1779, 3679, 1481, 684, 1456, 2214, 1655, 1715, 4488, 2838, 122, 4542, 945, 2342, 2761, 2758, 826, 261, 3281, 912, 1489, 895, 594,
  1102, 1739, 3380, 372, 5491, 658, 1762, 3014, 746, 600, 1941, 1846, 1095, 629, 1203, 690, 1307, 928, 1186, 30, 1161, 1513, 4950, 919, 1727, 706, 2652, 2357, 1272, 1013, 1116,
  779, 1309, 2864, 1870, 1276, 291, 4500, 383, 1343, 3742, 716, 1735, 1094, 3246, 3377, 790, 3483, 621, 482, 3612, 2052, 108, 1625, 3198, 23, 2531, 1469, 1601, 2536, 5631, 3825,
  1544, 4610, 1806, 2932, 1136, 2827, 1145, 251, 191, 1334, 142, 8, 1225, 2199, 1129, 2058, 915, 4561, 541, 1586, 461, 2712, 4844, 2347, 4479, 1809, 2644, 1007, 585, 3401, 2675,
  2437, 3010, 190, 2683, 1453, 1745, 5219, 2257, 2962, 2843, 2235, 298
)
if(do_spec){
  idx = spec_ids
}else if(do_lrd){
  idx = unlist(
     fread(
       paste0(main_dir, "/data/prospect/LRD_colour_selection_IDs.txt")
     )
   )
}

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
  }
  return(astro)
})
stopImplicitCluster()

df_combine$CATID = 1:dim(df_combine)[1]
fwrite(
  df_combine, 
  cat_stub
)