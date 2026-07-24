library(data.table)
library(foreach)
library(celestial)
library(matrixStats)
library(magicaxis)
library(doParallel)

LAGN2BHAR = ( 3.154e7 * 1e-7 * (0.1 * 3e8^2)**-1 / 2e30)

uni_age_fast = approxfun(
  x = seq(0, 100, 0.01),
  y = cosdistUniAgeAtz(z = seq(0, 100, 0.01), ref = "737"),
  rule = 2
)

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

registerDoParallel(cores = 8)
df_combine = data.frame(foreach(ii = 1:dim(prospect_cat)[1], .combine = rbind, .errorhandling = "stop") %dopar% {
  message("Running galID: ", ii)
  if(prospect_cat$FitType[ii] == "Stellar"){
    read_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy/", "/galID_",ii, ".rds"
    )
  }else{
    read_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/", "/galID_",ii, ".rds"
    )
  }
  
  if(file.exists(read_stub)){
    out = readRDS(read_stub)
    post = out$post
    
    z = post[,"z"]
    uniage = uni_age_fast(z)
     
    StellarMass = post[, "StellarMass"]
    SFR10 = post[, "SFR10"]
    
    # quench_thresh = 0.2 * (1e9 * uniage) * StellarMass
    
    # quenched_idx = SFR10 <= quench_thresh
    # SFR10[quenched_idx] = quench_thresh[quenched_idx]
    # 
    if(is.null(out$bestfit$SEDout$AGN)){
      BHAR_SFR = 0/SFR10
    }else{
      AGNlum = post[,"AGNlum"]
      BHAR_SFR = (AGNlum * LAGN2BHAR) / SFR10
    }
    
    df = data.frame(
      "BHARSFR10Q50" = quantile(BHAR_SFR, probs = 0.50, na.rm = TRUE),
      "BHARSFR10Q16" = quantile(BHAR_SFR, probs = 0.16, na.rm = TRUE),
      "BHARSFR10Q84" = quantile(BHAR_SFR, probs = 0.84, na.rm = TRUE)
    )
    df$BHARSFR10_lo = df$BHARSFR10Q50 - df$BHARSFR10Q16
    df$BHARSFR10_hi = df$BHARSFR10Q84 - df$BHARSFR10Q50
    
    rm(out)
    gc()
  }else{
    df = data.frame(
      "BHARSFR10Q50" = NA,
      "BHARSFR10Q16" = NA,
      "BHARSFR10Q84" = NA
    )
    df$BHARSFR10_lo = df$BHARSFR10Q50 - df$BHARSFR10Q16
    df$BHARSFR10_hi = df$BHARSFR10Q84 - df$BHARSFR10Q50
  }
  df$ID = prospect_cat$ID[ii]
  return(df)
})
stopImplicitCluster()

min_SFR10_values = prospect_cat$StellarMassQ50 * 0.2/(1e9*uni_age_fast(prospect_cat$zQ50))
quenched_idx = which(prospect_cat$SFR10Q50 < min_SFR10_values)
df_combine$lowerlim = FALSE
df_combine$lowerlimVal = 0
df_combine$lowerlim[quenched_idx] = TRUE
df_combine$lowerlimVal[quenched_idx] = prospect_cat$AGNlumQ50[quenched_idx] * LAGN2BHAR / min_SFR10_values[quenched_idx]


fwrite(
  df_combine, 
  "~/Documents/Epochs/data/prospect/prospect_robust_BHAR_SFR_ratio.csv"
)

# magplot(
#   prospect_cat$StellarMassQ50[],
#   df_combine$BHARSFR10Q50[prospect_cat$SFR10Q50 >= prospect_cat$StellarMassQ50 * 0.2/(1e9*uni_age_fast(prospect_cat$zQ50))],
#   log = "xy",
#   z = prospect_cat$zQ50[prospect_cat$SFR10Q50 >= prospect_cat$StellarMassQ50 * 0.2/(1e9*uni_age_fast(prospect_cat$zQ50))]
# )
