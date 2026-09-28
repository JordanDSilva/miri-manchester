library(data.table)
library(foreach)
library(celestial)
library(matrixStats)
library(magicaxis)
library(doParallel)

set.seed(666)

LAGN2BHAR = function(epsilon){( 3.154e7 * 1e-7 * (epsilon * 3e8^2)**-1 / 2e30)}
epsilon_samples = runif(n = 100000, min = 0.06, max = 0.4)

LAGN_thresh = 1e44

uni_age_fast = approxfun(
  x = seq(0, 100, 0.01),
  y = cosdistUniAgeAtz(z = seq(0, 100, 0.01), ref = "737"),
  rule = 2
)

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

registerDoParallel(cores = 8)
df_combine = data.frame(foreach(ii = 1:dim(prospect_cat)[1], .combine = rbind, .errorhandling = "stop") %dopar% {
  message("Running galID: ", ii)
  if(prospect_cat$FitType[ii] == "Stellar"){
    read_stub = paste0(
      "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy/", "/galID_",ii, ".rds"
    )
  }else{
    read_stub = paste0(
      "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/", "/galID_",ii, ".rds"
    )
  }
  
  if(file.exists(read_stub)){
    out = readRDS(read_stub)
    post = out$post
    
    z = post[,"z"]
    uniage = uni_age_fast(z)
     
    StellarMass = post[, "StellarMass"]
    SFR10 = post[, "SFR10"]
    
    Npost = length(SFR10)
    
    # quench_thresh = 0.2 * (1e9 * uniage) * StellarMass
    
    # quenched_idx = SFR10 <= quench_thresh
    # SFR10[quenched_idx] = quench_thresh[quenched_idx]
    # 
    if(is.null(out$bestfit$SEDout$AGN)){
      BHAR_SFR = rep(0,length(SFR10))
    }else{
      AGNlum = post[,"AGNlum"]
      BHAR_SFR = SFR10 / (AGNlum * LAGN2BHAR(epsilon_samples[sample(1:Npost)]))
      # BHAR_SFR = SFR10 / (AGNlum * LAGN2BHAR(0.1))
      
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
  return(df)
})
stopImplicitCluster()

min_SFR10_values = prospect_cat$StellarMassQ50 * 0.2/(1e9*uni_age_fast(prospect_cat$zQ50))
quenched_idx = which(prospect_cat$SFR10Q50 < min_SFR10_values)
df_combine$lowerlim = FALSE
# df_combine$lowerlimVal = 0
df_combine$lowerlim[quenched_idx] = TRUE
# df_combine$lowerlimVal[quenched_idx] = prospect_cat$AGNlumQ50[quenched_idx] * LAGN2BHAR / min_SFR10_values[quenched_idx]

df_combine[is.na(df_combine)] <- -999
fwrite(
  data.frame(df_combine), 
  "~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_BHAR_SFR_ratio.csv"
)

