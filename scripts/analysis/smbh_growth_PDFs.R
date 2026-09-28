library(celestial)
library(ProSpect)
library(data.table)
library(foreach)
library(stringr)
library(magicaxis)
library(Highlander)
library(ks)
library(MASS)
library(fields)
library(lamW)
library(rhdf5)
library(hyper.fit)
library(rgl)
library(matrixStats)
library(LaplacesDemon)

Msun = 1.98840987e+30
speedOfLight = 299792458
eddTime = ((1.26e38 * 1e-7) / (Msun * speedOfLight^2))^-1 / (60 * 60 * 24 * 365.2425) # in years

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

uni_age_fast = approxfun(
  x = seq(0, 100, 0.01),
  y = cosdistUniAgeAtz(z = seq(0, 100, 0.01), ref = "737"),
  rule = 2
)
Lagn_func = function(t, M0, lambda, epsilon, t0){
  lambda * 1.26e38 * M0 * exp(
    (t - t0)/(eddTime/1e9) * (lambda * (1 - epsilon))/epsilon
  )
}
SMBH_func = function(t, M0, lambda, epsilon, t0){
  M0 * exp(
    (t - t0)/(eddTime/1e9) * (lambda * (1 - epsilon))/epsilon
  )
}

t0 = cosdistUniAgeAtz(z = 30, ref = "737")

# do distributions of things
N = 100000
M0 = runif(n = N, min = 2, max = 5)
epsilon = runif(n = N, min = 0.06, max = 0.4)

results_fit = list()
for(kk in 1:nrow(prospect_cat)){
  if(prospect_cat$FitType[kk] == "StellarAGN"){
    
    # read in posterior
    read_gal = readRDS(
      paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", kk, ".rds")
    )
    
    uni_age = cosdistUniAgeAtz(z = prospect_cat$zQ50[kk], ref = "737")
    
    Lagn_post = log10( read_gal$post[,c("AGNlum")] )
    lagn_dens = density(Lagn_post)
    
    cdf_vals <- cumsum(lagn_dens$y) * diff(lagn_dens$x)[1]
    cdf_vals <- cdf_vals / max(cdf_vals)   # normalize to end exactly at 1
    
    # Invert: given u, find corresponding x
    inv_cdf <- approxfun(cdf_vals, lagn_dens$x, rule = 2)

    u <- runif(N)
    Lagn_sample <- 10^inv_cdf(u)
    
    # plot(density(Lagn_post))
    # plot(density(Lagn_sample))
    
    arg_eps = (epsilon / (1e9*(uni_age - t0) * (1 - epsilon)))
    arg1 = (eddTime * 1.26e38 * M0)^-1 * Lagn_sample * arg_eps^-1 
    lambda = eddTime * arg_eps * lambertW0(arg1)
    
    SMBH_post = SMBH_func(
      M0 = M0,
      t0 = t0,
      epsilon = epsilon,
      t = uni_age,
      lambda = lambda
    )

    SMBH_post = pmax(M0, SMBH_post)
    
    df = c(
      "Q50" = quantile(lambda, probs = 0.5, na.rm = TRUE),
      "Q16" = quantile(lambda, probs = 0.16, na.rm = TRUE),
      "Q84" = quantile(lambda, probs = 0.84, na.rm = TRUE),
      "lo" = quantile(lambda, probs = 0.5, na.rm = TRUE) - 
        quantile(lambda, probs = 0.16, na.rm = TRUE),
      "hi" = quantile(lambda, probs = 0.84, na.rm = TRUE) - 
        quantile(lambda, probs = 0.5, na.rm = TRUE),
      "SMBH_Q50" = quantile(SMBH_post, probs = 0.5, na.rm = TRUE),
      "SMBH_Q16" = quantile(SMBH_post, probs = 0.16, na.rm = TRUE),
      "SMBH_Q84" = quantile(SMBH_post, probs = 0.84, na.rm = TRUE),
      "SMBH_lo" = quantile(SMBH_post, probs = 0.5, na.rm = TRUE) - quantile(SMBH_post, probs = 0.16, na.rm = TRUE),
      "SMBH_hi" = quantile(SMBH_post, probs = 0.84, na.rm = TRUE) - quantile(SMBH_post, probs = 0.5, na.rm = TRUE)
    )
    
    names(df) = c(
      "ledd_Q50",
      "ledd_Q16",
      "ledd_Q84",
      "ledd_lo",
      "ledd_hi",
      "smbh_Q50",
      "smbh_Q16",
      "smbh_Q84",
      "smbh_lo",
      "smbh_hi"
    )
    
  }else{
    df = c(
      "ledd_Q50" = -999,
      "ledd_Q16" = -999,
      "ledd_Q84" = -999,
      "ledd_lo"  = 0, 
      "ledd_hi"  = 0, 
      "smbh_Q50" = -999,
      "smbh_Q16" = -999,
      "smbh_Q84" = -999,
      "smbh_lo"  = 0, 
      "smbh_hi"  = 0
    )
  }
  
  results_fit[[kk]] = list("df" = df)
}

SMBH_df = data.frame(do.call(rbind, lapply(results_fit, function(x)x$df)))
names(SMBH_df) = paste0(c(
  "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
  "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
))

format( eddTime / median((1-epsilon)/epsilon)*median(SMBH_df$LEddQ50[SMBH_df$LEddQ50 != -999]), scientific = T)

eddTime / median((1-0.1)/0.1)*median(0.5)

# SMBH_results = do_fit(M0 = 1e2, epsilon = epsilon, t0 = t0)
# SMBH_df_m02 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
# names(SMBH_df_m02) = paste0("m02_",c(
#   "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
#   "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
# ))
# 
# SMBH_results = do_fit(M0 = 1e3, epsilon = epsilon, t0 = t0)
# SMBH_df_m03 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
# names(SMBH_df_m03) = paste0("m03_",c(
#   "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
#   "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
# ))
# 
# SMBH_results = do_fit(M0 = 1e4, epsilon = epsilon, t0 = t0)
# SMBH_df_m04 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
# names(SMBH_df_m04) = paste0("m04_",c(
#   "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
#   "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
# ))
# 
# SMBH_df = dplyr::bind_cols(
#   SMBH_df_m02,
#   SMBH_df_m03,
#   SMBH_df_m04
# )

fwrite(
  SMBH_df,
  "~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_SMBH_results_PDF.csv"
)

# read in data, and also data from BHAR/SFR ratio

SMBH_df = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_SMBH_results_PDF.csv"))
SMBH_df_fixed = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_SMBH_results_fixed.csv"))
BHARSFR_df = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_BHAR_SFR_ratio.csv"))

SMBH_combine = cbind(
  SMBH_df,
  SMBH_df_fixed,
  BHARSFR_df
)
SMBH_combine[is.na(SMBH_combine)] <- -999
SMBH_combine$lowerlim = as.integer(SMBH_combine$lowerlim)

h5file <- path.expand("~/Documents/JWST-MIRI-AGN/data/analysis_data.h5")
# Create the file if it doesn't exist yet; otherwise we append into it
if (!file.exists(h5file)) {
  h5createFile(h5file)
}
if (any(h5ls(h5file)$group == "/SMBH" & h5ls(h5file)$name == "SMBH")) {
  h5delete(h5file, "SMBH/SMBH")
}
h5createGroup(h5file, "SMBH")
h5write(SMBH_combine, file = h5file, name = paste0("SMBH/SMBH"))
