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

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

#Reines and Volonetri+15
x  <- seq(1, 15.5, 0.01)
dx <- x - 11

line_ll <- (7.45 - 0.08) + (1.05 - 0.11) * dx
line_lh <- (7.45 - 0.08) + (1.05 + 0.11) * dx
line_hl <- (7.45 + 0.08) + (1.05 - 0.11) * dx
line_hh <- (7.45 + 0.08) + (1.05 + 0.11) * dx

upper <- pmax(line_ll, line_lh, line_hl, line_hh)
lower <- pmin(line_ll, line_lh, line_hl, line_hh)
fwrite(
  data.frame("mstar" = x, "smbh_lo" = lower, "smbh_hi" = upper),
  file = "~/Documents/Epochs/data/literature/AGN-Catalogues/reinesvolonteri15.csv"
)

x  <- seq(1, 15.5, 0.01)
dx <- log10(10^x/(3e10))

line_ll <- (7.56 - 0.09) + (1.39 - 0.13) * dx
line_lh <- (7.56 - 0.09) + (1.39 + 0.13) * dx
line_hl <- (7.56 + 0.09) + (1.39 - 0.13) * dx
line_hh <- (7.56 + 0.09) + (1.39 + 0.13) * dx

upper <- pmax(line_ll, line_lh, line_hl, line_hh)
lower <- pmin(line_ll, line_lh, line_hl, line_hh)
fwrite(
  data.frame("mstar" = x, "smbh_lo" = lower, "smbh_hi" = upper),
  file = "~/Documents/Epochs/data/literature/AGN-Catalogues/greene20.csv"
)

ding20 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/ding20_smbh.csv"))
maiolino24 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/maiolino24_smbh.csv"))
juodzbalis26 = cbind(
  data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_agn_data_raw.csv")),
  data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_sed_fits_raw.csv"))
)

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

t0 = 0.09850496
M0 = 1e4
epsilon = 0.1

do_fit = function(M0, epsilon = 0.1, t0 = 0){
  results_fit = list()
  for(kk in 1:nrow(prospect_cat)){
    if(prospect_cat$FitType[kk] == "StellarAGN"){
      
      # read in posterior
      read_gal = readRDS(
        paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", kk, ".rds")
      )
      
      uni_age = cosdistUniAgeAtz(z = prospect_cat$zQ50[kk], ref = "737")
      
      Lagn_post = read_gal$post[,c("AGNlum")]
      
      arg_eps = (epsilon / (1e9*uni_age * (1 - epsilon)))
      arg1 = (eddTime * 1.26e38 * M0)^-1 * Lagn_post * arg_eps^-1 
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
      highout = NULL
    }
    
    results_fit[[kk]] = list("df" = df, "highout" = highout)
  }
  results_fit
}

SMBH_results = do_fit(M0 = 1e2, epsilon = epsilon, t0 = t0)
SMBH_df_m02 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
names(SMBH_df_m02) = paste0("m02_",c(
  "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
  "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
))

SMBH_results = do_fit(M0 = 1e3, epsilon = epsilon, t0 = t0)
SMBH_df_m03 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
names(SMBH_df_m03) = paste0("m03_",c(
  "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
  "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
))

SMBH_results = do_fit(M0 = 1e4, epsilon = epsilon, t0 = t0)
SMBH_df_m04 = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))
names(SMBH_df_m04) = paste0("m04_",c(
  "LEddQ50","LEddQ16","LEddQ84","LEdd_lo","LEdd_hi",
  "SMBHQ50","SMBHQ16","SMBHQ84","SMBH_lo","SMBH_hi"
))

SMBH_df = dplyr::bind_cols(
  SMBH_df_m02,
  SMBH_df_m03,
  SMBH_df_m04
)

fwrite(
  SMBH_df,
  "~/Documents/Epochs/data/prospect/prospect_robust_SMBH_results.csv"
)

# read in data, and also data from BHAR/SFR ratio

SMBH_df = data.frame(fread("~/Documents/Epochs/data/prospect/prospect_robust_SMBH_results.csv"))
BHARSFR_df = data.frame(fread("~/Documents/Epochs/data/prospect/prospect_robust_BHAR_SFR_ratio.csv"))

SMBH_combine = cbind(
  SMBH_df,
  BHARSFR_df
)
SMBH_combine[is.na(SMBH_combine)] <- -999
SMBH_combine$lowerlim = as.integer(SMBH_combine$lowerlim)

h5file <- path.expand("~/Documents/Epochs/data/analysis_data.h5")
# Create the file if it doesn't exist yet; otherwise we append into it
if (!file.exists(h5file)) {
  h5createFile(h5file)
}
h5delete(h5file, "SMBH/SMBH")
h5createGroup(h5file, "SMBH")
h5write(SMBH_combine, file = h5file, name = paste0("SMBH/SMBH"))
