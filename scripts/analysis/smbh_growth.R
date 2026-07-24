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

Msun = 1.98840987e+30
speedOfLight = 299792458
eddTime = ((1.26e38 * 1e-7) / (Msun * speedOfLight^2))^-1 / (60 * 60 * 24 * 365.2425) # in years

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

#Reines and Volonetri+15
x  <- seq(4, 11.5, 0.01)
dx <- x - 11

line_ll <- (7.45 - 0.08) + (1.05 - 0.11) * dx
line_lh <- (7.45 - 0.08) + (1.05 + 0.11) * dx
line_hl <- (7.45 + 0.08) + (1.05 - 0.11) * dx
line_hh <- (7.45 + 0.08) + (1.05 + 0.11) * dx

upper <- pmax(line_ll, line_lh, line_hl, line_hh)
lower <- pmin(line_ll, line_lh, line_hl, line_hh)

x  <- seq(4, 11.5, 0.01)
dx <- log10(10^x/(3e10))

line_ll <- (7.56 - 0.09) + (1.39 - 0.13) * dx
line_lh <- (7.56 - 0.09) + (1.39 + 0.13) * dx
line_hl <- (7.56 + 0.09) + (1.39 - 0.13) * dx
line_hh <- (7.56 + 0.09) + (1.39 + 0.13) * dx

upper <- pmax(line_ll, line_lh, line_hl, line_hh)
lower <- pmin(line_ll, line_lh, line_hl, line_hh)

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
M0 = 2
epsilon = 0.1

CMANiters = 10000
LDNiters = 50000
do_fit = function(M0, epsilon = 0.1, t0 = 0){
  results_fit = list()
  for(kk in 1:nrow(prospect_cat)){
    if(prospect_cat$FitType[kk] == "StellarAGN"){
      
      # read in posterior
      read_gal = readRDS(
        paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", kk, ".rds")
      )
      
      uni_age = cosdistUniAgeAtz(z = prospect_cat$zQ50[kk], ref = "737")
      # uni_age_lo = uni_age - cosdistUniAgeAtz(z = prospect_cat$zQ84[kk], ref = "737")
      # uni_age_hi = cosdistUniAgeAtz(z = prospect_cat$zQ16[kk], ref = "737") - uni_age
      
      Lagn_level = log10(prospect_cat$AGNlumQ50[kk])
      yerr_hi = log10(prospect_cat$AGNlumQ84[kk]) - Lagn_level
      yerr_lo = Lagn_level- log10(prospect_cat$AGNlumQ16[kk])
      
      Lagn_post = log10(read_gal$post[,c("AGNlum")])
      dens1d <- density(Lagn_post, bw = "SJ")
      eval_dens_1d <- function(x_query) {
        approx(x = dens1d$x, y = dens1d$y, xout = x_query, yleft = 0, yright = 0)$y
      }
      
      Data = list(
        "yy" = Lagn_level,
        "uni_age" = uni_age, 
        "M0" = M0,
        "t0" = t0,
        "epsilon" = epsilon
      )

      like_fun = function(parm, Data){
        
        Lagn_model = log10(
          Lagn_func(
            t = Data$uni_age, 
            M0 = 10^Data$M0, 
            lambda = 10^parm[1], 
            epsilon = Data$epsilon, 
            t0 = Data$t0
          )
        )
        
        SMBH_model = SMBH_func(
          t = Data$uni_age, 
          M0 = 10^Data$M0, 
          lambda = 10^parm[1], 
          epsilon = Data$epsilon, 
          t0 = Data$t0
        )
        
        jd = eval_dens_1d(
          Lagn_model
        )
        
        jd = pmax(jd, 1e-300)
        
        LL = sum(log(jd), na.rm = TRUE)
        
        Lprior = erf(
          SMBH_model/Data$M0
        )*100 - 100 # require growth from init mass
          
        return(
          LL + Lprior
        )
      }
      
      highout = Highlander(
        parm = c(-2),
        Data = Data,
        lower = c(-6),
        upper = c(4),
        liketype = "max",
        likefunc = like_fun,
        Niters = c(CMANiters, LDNiters)
      )

      SMBH_post = SMBH_func(
        M0 = 10^Data$M0,
        t0 = Data$t0,
        epsilon = Data$epsilon,
        t = Data$uni_age,
        lambda = 10^highout$LD_last$Posterior1[,1]
      )
      
      df = c(
        "Q50" = quantile(10^highout$LD_last$Posterior1, probs = 0.5, na.rm = TRUE),
        "Q16" = quantile(10^highout$LD_last$Posterior1, probs = 0.16, na.rm = TRUE),
        "Q84" = quantile(10^highout$LD_last$Posterior1, probs = 0.84, na.rm = TRUE),
        "lo" = quantile(10^highout$LD_last$Posterior1, probs = 0.5, na.rm = TRUE) - 
          quantile(10^highout$LD_last$Posterior1, probs = 0.16, na.rm = TRUE),
        "hi" = quantile(10^highout$LD_last$Posterior1, probs = 0.84, na.rm = TRUE) - 
          quantile(10^highout$LD_last$Posterior1, probs = 0.5, na.rm = TRUE),
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
          "ledd_Q50" = 0,
          "ledd_Q16" = 0,
          "ledd_Q84" = 0,
          "ledd_lo"  = 0, 
          "ledd_hi"  = 0, 
          "smbh_Q50" = 0,
          "smbh_Q16" = 0,
          "smbh_Q84" = 0,
          "smbh_lo"  = 0, 
          "smbh_hi"  = 0
        )
      }
      results_fit[[kk]] = list("df" = df, "highout" = highout)
  }
  results_fit
}

SMBH_results = do_fit(M0 = M0, epsilon = epsilon, t0 = t0)

saveRDS(
  SMBH_results, "~/Desktop/smbh_results.rds"
)

SMBH_df = data.frame(do.call(rbind, lapply(SMBH_results, function(x)x$df)))

zidx = prospect_cat$zQ50 > 1 & prospect_cat$zQ50 < 10
magplot(
  (prospect_cat$StellarMassQ50[zidx]),
  (SMBH_df$smbh_Q50[zidx]),
  z = prospect_cat$zQ50[zidx],
  xlim = 10^c(1,12),
  ylim = 10^c(1,9.5),
  pch = 16,
  log = "xy"
)
magerr(
  (prospect_cat$StellarMassQ50[zidx]),
  (SMBH_df$smbh_Q50[zidx]),
  xlo = prospect_cat$StellarMass_lo[zidx],
  xhi = prospect_cat$StellarMass_hi[zidx],
  ylo = SMBH_df$smbh_lo[zidx],
  hi = SMBH_df$smbh_hi[zidx]
)
abline(0,1)
abline(-3, 1)
polygon(
  10^c(x, rev(x)),
  10^c(lower, rev(upper)),
  col = adjustcolor("grey", 0.3)
)

maghist(
  SMBH_df$smbh_Q50[prospect_cat$FitType == "StellarAGN" & zidx]/prospect_cat$StellarMassQ50[prospect_cat$FitType == "StellarAGN" & zidx],
  log = 'xy'
)
abline(v = 0.002, col = "blue")
