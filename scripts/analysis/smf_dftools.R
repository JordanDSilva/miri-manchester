library(celestial)
library(ProSpect)
library(Rfits)
library(data.table)
library(magicaxis)
library(foreach)
library(Highlander)
library(matrixStats)

main_dir = "~/Documents/Epochs/"
fname = "smiles"

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

combined_cat = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
))
eazy_stats = data.frame(fread(
  "~/Documents/Epochs/data/eazy/smiles_agn_blue_sfhz/smiles_eazy_stats_agn_blue_sfhz.csv"
))

Vmax_data = data.frame(fread(paste0(main_dir, "/data/depths/", "zmax_5sigma_data.csv")))

shuntov25_smf = data.frame(fread(
  "~/Downloads/SMF_fits.ecsv", skip = 12
))
shuntov25_smf_points = data.frame(fread(
  "~/Downloads/SMF_data_points.ecsv", skip = 12
))
wang35_smf_parms = data.frame(
  "z" = c(3,4,5,6,7,8),
  "M" = c(11, 10.94, 10.81, 10.54, 10.74, 10.12),
  "phi" = c(13.6, 4.49, 1.44, 1.06, 0.41, 1.51)*1e-5,
  "alpha" = c(-1.57, -1.77, -2.0, -2.1, -2.13, -2.0)
)
weibel24_smf_parms = data.frame(
  "z" = c(4,5,6,7,8),
  "M" = c(11.01, 10.26, 10.01, 10.0, 10.0),
  "phi" = c(-4.52, -4.07, -4.26, -4.36, -4.86),
  "alpha" = c(-1.79, -1.86, -1.95, -1.93, -2.16)
)

single_schechter <- function(x, parm) {
  M0    <- parm[1]
  phi1  <- 10^parm[2]
  alpha <- parm[3]
  mu = 10^(x - M0)
  phi = log(10)*phi1*mu^(alpha+1)*exp(-mu)
  
  return(phi)
}
double_schechter = function(x, p){
  
  ## best starting  parm = c(11, -1, -1, -2, -3), 
  
  Mstar = p[1]
  alpha1 = p[2]
  alpha2 = p[3]
  phi1 = 10^p[4]
  phi2 = 10^p[5]
  
  mu = 10^x / 10^Mstar
  
  ret = log(10) * exp(-1*mu) * (phi1 * mu^(alpha1 + 1) + phi2 * mu^(alpha2 + 1))
  return( ret )
}

lbt_bins = c(seq(8.0, 13.0, 0.5))
zbins = approx(
  x = cosdistTravelTime(z = seq(0,40,0.001), ref = "737"),
  y = seq(0,40,0.001),
  xout = lbt_bins
)$y
# zbins = seq(0.0, 12.0, 0.5)
zbins = c(0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10.0, 12.0)
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

mstar_bins = seq(6.5, 12.5, 0.5)
mstar_mids = mstar_bins[-1] - diff(mstar_bins)*0.5

fit_parm = list()
for(kk in 1:(length(zbins)-1)){
  zidx = which( combined_cat$zQ50 >= zbins[kk] & combined_cat$zQ50 < zbins[kk + 1])
  Vmax = (cosdistCoVol(z = pmin(zbins[kk+1], Vmax_data$zmax[zidx]), ref = "737")
          - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = (cosdistCoVol(z = zbins[kk+1], ref = "737")
  #         - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = rep(Vmax, length(zidx))
  
  mstar_zbin = log10(combined_cat$StellarMassQ50[zidx])
  mstar_err_zbin = 0.5*(log10(combined_cat$StellarMassQ84[zidx]) - log10(combined_cat$StellarMassQ16[zidx]))
  
  mlim = quantile(
    mstar_zbin + 0.4 * (-2.5*log10(phot_cat$Phot32$F444W_scaled_flux_err_app_0.32[zidx]/frame_info$F444W_5sigma_depth_32asec)), 
    probs = 0.6
  )
  mlim_up = 100.0
  obs.selection = function(x){
    ifelse(
      x <= mlim | x >= mlim_up,
      0,
      1
    )
  }
  
  idx = mstar_zbin >= mlim & mstar_zbin < mlim_up
  
  if(kk >= 8){
    pfunc = function(parm){
      logprior = sum(
        c(
          # ifelse(parm[1] < 9.0, -9999999, 0.0),
          # ifelse(parm[1] > 12.0, -9999999, 0.0),
          dnorm(parm[1], mean = 11.0, sd = 1.0, log = TRUE),
          ifelse(parm[2] < -8.0, -9999999, 0.0),
          ifelse(parm[2] > -1.0, -9999999, 0.0),
          dnorm(parm[3], mean = -2.0, sd = 0.3, log = TRUE),
          0
        ),
        na.rm = TRUE
      )
      return(logprior)
    }
  }else{
    pfunc = function(parm){
      logprior = sum(
        c(
          ifelse(parm[1] < 9.0, -9999999, 0.0),
          ifelse(parm[1] > 12.0, -9999999, 0.0),
          ifelse(parm[2] < -8.0, -9999999, 0.0),
          ifelse(parm[2] > -1.0, -9999999, 0.0),
          ifelse(parm[3] < -2.5, -9999999, 0.0),
          ifelse(parm[3] > -1.0, -9999999, 0.0),
          0
        ),
        na.rm = TRUE
      )
      return(logprior)
    }
  }
  
  df = dftools::dffit(
    x = mstar_zbin[idx], 
    x.err = mstar_err_zbin[idx],
    selection = Vmax[idx], 
    obs.selection = obs.selection,
    gdf = single_schechter,
    dx = 0.01, 
    # xmin = 6,
    # xmax = 12,
    n.bootstrap = NULL,
    n.iterations = 1000,
    p.initial = c(10, -4, -0.5), 
    method = "Nelder-Mead", 
    prior = pfunc
  )
  dftools::dfplot(df, xlim = c(5, 13), ylim = c(1e-7, 1e-1))
  df$fit$p.best
  fit_parm[[kk]] = df$fit$p.best
  
  lines(
    shuntov25_smf$M_star[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    col = "forestgreen",
    lty = 2,
    lwd = 3
  )
  lines(
    seq(6.5, 12.5, 0.01),
    single_schechter(seq(6.5, 12.5, 0.01), parm = c(wang35_smf_parms[1, "M"], log10(wang35_smf_parms[1, "phi"]), wang35_smf_parms[1, "alpha"])),
    col = "red",
    lwd = 3
  )
  
}

parms = do.call(rbind, fit_parm)





fit_parm = list()
for(kk in 1:(length(zbins)-1)){
  zidx = which( combined_cat$zQ50 >= zbins[kk] & combined_cat$zQ50 < zbins[kk + 1])
  Vmax = (cosdistCoVol(z = pmin(zbins[kk+1], Vmax_data$zmax[zidx]), ref = "737")
          - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = (cosdistCoVol(z = zbins[kk+1], ref = "737")
  #         - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = rep(Vmax, length(zidx))
  
  mstar_zbin = log10(combined_cat$AGNlumQ50[zidx])
  mstar_err_zbin = 0.5*(log10(combined_cat$AGNlumQ84[zidx]) - log10(combined_cat$AGNlumQ16[zidx]))
  
  mlim = quantile(
    mstar_zbin + 0.4 * (-2.5*log10(phot_cat$Phot32$F444W_scaled_flux_err_app_0.32[zidx]/frame_info$F444W_5sigma_depth_32asec)), 
    probs = 0.5, 
    na.rm = TRUE
  )
  mlim_up = 100.0
  obs.selection = function(x){
    ifelse(
      x <= mlim | x >= mlim_up,
      0,
      1
    )
  }
  
  idx = mstar_zbin >= mlim & mstar_zbin < mlim_up & !is.na(mstar_zbin)
  
  if(kk >= 8){
    pfunc = function(parm){
      logprior = sum(
        c(
          # ifelse(parm[1] < 9.0, -9999999, 0.0),
          # ifelse(parm[1] > 12.0, -9999999, 0.0),
          dnorm(parm[1], mean = 43, sd = 1.0, log = TRUE),
          ifelse(parm[2] < -8.0, -9999999, 0.0),
          ifelse(parm[2] > -1.0, -9999999, 0.0),
          dnorm(parm[3], mean = -2.0, sd = 0.3, log = TRUE),
          0
        ),
        na.rm = TRUE
      )
      return(logprior)
    }
  }else{
    pfunc = function(parm){
      logprior = sum(
        c(
          ifelse(parm[1] < 35, -9999999, 0.0),
          ifelse(parm[1] > 48, -9999999, 0.0),
          ifelse(parm[2] < -8.0, -9999999, 0.0),
          ifelse(parm[2] > -1.0, -9999999, 0.0),
          ifelse(parm[3] < -2.5, -9999999, 0.0),
          ifelse(parm[3] > -1.0, -9999999, 0.0),
          0
        ),
        na.rm = TRUE
      )
      return(logprior)
    }
  }
  
  if(sum(idx) >= 2){
    df = dftools::dffit(
      x = mstar_zbin[idx], 
      x.err = mstar_err_zbin[idx],
      selection = Vmax[idx], 
      obs.selection = obs.selection,
      gdf = single_schechter,
      dx = 0.01, 
      # xmin = 6,
      # xmax = 12,
      n.bootstrap = NULL,
      n.iterations = 1000,
      p.initial = c(44, -7, -0.5), 
      method = "Nelder-Mead", 
      prior = pfunc
    )
    dftools::dfplot(df, xlim = c(41, 49), ylim = c(1e-7, 1e-1))
    fit_parm[[kk]] = df$fit$p.best
  }else{
    fit_parm[[kk]] = rep(0,3)
  }
}
