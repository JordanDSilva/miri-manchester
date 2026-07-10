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

mstar_bins = seq(6.5, 12.5, 0.5)
mstar_mids = mstar_bins[-1] - diff(mstar_bins)*0.5

zbins = c(3.5, 4.5, 5.5, 6.5, 7.5, 8.5)
zmids = zbins[-1] - diff(zbins)*0.5

png("~/Documents/Epochs/plots/SMF_Shuntov25_highz_vmax.png", width = 15, height = 10, units = "in", res = 240)
par(mfrow = c(2,3), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(2.5, 2.5, 1.5, 1.5))
for(kk in 1:(length(zbins)-1)){
  zidx = which( combined_cat$zQ50 >= zbins[kk] & combined_cat$zQ50 < zbins[kk + 1])
  Vmax = (cosdistCoVol(z = pmin(zbins[kk+1], Vmax_data$zmax[zidx]), ref = "737")
          - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = (cosdistCoVol(z = zbins[kk+1], ref = "737")
  #         - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$detect_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  
  mstar_zbin = log10(combined_cat$StellarMassQ50[zidx])
  mstar_err_zbin = 0.5*(log10(combined_cat$StellarMassQ84[zidx]) - log10(combined_cat$StellarMassQ16[zidx]))
  
  mlim = quantile(
    mstar_zbin + 0.4 * (-2.5*log10(phot_cat$Phot32$F444W_scaled_flux_err_app_0.32[zidx]/frame_info$F444W_5sigma_depth_32asec)), 
    probs = 0.6
  )
  
  # hh = maghist(
  #   mstar_zbin,
  #   breaks = mstar_bins,
  #   plot = FALSE
  # )
  
  Veff_err = vector(length = length(mstar_mids))
  Veff = vector(length = length(mstar_mids))
  for(i in 1:length(mstar_mids)){
    midx = mstar_zbin >= mstar_bins[i] & mstar_zbin < mstar_bins[i+1] & mstar_zbin >= mlim
    y = sum(
      1/Vmax[midx], na.rm = TRUE
    )
    y_err = sqrt(sum(
      (1/Vmax[midx])^2, na.rm = TRUE
    ))
    Veff[i] = y
    Veff_err[i] = y_err
  }
  magplot(
    NA,
    log = "y",
    xlim = c(7.0, 12.0),
    ylim = c(1e-6, 1)
  )
  lines(
    shuntov25_smf$M_star[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    col = "forestgreen",
    lty = 2,
    lwd = 3
  )
  # points(
  #   shuntov25_smf_points$V2[shuntov25_smf_points$V1 == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
  #   shuntov25_smf_points$V3[shuntov25_smf_points$V1 == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
  #   pch = 1, 
  #   col = "forestgreen",
  #   cex = 1.5
  # )
  # points(foo$mstar, foo$phi, pch = 8, col = "red", cex=2.0)
  # points(
  #   c(7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75),
  #   10^c(-1.57, -1.97, -2.38, -2.74, -3.17, -3.68, -4.23, -4.78, -5.91), 
  #   pch = 2, col = "blue"
  # )
  # lines(
  #   seq(6.5, 12.5, 0.01),
  #   single_schechter(seq(6.5, 12.5, 0.01), parm = c(11.01, -4.52, -1.79)),
  #   col = "blue"
  # )
  lines(
    seq(6.5, 12.5, 0.01),
    single_schechter(seq(6.5, 12.5, 0.01), parm = c(wang35_smf_parms[kk, "M"], log10(wang35_smf_parms[kk, "phi"]), wang35_smf_parms[kk, "alpha"])),
    col = "red",
    lwd = 3
  )
  lines(
    seq(6.5, 12.5, 0.01),
    single_schechter(seq(6.5, 12.5, 0.01), parm = c(weibel24_smf_parms[kk, "M"], (weibel24_smf_parms[kk, "phi"]), weibel24_smf_parms[kk, "alpha"])),
    col = "cornflowerblue",
    lwd = 3
  )
  # magerr(
  #   hh$mids,
  #   hh$counts/(Vmax*diff(hh$breaks)*0.5),
  #   xlo = 0.25,
  #   ylo = sqrt(hh$counts)/(Vmax*diff(hh$breaks)*0.5)
  # )
  # points(
  #   hh$mids,
  #   hh$counts/(Vmax*diff(hh$breaks)*0.5),
  #   cex = 1.5, 
  #   pch = 16, 
  #   col = adjustcolor("black", 0.8)
  # )
  magerr(
    mstar_mids,
    Veff/(diff(mstar_bins)*0.5),
    xlo = unique(diff(mstar_bins)*0.5),
    ylo = Veff_err/(diff(mstar_bins)*0.5)
  )
  points(
    mstar_mids,
    Veff/(diff(mstar_bins)*0.5),
    cex = 1.5, 
    pch = 16, 
    col = adjustcolor("black", 0.8)
  )
  legend(
    x = "topright", 
    legend = c(
      paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))
    ),
    cex = 1.2
  )
  legend(
    x = "topleft", 
    lty = c(NA, 2, 1, 1),
    pch = c(16, NA, NA, NA),
    col = c("black", "forestgreen", "red", "cornflowerblue"),
    legend = c("ProSpect", "Fit from Shuntov+25 \n COSMOS-WEB", "Fit from Wang+25 \n PRIMER MIRI", "Fit from Weibel+24 \n NIRCam"),
    cex = 1.0
  )
}
mtext(
  "log10(Mstar)", side = 1, outer = TRUE, cex = 2.0
)
mtext(
  "Phi(Mstar)", side = 2, outer = TRUE, cex = 2.0
)
dev.off()
 
png("~/Documents/Epochs/plots/SMF_Shuntov25_all.png", width = 15, height = 10, units = "in", res = 240)
zbins = c(0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10.0, 12.0)
par(mfrow = c(3,5), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(2.5, 2.5, 1.5, 1.5))
for(kk in 1:(length(zbins)-1)){
  zidx = which( combined_cat$zQ50 >= zbins[kk] & combined_cat$zQ50 < zbins[kk + 1])
  Vmax = (cosdistCoVol(z = pmin(zbins[kk+1], Vmax_data$zmax[zidx]), ref = "737")
          - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  # Vmax = (cosdistCoVol(z = zbins[kk+1], ref = "737")
  #         - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)
  
  mstar_zbin = log10(combined_cat$StellarMassQ50[zidx])
  mstar_err_zbin = 0.5*(log10(combined_cat$StellarMassQ84[zidx]) - log10(combined_cat$StellarMassQ16[zidx]))
  
  mlim = quantile(
    mstar_zbin + 0.4 * (-2.5*log10(phot_cat$Phot32$F444W_scaled_flux_err_app_0.32[zidx]/frame_info$F444W_5sigma_depth_32asec)), 
    probs = 0.6
  )
  
  Veff_err = vector(length = length(mstar_mids))
  Veff = vector(length = length(mstar_mids))
  for(i in 1:length(mstar_mids)){
    midx = mstar_zbin >= mstar_bins[i] & mstar_zbin < mstar_bins[i+1] & mstar_zbin >= mlim
    y = sum(
      1/Vmax[midx], na.rm = TRUE
    )
    y_err = sqrt(sum(
      (1/Vmax[midx])^2, na.rm = TRUE
    ))
    Veff[i] = y
    Veff_err[i] = y_err
  }
  
  hh = maghist(
    mstar_zbin,
    breaks = mstar_bins,
    plot = FALSE
  )
  magplot(
    NA,
    log = "y",
    xlim = c(6.5, 12.5),
    ylim = c(1e-8, 1e0)
  )
  lines(
    shuntov25_smf$M_star[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    col = "blue", 
    lwd = 3
  )
  magerr(
    shuntov25_smf$M_star[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    ylo = shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))] - 
      shuntov25_smf$Phi_16[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    yhi = shuntov25_smf$Phi_84[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))] - 
      shuntov25_smf$Phi[shuntov25_smf$Redshift == paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))],
    col = adjustcolor("blue", 0.3),
    poly = TRUE,
    border = NA
  )
  lines(
    seq(6.5, 12.5, 0.01),
    double_schechter(
      x = seq(6.5, 12.5, 0.01),
      p = c(
        10.963 + 0.211 * zbins[kk] - 0.094 * zbins[kk]^2,
        -1.627 + 0.428 * zbins[kk] - 0.154 * zbins[kk]^2,
        -0.539 + 0.053 * zbins[kk] - 0.017 * zbins[kk]^2,
        -3.212 + 0.059 * zbins[kk] - 0.105 * zbins[kk]^2,
        -2.592 - 0.515 * zbins[kk] + 0.045 * zbins[kk]^2
      )
    )
  )
  # magerr(
  #   hh$mids,
  #   hh$counts/(Vmax*diff(hh$breaks)*0.5),
  #   xlo = 0.25,
  #   ylo = sqrt(hh$counts)/(Vmax*diff(hh$breaks)*0.5)
  # )
  # points(
  #   hh$mids,
  #   hh$counts/(Vmax*diff(hh$breaks)*0.5),
  #   cex = 1.5, 
  #   pch = 16, 
  #   col = adjustcolor("black", 0.8)
  # )
  magerr(
    mstar_mids,
    Veff/(diff(mstar_bins)*0.5),
    xlo = unique(diff(mstar_bins)*0.5),
    ylo = Veff_err/(diff(mstar_bins)*0.5)
  )
  points(
    mstar_mids,
    Veff/(diff(mstar_bins)*0.5),
    cex = 1.5, 
    pch = 16, 
    col = adjustcolor("black", 0.8)
  )
  legend(
    x = "topright", 
    legend = c(
      paste0(sprintf("%.1f", zbins[kk]), " < z < ", sprintf("%.1f", zbins[kk+1]))
    ),
    cex = 1.2
  )
  legend(
    x = "bottomright", 
    lty = c(NA, 1),
    pch = c(16, NA),
    col = c("black", "blue"),
    legend = c("ProSpect", "Fit from COSMOS-WEB \n Shuntov+25"),
    cex = 1.2
  )
}
mtext(
  "log10(Mstar)", side = 1, outer = TRUE, cex = 2.0
)
mtext(
  "Phi(Mstar)", side = 2, outer = TRUE, cex = 2.0
)
dev.off()