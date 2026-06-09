library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
field_stats = data.frame(fread("~/Documents/Epochs/data/depths/smiles_profound_frame_stats.csv"))

lbt_bins = seq(cosdistTravelTime(z = 2, ref = "737"), cosdistTravelTime(z = 20, ref = "737"), 0.3)
zbins = approx(
  x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
  y = seq(0,40,1),
  xout = lbt_bins
)$y
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2
print(zmids)
length(zmids)

sm_bins = seq(6.5, 12.5, 0.5)
sm_mids = sm_bins[1:(length(sm_bins)-1)] + diff(sm_bins)/2

pdf("~/Desktop/stuff.pdf", width = 10, height = 14)
par(
  mfrow = c(5,2), oma = rep(1.5, 4), mar = rep(2.5, 4)
)
for(i in 1:length(zmids)){
  zmin = zbins[i]
  zmax = zbins[i+1]
  zidx = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50 >= zmin & photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50 < zmax
  
  magplot(
    (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[zidx]),
    (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10Q50[zidx]), 
    # photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10Q50[zidx]/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR100Q50[zidx],
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[zidx]),
    pch = 16,
    xlim = 10^c(4.5, 12.5),
    ylim = 10^c(-5.5, 2.5),
    log = "xy", 
    title = "log10(AGNlum)",
    xlab = "Mstar",
    ylab = "SFR10", 
    scale = c(0.9, 1/20),
    titleshift = 1.5
  )
  magerr(
    (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[zidx]),
    (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10Q50[zidx]), 
    xlo = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMass_lo[zidx],
    xhi = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMass_hi[zidx],
    ylo = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10_lo[zidx],
    yhi = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10_hi[zidx],
    col = adjustcolor("black", 0.3)
  )
  legend(
    x = "topleft",
    legend = paste0(round(zmin,3), "< z <", round(zmax,3))
  )
}

dev.off()


magplot(
  photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$UV1500Q50[zidx],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$SFR10Q50[zidx]),
  xlim = c(-24, -16),
  ylim = c(-5.5, 2.5),
)
legend(
  x = "topleft",
  legend = paste0(round(zmin,3), "< z <", round(zmax,3))
)

jvolume = (cosdistCoVol(z = zmax, ref = "737") - cosdistCoVol(z = zmin, ref = "737")) * 1e9 * (field_stats$detect_effarea_amin2)/(4 * pi * (180/pi)^2 * 60^2)

# hh = maghist(
#   log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[zidx]),
#   breaks = sm_bins, 
#   plot = FALSE
# )
# magplot(
#   sm_mids,
#   hh$counts/volume, 
#   log = 'y'
# )


