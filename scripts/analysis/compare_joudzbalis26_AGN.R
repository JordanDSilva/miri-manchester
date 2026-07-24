library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)
h5file = "~/Documents/Epochs/data/analysis_data.h5"

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
AGN_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))

juodzbalis26 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_agn_data_raw.csv"))
juodzbalis26_sed = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_sed_fits_raw.csv"))
juodzbalis26_combine = cbind(juodzbalis26, juodzbalis26_sed)

match_juodzbalis26 = coordmatch(
  coordref = juodzbalis26[,c("RA", "Dec")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 0.3
)

magplot(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  (prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID]), 
  xlab = "Juodzbalis+26 Ha broadline spectroscopic redshift",
  ylab = "ProSpect Photoz", 
  xlim = c(1,8),
  ylim = c(1,8), 
  grid = FALSE,
  pch = 16,
  cex = 1.5,
  col = adjustcolor("black", 0.5),
)
magerr(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  unname(prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID]),
  ylo = unname(prospect_cat$z_lo[match_juodzbalis26$bestmatch$compareID]),
  yhi = unname(prospect_cat$z_hi[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)

redshift_outlier_idx = abs(prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID] - juodzbalis26$z[match_juodzbalis26$bestmatch$refID]) > 1.0
stellar_better_idx = prospect_cat$FitType[match_juodzbalis26$bestmatch$compareID] == "Stellar"

magplot(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlim = c(43, 46),
  ylim = c(36, 46),
  xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
  ylab = "ProSpect AGN Lbol [erg/s]", 
  grid = FALSE,
  col = adjustcolor("black", 0.7),
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26$log_Lbol_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26$log_Lbol_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(prospect_cat$AGNlumQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(prospect_cat$AGNlumQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])
)
points(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID][redshift_outlier_idx],
  log10(AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])[redshift_outlier_idx],
  pch = 1, cex = 2.5, col = "red", lwd = 1.5
)
points(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID][stellar_better_idx],
  log10(AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])[stellar_better_idx],
  pch = 8, cex = 2.5, col = "cornflowerblue", lwd = 1.5
)
abline(0,1, lty = 2)
legend(
  x = "bottomright",
  pch = c(16,8,1),
  col = c("black", "cornflowerblue", "red"),
  legend = c("Stellar+AGN better fit", "Stellar better fit \n (Photometric upper limit)", "Redshift error"),
  cex = 2.0
)

magplot(
  juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID],
  log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  xlim = c(5, 11),
  ylim = c(5, 11),
  xlab = "Juodzbalis+26 CIGALE Mstar from ForcePho [Msun]",
  ylab = "ProSpect Mstar [Msun]", 
  grid = FALSE,
  col = adjustcolor("black", 0.7),
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID],
  log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26_combine$log_Mstar_CIGALE_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26_combine$log_Mstar_CIGALE_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(prospect_cat$StellarMassQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(prospect_cat$StellarMassQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])
)
points(
  juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID][redshift_outlier_idx],
  log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])[redshift_outlier_idx],
  pch = 1, cex = 2.5, col = "red", lwd = 1.5
)
points(
  juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID][stellar_better_idx],
  log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])[stellar_better_idx],
  pch = 8, cex = 2.5, col = "cornflowerblue", lwd = 1.5
)
abline(0,1, lty = 2)
legend(
  x = "bottomleft",
  pch = c(16,8,1),
  col = c("black", "cornflowerblue", "red"),
  legend = c("Stellar+AGN better fit", "Stellar better fit", "Redshift error"),
  cex = 2.0
)

# h5createFile(h5file)


albert_cat = data.frame(fread("~/Desktop/Ha_and_NII_broad_line_dataV2.csv"))
match_albert = coordmatch(
  coordref = albert_cat[,c("ALPHA_J2000", "DELTA_J2000")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 0.3
)
magplot(
  albert_cat$REDSHIFT[match_albert$bestmatch$refID],
  (prospect_cat$zQ50[match_albert$bestmatch$compareID]), 
  xlab = "Juodzbalis+26 Ha broadline spectroscopic redshift",
  ylab = "ProSpect Photoz", 
  xlim = c(0,8),
  ylim = c(0,8), 
  grid = FALSE,
  pch = 16,
  cex = 1.5,
  col = adjustcolor("black", 0.5),
)
magerr(
  albert_cat$REDSHIFT[match_albert$bestmatch$refID],
  unname(prospect_cat$zQ50[match_albert$bestmatch$compareID]),
  ylo = unname(prospect_cat$z_lo[match_albert$bestmatch$compareID]),
  yhi = unname(prospect_cat$z_hi[match_albert$bestmatch$compareID]),
)
abline(0,1, lty = 2)

magplot(
  (albert_cat$L_BH_BOL[match_albert$bestmatch$refID]),
  (AGN_cat$AGNlumQ50[match_albert$bestmatch$compareID]),
  z = log10(AGN_cat$StellarMassQ50[match_albert$bestmatch$compareID]),
  xlim = 10^c(36, 47),
  ylim = 10^c(36, 47),
  xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
  ylab = "ProSpect AGN Lbol [erg/s]", 
  grid = FALSE,
  col = adjustcolor("black", 0.7),
  pch = 16, 
  cex = 1.5,
  log = "xy"
)
magerr(
  albert_cat$L_BH_BOL[match_albert$bestmatch$refID],
  AGN_cat$AGNlumQ50[match_albert$bestmatch$compareID],
  xlo = albert_cat$log10_BH_MASS_LOWER_ERR[match_albert$bestmatch$refID],
  xhi = albert_cat$log10_BH_MASS_UPPER_ERR[match_albert$bestmatch$refID],
  ylo = AGN_cat$AGNlum_lo[match_albert$bestmatch$compareID],
  yhi = AGN_cat$AGNlum_lo[match_albert$bestmatch$compareID]
)
abline(0, 1, lty = 2)
