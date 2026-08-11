library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)
library(bit64)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)
h5file = "~/Documents/Epochs/data/analysis_data.h5"

prospect_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
AGN_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
Stellar_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

# read in literature stuff
stemo20 = data.frame(Rfits_read_table("~/Documents/Epochs/data/literature/AGN-Catalogues/villforth20_ACS_AGN_Catalog.fits"))

barrows21 = data.frame(fread("/Users/jordan/Documents/Epochs/data/literature/AGN-Catalogues/barrows21_wiseAGN.csv"))

juodzbalis26 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_agn_data_raw.csv"))
juodzbalis26_sed = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_sed_fits_raw.csv"))
juodzbalis26_combine = cbind(juodzbalis26, juodzbalis26_sed)


jades_dr4 = Rfits_read_table("~/Documents/Epochs/data/literature/specz-Catalogues/jades_Combined_DR4_external_v1.2.1.fits")
jades_deephst = jades_dr4[jades_dr4$PID == 1210 | jades_dr4$PID == 3215, ]
scholtz25 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/scholtz25_agn.csv"))
head(jades_deephst$NIRSpec_ID)
head(scholtz25$ID)

scholtz25 = scholtz25 |>
  mutate(NIRSpec_ID = as.integer64(sub(".*?(\\d+)$", "\\1", ID)))

scholtz25_merged = jades_deephst |>
  inner_join(scholtz25, by = "NIRSpec_ID")

literature_AGN = data.frame(
  "RA" = c(stemo20$RA, barrows21$RA, juodzbalis26_combine$RA, scholtz25_merged$RA_TARG),
  "DEC" = c(stemo20$DEC, barrows21$DEC, juodzbalis26_combine$Dec, scholtz25_merged$Dec_TARG),
  "z" = c(stemo20$Z, barrows21$zbest, juodzbalis26_combine$z, scholtz25_merged$z),
  
  "Lagn" = c(
    log10(stemo20$L_bol_x_md), 
    log10(barrows21$Lbolagn), 
    juodzbalis26_combine$log_Lbol,
    scholtz25_merged$log10_Lbol
    ),
  "Lagn_lo" = c(
    log10(stemo20$L_bol_x_md)-log10(stemo20$L_bol_x_lo), 
    log10(barrows21$Lbolagn)-log10(barrows21$e_Lbolagn),
    juodzbalis26_combine$log_Lbol_lo,
    rep(0.25, nrow(scholtz25_merged))
    ),
  "Lagn_hi" = c(
    log10(stemo20$L_bol_x_hi)-log10(stemo20$L_bol_x_md), 
    log10(barrows21$E_Lbolagn)-log10(barrows21$Lbolagn),
    juodzbalis26_combine$log_Lbol_hi,
    rep(0.25, nrow(scholtz25_merged))
    ),
  
  "Mstar" = c(
    log10(stemo20$M_star_md), 
    log10(barrows21$Mstar),
    juodzbalis26_combine$log_Mstar_FPHO,
    scholtz25_merged$log10_Mstar
    ),
  "Mstar_lo" = c(
    log10(stemo20$L_bol_x_md)-log10(stemo20$L_bol_x_lo), 
    log10(barrows21$Mstar)-log10(barrows21$e_Mstar),
    juodzbalis26_combine$log_Mstar_FPHO_lo,
    scholtz25_merged$log10_Mstar_lo
    ),
  "Mstar_hi" = c(
    log10(stemo20$L_bol_x_hi)-log10(stemo20$L_bol_x_md), 
    log10(barrows21$E_Mstar)-log10(barrows21$Mstar),
    juodzbalis26_combine$log_Mstar_FPHO_hi,
    scholtz25_merged$log10_Mstar_hi
    ),
  
  "survey" = c(
    rep("Stemo20", nrow(stemo20)), 
    rep("Barrows21", nrow(barrows21)), 
    rep("Juodzbalis26", nrow(juodzbalis26_combine)),
    rep("Scholtz26", nrow(scholtz25_merged))
    )
)

match_cat_idx = coordmatch(
  coordref = literature_AGN[, c("RA", "DEC")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 0.3
)

df = data.frame(
  "ID" = prospect_cat$ID[match_cat_idx$bestmatch$compareID],
  "RA" = prospect_cat$RA[match_cat_idx$bestmatch$compareID],
  "Dec" = prospect_cat$Dec[match_cat_idx$bestmatch$compareID],
  "zQ50" = prospect_cat$zQ50[match_cat_idx$bestmatch$compareID],
  "zsurvey" = literature_AGN$z[match_cat_idx$bestmatch$refID],
  "matchIDs" = match_cat_idx$bestmatch$compareID,
  "survey" = literature_AGN$survey[match_cat_idx$bestmatch$refID],
  
  "FitType" = prospect_cat$FitType[match_cat_idx$bestmatch$compareID],
  
  "AGNlumQ50" = prospect_cat$AGNlumQ50[match_cat_idx$bestmatch$compareID],
  "AGNlumQ16" = prospect_cat$AGNlumQ16[match_cat_idx$bestmatch$compareID],
  "AGNlumQ84" = prospect_cat$AGNlumQ84[match_cat_idx$bestmatch$compareID],
  "AGNlum_lo" = prospect_cat$AGNlum_lo[match_cat_idx$bestmatch$compareID],
  "AGNlum_hi" = prospect_cat$AGNlum_hi[match_cat_idx$bestmatch$compareID],
  
  "StellarMassQ50" = prospect_cat$StellarMassQ50[match_cat_idx$bestmatch$compareID],
  "StellarMassQ16" = prospect_cat$StellarMassQ16[match_cat_idx$bestmatch$compareID],
  "StellarMassQ84" = prospect_cat$StellarMassQ84[match_cat_idx$bestmatch$compareID],
  "StellarMass_lo" = prospect_cat$StellarMass_lo[match_cat_idx$bestmatch$compareID],
  "StellarMass_hi" = prospect_cat$StellarMass_hi[match_cat_idx$bestmatch$compareID],
  
  "StellarAGN_AGNlumQ50" = AGN_cat$AGNlumQ50[match_cat_idx$bestmatch$compareID],
  "StellarAGN_AGNlumQ16" = AGN_cat$AGNlumQ16[match_cat_idx$bestmatch$compareID],
  "StellarAGN_AGNlumQ84" = AGN_cat$AGNlumQ84[match_cat_idx$bestmatch$compareID],
  "StellarAGN_AGNlum_lo" = AGN_cat$AGNlum_lo[match_cat_idx$bestmatch$compareID],
  "StellarAGN_AGNlum_hi" = AGN_cat$AGNlum_hi[match_cat_idx$bestmatch$compareID],
  
  "StellarAGN_StellarMassQ50" = AGN_cat$StellarMassQ50[match_cat_idx$bestmatch$compareID],
  "StellarAGN_StellarMassQ16" = AGN_cat$StellarMassQ16[match_cat_idx$bestmatch$compareID],
  "StellarAGN_StellarMassQ84" = AGN_cat$StellarMassQ84[match_cat_idx$bestmatch$compareID],
  "StellarAGN_StellarMass_lo" = AGN_cat$StellarMass_lo[match_cat_idx$bestmatch$compareID],
  "StellarAGN_StellarMass_hi" = AGN_cat$StellarMass_hi[match_cat_idx$bestmatch$compareID],
  
  "Stellar_StellarMassQ50" = Stellar_cat$StellarMassQ50[match_cat_idx$bestmatch$compareID],
  "Stellar_StellarMassQ16" = Stellar_cat$StellarMassQ16[match_cat_idx$bestmatch$compareID],
  "Stellar_StellarMassQ84" = Stellar_cat$StellarMassQ84[match_cat_idx$bestmatch$compareID],
  "Stellar_StellarMass_lo" = Stellar_cat$StellarMass_lo[match_cat_idx$bestmatch$compareID],
  "Stellar_StellarMass_hi" = Stellar_cat$StellarMass_hi[match_cat_idx$bestmatch$compareID],
  
  "AGNlumQ50_lit" = 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID],
  "AGNlumQ16_lit" = 10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - literature_AGN$Lagn_lo[match_cat_idx$bestmatch$refID]),
  "AGNlumQ84_lit" = 10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] + literature_AGN$Lagn_hi[match_cat_idx$bestmatch$refID]),
  "AGNlum_lo_lit" = 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - (10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - literature_AGN$Lagn_lo[match_cat_idx$bestmatch$refID])),
  "AGNlum_hi_lit" = (10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] + literature_AGN$Lagn_hi[match_cat_idx$bestmatch$refID])) - 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID],
  
  "StellarMassQ50_lit" = 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID],
  "StellarMassQ16_lit" = 10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - literature_AGN$Mstar_lo[match_cat_idx$bestmatch$refID]),
  "StellarMassQ84_lit" = 10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] + literature_AGN$Mstar_hi[match_cat_idx$bestmatch$refID]),
  "StellarMass_lo_lit" = 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - (10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - literature_AGN$Mstar_lo[match_cat_idx$bestmatch$refID])),
  "StellarMass_hi_lit" = (10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] + literature_AGN$Mstar_hi[match_cat_idx$bestmatch$refID])) - 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID]
)

h5file = "~/Documents/Epochs/data/analysis_data.h5"
h5createFile(h5file)
h5delete(
  h5file, 
  "literatureAGN"
)
h5createGroup(
  file = h5file,
  "literatureAGN"
)
h5write(
  obj = df,
  file = h5file,
  name = "literatureAGN/compare"
)

magplot(
  df$zQ50,
  df$zsurvey
)
magplot(
  df$AGNlumQ50_lit, 
  df$StellarAGN_AGNlumQ50,
  log = "xy"
)
abline(0,1)

# match_juodzbalis26 = coordmatch(
#   coordref = bind_rows(
#     juodzbalis26[, c("RA", "Dec")] |> rename(DEC = Dec),
#     stemo20[, c("RA", "DEC")]
#   ),
#   coordcompare = prospect_cat[, c("RA", "Dec")],
#   rad = 0.3
# )
# 
# magplot(
#   c(juodzbalis26$z, stemo20$Z)[match_juodzbalis26$bestmatch$refID],
#   (prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID]), 
#   xlab = "Juodzbalis+26 Ha broadline spectroscopic redshift",
#   ylab = "ProSpect Photoz", 
#   xlim = c(0,8),
#   ylim = c(0,8), 
#   grid = FALSE,
#   pch = 16,
#   cex = 1.5,
#   col = adjustcolor("black", 0.5),
# )
# 
# 
# 
# 
# match_juodzbalis26 = coordmatch(
#   coordref = juodzbalis26[,c("RA", "Dec")],
#   coordcompare = prospect_cat[, c("RA", "Dec")],
#   rad = 0.3
# )
# 
# magplot(
#   juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
#   (prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID]), 
#   xlab = "Juodzbalis+26 Ha broadline spectroscopic redshift",
#   ylab = "ProSpect Photoz", 
#   xlim = c(1,8),
#   ylim = c(1,8), 
#   grid = FALSE,
#   pch = 16,
#   cex = 1.5,
#   col = adjustcolor("black", 0.5),
# )
# magerr(
#   juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
#   unname(prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID]),
#   ylo = unname(prospect_cat$z_lo[match_juodzbalis26$bestmatch$compareID]),
#   yhi = unname(prospect_cat$z_hi[match_juodzbalis26$bestmatch$compareID]),
# )
# abline(0,1, lty = 2)
# 
# redshift_outlier_idx = abs(prospect_cat$zQ50[match_juodzbalis26$bestmatch$compareID] - juodzbalis26$z[match_juodzbalis26$bestmatch$refID]) > 1.0
# stellar_better_idx = prospect_cat$FitType[match_juodzbalis26$bestmatch$compareID] == "Stellar"
# 
# magplot(
#   juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
#   log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
#   xlim = c(43, 46),
#   ylim = c(36, 46),
#   xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
#   ylab = "ProSpect AGN Lbol [erg/s]", 
#   grid = FALSE,
#   col = adjustcolor("black", 0.7),
#   pch = 16, 
#   cex = 1.5
# )
# magerr(
#   juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
#   log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
#   xlo = juodzbalis26$log_Lbol_lo[match_juodzbalis26$bestmatch$refID],
#   xhi = juodzbalis26$log_Lbol_hi[match_juodzbalis26$bestmatch$refID],
#   ylo = log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]) -
#     log10(prospect_cat$AGNlumQ16[match_juodzbalis26$bestmatch$compareID]),
#   yhi = log10(prospect_cat$AGNlumQ84[match_juodzbalis26$bestmatch$compareID]) -
#     log10(prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])
# )
# points(
#   juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID][redshift_outlier_idx],
#   log10(AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])[redshift_outlier_idx],
#   pch = 1, cex = 2.5, col = "red", lwd = 1.5
# )
# points(
#   juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID][stellar_better_idx],
#   log10(AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID])[stellar_better_idx],
#   pch = 8, cex = 2.5, col = "cornflowerblue", lwd = 1.5
# )
# abline(0,1, lty = 2)
# legend(
#   x = "bottomright",
#   pch = c(16,8,1),
#   col = c("black", "cornflowerblue", "red"),
#   legend = c("Stellar+AGN better fit", "Stellar better fit \n (Photometric upper limit)", "Redshift error"),
#   cex = 2.0
# )
# 
# prospect_cat$ID[match_juodzbalis26$bestmatch$compareID]
# prospect_cat$FitType[match_juodzbalis26$bestmatch$compareID]
# prospect_cat$AGN_Dale_frac_Q50[match_juodzbalis26$bestmatch$compareID]
# AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]
# 
# magplot(
#   juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID],
#   log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
#   xlim = c(5, 11),
#   ylim = c(5, 11),
#   xlab = "Juodzbalis+26 CIGALE Mstar from ForcePho [Msun]",
#   ylab = "ProSpect Mstar [Msun]", 
#   grid = FALSE,
#   col = adjustcolor("black", 0.7),
#   pch = 16, 
#   cex = 1.5
# )
# magerr(
#   juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID],
#   log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
#   xlo = juodzbalis26_combine$log_Mstar_CIGALE_lo[match_juodzbalis26$bestmatch$refID],
#   xhi = juodzbalis26_combine$log_Mstar_CIGALE_hi[match_juodzbalis26$bestmatch$refID],
#   ylo = log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]) -
#     log10(prospect_cat$StellarMassQ16[match_juodzbalis26$bestmatch$compareID]),
#   yhi = log10(prospect_cat$StellarMassQ84[match_juodzbalis26$bestmatch$compareID]) -
#     log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])
# )
# points(
#   juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID][redshift_outlier_idx],
#   log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])[redshift_outlier_idx],
#   pch = 1, cex = 2.5, col = "red", lwd = 1.5
# )
# points(
#   juodzbalis26_combine$log_Mstar_CIGALE[match_juodzbalis26$bestmatch$refID][stellar_better_idx],
#   log10(prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID])[stellar_better_idx],
#   pch = 8, cex = 2.5, col = "cornflowerblue", lwd = 1.5
# )
# abline(0,1, lty = 2)
# legend(
#   x = "bottomleft",
#   pch = c(16,8,1),
#   col = c("black", "cornflowerblue", "red"),
#   legend = c("Stellar+AGN better fit", "Stellar better fit", "Redshift error"),
#   cex = 2.0
# )
# 
# df = data.frame(
#   "ID" = prospect_cat$ID[match_juodzbalis26$bestmatch$compareID],
#   "RA" = prospect_cat$RA[match_juodzbalis26$bestmatch$compareID],
#   "Dec" = prospect_cat$Dec[match_juodzbalis26$bestmatch$compareID],
#   "matchIDs" = match_juodzbalis26$bestmatch$compareID,
#   
#   "FitType" = prospect_cat$FitType[match_juodzbalis26$bestmatch$compareID],
#   
#   "AGNlumQ50" = prospect_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID],
#   "AGNlumQ16" = prospect_cat$AGNlumQ16[match_juodzbalis26$bestmatch$compareID],
#   "AGNlumQ84" = prospect_cat$AGNlumQ84[match_juodzbalis26$bestmatch$compareID],
#   "AGNlum_lo" = prospect_cat$AGNlum_lo[match_juodzbalis26$bestmatch$compareID],
#   "AGNlum_hi" = prospect_cat$AGNlum_hi[match_juodzbalis26$bestmatch$compareID],
#   
#   "StellarMassQ50" = prospect_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID],
#   "StellarMassQ16" = prospect_cat$StellarMassQ16[match_juodzbalis26$bestmatch$compareID],
#   "StellarMassQ84" = prospect_cat$StellarMassQ84[match_juodzbalis26$bestmatch$compareID],
#   "StellarMass_lo" = prospect_cat$StellarMass_lo[match_juodzbalis26$bestmatch$compareID],
#   "StellarMass_hi" = prospect_cat$StellarMass_hi[match_juodzbalis26$bestmatch$compareID],
#   
#   "StellarAGN_AGNlumQ50" = AGN_cat$AGNlumQ50[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_AGNlumQ16" = AGN_cat$AGNlumQ16[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_AGNlumQ84" = AGN_cat$AGNlumQ84[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_AGNlum_lo" = AGN_cat$AGNlum_lo[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_AGNlum_hi" = AGN_cat$AGNlum_hi[match_juodzbalis26$bestmatch$compareID],
#   
#   "StellarAGN_StellarMassQ50" = AGN_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_StellarMassQ16" = AGN_cat$StellarMassQ16[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_StellarMassQ84" = AGN_cat$StellarMassQ84[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_StellarMass_lo" = AGN_cat$StellarMass_lo[match_juodzbalis26$bestmatch$compareID],
#   "StellarAGN_StellarMass_hi" = AGN_cat$StellarMass_hi[match_juodzbalis26$bestmatch$compareID],
# 
#   "Stellar_StellarMassQ50" = Stellar_cat$StellarMassQ50[match_juodzbalis26$bestmatch$compareID],
#   "Stellar_StellarMassQ16" = Stellar_cat$StellarMassQ16[match_juodzbalis26$bestmatch$compareID],
#   "Stellar_StellarMassQ84" = Stellar_cat$StellarMassQ84[match_juodzbalis26$bestmatch$compareID],
#   "Stellar_StellarMass_lo" = Stellar_cat$StellarMass_lo[match_juodzbalis26$bestmatch$compareID],
#   "Stellar_StellarMass_hi" = Stellar_cat$StellarMass_hi[match_juodzbalis26$bestmatch$compareID],
# 
#   "AGNlumQ50_j26" = 10^juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID],
#   "AGNlumQ16_j26" = 10^(juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID] - juodzbalis26_combine$log_Lbol_lo[match_juodzbalis26$bestmatch$refID]),
#   "AGNlumQ84_j26" = 10^(juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID] + juodzbalis26_combine$log_Lbol_hi[match_juodzbalis26$bestmatch$refID]),
#   "AGNlum_lo_j26" = 10^juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID] - (10^(juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID] - juodzbalis26_combine$log_Lbol_lo[match_juodzbalis26$bestmatch$refID])),
#   "AGNlum_hi_j26" = (10^(juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID] + juodzbalis26_combine$log_Lbol_hi[match_juodzbalis26$bestmatch$refID])) - 10^juodzbalis26_combine$log_Lbol[match_juodzbalis26$bestmatch$refID],
#   
#   "StellarMassQ50_j26" = 10^juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID],
#   "StellarMassQ16_j26" = 10^(juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID] - juodzbalis26_combine$log_Mstar_FPHO_lo[match_juodzbalis26$bestmatch$refID]),
#   "StellarMassQ84_j26" = 10^(juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID] + juodzbalis26_combine$log_Mstar_FPHO_hi[match_juodzbalis26$bestmatch$refID]),
#   "StellarMass_lo_j26" = 10^juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID] - (10^(juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID] - juodzbalis26_combine$log_Mstar_FPHO_lo[match_juodzbalis26$bestmatch$refID])),
#   "StellarMass_hi_j26" = (10^(juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID] + juodzbalis26_combine$log_Mstar_FPHO_hi[match_juodzbalis26$bestmatch$refID])) - 10^juodzbalis26_combine$log_Mstar_FPHO[match_juodzbalis26$bestmatch$refID]
# )
# 
# h5file = "~/Documents/Epochs/data/analysis_data.h5"
# h5createFile(h5file)
# h5createGroup(
#   file = h5file,
#   "juodzbalis26"
# )
# h5delete(
#   h5file, 
#   "juodzbalis26/compare"
# )
# h5write(
#   obj = df,
#   file = h5file,
#   name = "juodzbalis26/compare"
# )
# 
# # albert_cat = data.frame(fread("~/Desktop/Ha_and_NII_broad_line_dataV2.csv"))
# # match_albert = coordmatch(
# #   coordref = albert_cat[,c("ALPHA_J2000", "DELTA_J2000")],
# #   coordcompare = prospect_cat[, c("RA", "Dec")],
# #   rad = 0.3
# # )
# # magplot(
# #   albert_cat$REDSHIFT[match_albert$bestmatch$refID],
# #   (prospect_cat$zQ50[match_albert$bestmatch$compareID]), 
# #   xlab = "Juodzbalis+26 Ha broadline spectroscopic redshift",
# #   ylab = "ProSpect Photoz", 
# #   xlim = c(0,8),
# #   ylim = c(0,8), 
# #   grid = FALSE,
# #   pch = 16,
# #   cex = 1.5,
# #   col = adjustcolor("black", 0.5),
# # )
# # magerr(
# #   albert_cat$REDSHIFT[match_albert$bestmatch$refID],
# #   unname(prospect_cat$zQ50[match_albert$bestmatch$compareID]),
# #   ylo = unname(prospect_cat$z_lo[match_albert$bestmatch$compareID]),
# #   yhi = unname(prospect_cat$z_hi[match_albert$bestmatch$compareID]),
# # )
# # abline(0,1, lty = 2)
# # 
# # magplot(
# #   (albert_cat$L_BH_BOL[match_albert$bestmatch$refID]),
# #   (AGN_cat$AGNlumQ50[match_albert$bestmatch$compareID]),
# #   z = log10(AGN_cat$StellarMassQ50[match_albert$bestmatch$compareID]),
# #   xlim = 10^c(36, 47),
# #   ylim = 10^c(36, 47),
# #   xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
# #   ylab = "ProSpect AGN Lbol [erg/s]", 
# #   grid = FALSE,
# #   col = adjustcolor("black", 0.7),
# #   pch = 16, 
# #   cex = 1.5,
# #   log = "xy"
# # )
# # magerr(
# #   albert_cat$L_BH_BOL[match_albert$bestmatch$refID],
# #   AGN_cat$AGNlumQ50[match_albert$bestmatch$compareID],
# #   xlo = albert_cat$log10_BH_MASS_LOWER_ERR[match_albert$bestmatch$refID],
# #   xhi = albert_cat$log10_BH_MASS_UPPER_ERR[match_albert$bestmatch$refID],
# #   ylo = AGN_cat$AGNlum_lo[match_albert$bestmatch$compareID],
# #   yhi = AGN_cat$AGNlum_lo[match_albert$bestmatch$compareID]
# # )
# # abline(0, 1, lty = 2)
