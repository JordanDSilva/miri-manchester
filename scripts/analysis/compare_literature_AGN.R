library(data.table)
library(Rfits)
library(Rwcs)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)
library(bit64)
library(RANN)

filter_info = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)
h5file = "~/Documents/JWST-MIRI-AGN/data/analysis_data.h5"

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
prospect_cat_DIC = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robustDIC_catalogue.csv"))
AGN_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
Stellar_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

SMBH_df = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/smbh/prospect_robust_SMBH_results_PDF.csv"))

sample_idx = prospect_cat$zQ50 >= 1.07467287 & prospect_cat$zQ50 < 9.98276084 & !is.na(prospect_cat$zQ50)

prospect_cat = prospect_cat[sample_idx, ]
prospect_cat_DIC = prospect_cat_DIC[sample_idx, ]

AGN_cat = AGN_cat[sample_idx, ]
Stellar_cat = Stellar_cat[sample_idx, ]
SMBH_df = SMBH_df[sample_idx, ]

smiles_footprint = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")
midis_footprint = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/midis/midis_midis_midis_profound_stack.fits")
nircam_footprint = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/smiles/warped/JADES-DR3-GS_F277W_smiles_warp.fits")


# read in literature stuff
stemo20 = data.frame(Rfits_read_table("~/Documents/JWST-MIRI-AGN/data/literature/AGN-Catalogues/villforth20_ACS_AGN_Catalog.fits"))
stemo20 = stemo20[stemo20$L_bol_x_md > 0, ]

barrows21 = data.frame(fread("/Users/jordan/Documents/JWST-MIRI-AGN/data/literature/AGN-Catalogues/barrows21_wiseAGN.csv"))

juodzbalis26 = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/literature/AGN-Catalogues/juodzbalis26_agn_data_raw.csv"))
juodzbalis26_sed = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/literature/AGN-Catalogues/juodzbalis26_sed_fits_raw.csv"))
juodzbalis26_combine = cbind(juodzbalis26, juodzbalis26_sed)

jades_dr4 = Rfits_read_table("~/Documents/JWST-MIRI-AGN/data/literature/specz-Catalogues/jades_Combined_DR4_external_v1.2.1.fits")
jades_deephst = jades_dr4[jades_dr4$PID == 1210 | jades_dr4$PID == 3215, ]
scholtz25 = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/literature/AGN-Catalogues/scholtz25_agn.csv"))
head(jades_deephst$NIRSpec_ID)
head(scholtz25$ID)

scholtz25 = scholtz25 |>
  mutate(NIRSpec_ID = as.integer64(sub(".*?(\\d+)$", "\\1", ID)))

scholtz25_merged = jades_deephst |>
  inner_join(scholtz25, by = "NIRSpec_ID")


# vadim_catalogue = data.frame(fread("~/Downloads/sample_final.csv"))
vadim_catalogue = data.frame(fread("~/Downloads/catalog_agn_classification.csv"))
vadim_own = vadim_catalogue
# vadim_own = vadim_catalogue[vadim_catalogue$flag == "OWN", ]

match_cat_vadim_idx = coordmatch(
  coordref = vadim_own[, c("ra", "dec")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 1.0
)

vadim_match = vadim_own[match_cat_vadim_idx$bestmatch$refID, ]
prospect_vadim_match = prospect_cat[match_cat_vadim_idx$bestmatch$compareID, ]
vadim_match_combine = cbind(vadim_match[, c("ra", "dec", "class_own", "class", "flag", "redshift")], prospect_vadim_match[, c("ID", "CATID", "FitType", "SFR10Q50", "StellarMassQ50", "AGNlumQ50", "zQ50", "DIC")], "DIC_Stellar" = Stellar_cat$DIC[match_cat_vadim_idx$bestmatch$compareID], "DIC_StellarAGN" = AGN_cat$DIC[match_cat_vadim_idx$bestmatch$compareID], "AGNlumQ50_AGN" = AGN_cat$AGNlumQ50[match_cat_vadim_idx$bestmatch$compareID])

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
    rep("Scholtz25", nrow(scholtz25_merged))
    ),
  
  "ledd" = c(
    rep(NA, nrow(stemo20)), 
    rep(NA, nrow(barrows21)), 
    juodzbalis26_combine$lambda_Edd,
    rep(NA, nrow(scholtz25_merged))
  ),
  
  "SMBHmass" = c(
    rep(NA, nrow(stemo20)), 
    rep(NA, nrow(barrows21)), 
    juodzbalis26_combine$log_MBH,
    rep(NA, nrow(scholtz25_merged))
  )
)

# ## match with Above AGN table
# vadim_match_literature_idx = coordmatch(
#   coordref = vadim_own[, c("ra", "dec")],
#   coordcompare = literature_AGN[, c("RA", "DEC")],
#   rad = 1.0
# )

# fwrite(
#   vadim_own[vadim_match_literature_idx$bestmatch$refID, ], 
#   "~/Desktop/final_sample_matched_miri.csv"
# )

literature_AGN$vadim_class = rep(NA, dim(literature_AGN)[1])
literature_AGN$vadim_class[vadim_match_literature_idx$bestmatch$compareID] = vadim_own$class[vadim_match_literature_idx$bestmatch$refID]

match_cat_idx = coordmatch(
  coordref = literature_AGN[, c("RA", "DEC")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 1.0
)

df = data.frame(
  "ID" = prospect_cat$ID[match_cat_idx$bestmatch$compareID],
  "CATID" = prospect_cat$CATID[match_cat_idx$bestmatch$compareID],
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
  
  "Stellar_DIC" = Stellar_cat$DIC[match_cat_idx$bestmatch$compareID],
  "StellarAGN_DIC" = AGN_cat$DIC[match_cat_idx$bestmatch$compareID],
  
  "AGNlumQ50_lit" = 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID],
  "AGNlumQ16_lit" = 10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - literature_AGN$Lagn_lo[match_cat_idx$bestmatch$refID]),
  "AGNlumQ84_lit" = 10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] + literature_AGN$Lagn_hi[match_cat_idx$bestmatch$refID]),
  "AGNlum_lo_lit" = 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - (10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] - literature_AGN$Lagn_lo[match_cat_idx$bestmatch$refID])),
  "AGNlum_hi_lit" = (10^(literature_AGN$Lagn[match_cat_idx$bestmatch$refID] + literature_AGN$Lagn_hi[match_cat_idx$bestmatch$refID])) - 10^literature_AGN$Lagn[match_cat_idx$bestmatch$refID],
  
  "Ledd_lit" = literature_AGN$ledd[match_cat_idx$bestmatch$refID],
  "SMBHmass_lit" = literature_AGN$SMBHmass[match_cat_idx$bestmatch$refID],
  "LeddQ50" = SMBH_df$LEddQ50[match_cat_idx$bestmatch$compareID],
  "SMBHQ50" = SMBH_df$SMBHQ50[match_cat_idx$bestmatch$compareID],
  
  "StellarMassQ50_lit" = 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID],
  "StellarMassQ16_lit" = 10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - literature_AGN$Mstar_lo[match_cat_idx$bestmatch$refID]),
  "StellarMassQ84_lit" = 10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] + literature_AGN$Mstar_hi[match_cat_idx$bestmatch$refID]),
  "StellarMass_lo_lit" = 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - (10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] - literature_AGN$Mstar_lo[match_cat_idx$bestmatch$refID])),
  "StellarMass_hi_lit" = (10^(literature_AGN$Mstar[match_cat_idx$bestmatch$refID] + literature_AGN$Mstar_hi[match_cat_idx$bestmatch$refID])) - 10^literature_AGN$Mstar[match_cat_idx$bestmatch$refID],
  
  "class_from_vadim" = literature_AGN$vadim_class[match_cat_idx$bestmatch$refID]
)

sum(df$FitType == "StellarAGN" & df$class_from_vadim == "LRD", na.rm = TRUE)/sum(df$class_from_vadim == "LRD", na.rm = TRUE)
sum(df$FitType == "StellarAGN" & df$class_from_vadim == "LBD", na.rm = TRUE)/sum(df$class_from_vadim == "LBD", na.rm = TRUE)


magplot(
  df$Ledd_lit,
  df$LeddQ50,
  xlim = c(0, 2.5),
  ylim = c(0, 2.5)
)
abline(0,1)


magplot(
  df$LeddQ50,
  df$AGNlumQ50_lit  / (10^df$SMBHmass_lit*1.26e38),
  xlim = c(0, 2.5),
  ylim = c(0, 2.5)
)
abline(0,1)

magplot(
  df$SMBHmass_lit,
  log10(df$SMBHQ50),
  xlim = c(5, 8),
  ylim = c(5, 8)
)
abline(0,1)

df$StellarAGN_AGNlumQ50[df$FitType == "StellarAGN" & df$class_from_vadim == "LBD"]/df$StellarAGN_StellarMassQ50[df$FitType == "StellarAGN" & df$class_from_vadim == "LBD"]

df$Stellar_DIC[df$FitType == "StellarAGN" & df$class_from_vadim == "LRD"]
df$StellarAGN_DIC[df$FitType == "StellarAGN" & df$class_from_vadim == "LRD"]

h5file = "~/Documents/JWST-MIRI-AGN/data/analysis_data.h5"
# Only create the file if it doesn't already exist
if (!file.exists(h5file)) {
  h5createFile(h5file)
}

# Only delete the group if it already exists inside the file
if ("literatureAGN" %in% h5ls(h5file)$name) {
  h5delete(h5file, "literatureAGN")
}

h5createGroup(
  file = h5file,
  "literatureAGN"
)

h5write(
  obj = df,
  file = h5file,
  name = "literatureAGN/compare"
)

write.table(
  df$CATID[df$FitType == "StellarAGN"],
  "~/Documents/JWST-MIRI-AGN/data/prospect/match_better_fit_AGN.txt"
)


colors = rep(adjustcolor("yellow", 0.7), 19755)

# Have to see how many of the actual candidates from all of these catalogues are actually within the smiles footprint, that we could detect
# png("~/Desktop/smiles_footprint.png", width = nrow(smiles_footprint$image[,]), height= ncol(smiles_footprint$image[,]), res = 240)
# plot(smiles_footprint$image[,], coord.type = "deg", sparse = 1, locut = 0.5, hicut = 0.999)
# ProFound::profoundSegimPlot(
#   image = matrix(0, nrow(smiles_footprint$image[,]), ncol(smiles_footprint$image[,])),
#   segim = smiles_footprint$segim[,]$imDat,
#   sparse = 1, 
#   col = rep(adjustcolor("yellow", 0.7), 19755),
#   add = TRUE
# )
# points(
#   Rwcs_s2p(RA = juodzbalis26_combine$RA, Dec = juodzbalis26_combine$Dec, keyvalues = smiles_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("red", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Juodzbalis26"], Dec = df$Dec[df$survey == "Juodzbalis26"], keyvalues = smiles_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "red", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = stemo20$RA, Dec = stemo20$DEC, keyvalues = smiles_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("green", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Stemo20"], Dec = df$Dec[df$survey == "Stemo20"], keyvalues = smiles_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "green", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = scholtz25_merged$RA_TARG, Dec = scholtz25_merged$Dec_TARG, keyvalues = smiles_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("cyan", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Scholtz26"], Dec = df$Dec[df$survey == "Scholtz26"], keyvalues = smiles_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "cyan", 
#   cex = 2.0
# )
# dev.off()
# 
# png("~/Desktop/midis_footprint.png", width = nrow(midis_footprint$image[,]), height= ncol(midis_footprint$image[,]), res = 240)
# plot(midis_footprint$image[,], coord.type = "deg", sparse = 1, locut = 0.5, hicut = 0.999)
# ProFound::profoundSegimPlot(
#   image = matrix(0, nrow(midis_footprint$image[,]), ncol(midis_footprint$image[,])),
#   segim = midis_footprint$segim[,]$imDat,
#   sparse = 1, 
#   col = rep(adjustcolor("yellow", 0.7), 19755),
#   add = TRUE
# )
# points(
#   Rwcs_s2p(RA = juodzbalis26_combine$RA, Dec = juodzbalis26_combine$Dec, keyvalues = midis_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("red", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Juodzbalis26"], Dec = df$Dec[df$survey == "Juodzbalis26"], keyvalues = midis_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "red", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = stemo20$RA, Dec = stemo20$DEC, keyvalues = midis_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("green", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Stemo20"], Dec = df$Dec[df$survey == "Stemo20"], keyvalues = midis_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "green", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = scholtz25_merged$RA_TARG, Dec = scholtz25_merged$Dec_TARG, keyvalues = midis_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("cyan", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Scholtz26"], Dec = df$Dec[df$survey == "Scholtz26"], keyvalues = midis_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "cyan", 
#   cex = 2.0
# )
# dev.off()
# 
# png("~/Desktop/nircam_footprint.png", width = nrow(nircam_footprint$image[,]), height= ncol(nircam_footprint$image[,]), res = 240)
# plot(nircam_footprint$image[,], coord.type = "deg", sparse = 1, locut = 0.5, hicut = 0.999)
# ProFound::profoundSegimPlot(
#   image = matrix(0, nrow(smiles_footprint$image[,]), ncol(smiles_footprint$image[,])),
#   segim = smiles_footprint$segim[,]$imDat,
#   sparse = 1, 
#   col = rep(adjustcolor("yellow", 0.7), 19755),
#   add = TRUE
# )
# points(
#   Rwcs_s2p(RA = juodzbalis26_combine$RA, Dec = juodzbalis26_combine$Dec, keyvalues = nircam_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("red", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Juodzbalis26"], Dec = df$Dec[df$survey == "Juodzbalis26"], keyvalues = nircam_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "red", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = stemo20$RA, Dec = stemo20$DEC, keyvalues = nircam_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("green", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Stemo20"], Dec = df$Dec[df$survey == "Stemo20"], keyvalues = nircam_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "green", 
#   cex = 2.0
# )
# points(
#   Rwcs_s2p(RA = scholtz25_merged$RA_TARG, Dec = scholtz25_merged$Dec_TARG, keyvalues = nircam_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("cyan", 0.4)
# )
# points(
#   Rwcs_s2p(RA = df$RA[df$survey == "Scholtz26"], Dec = df$Dec[df$survey == "Scholtz26"], keyvalues = nircam_footprint$image$keyvalues), 
#   pch = "*", 
#   col = "cyan", 
#   cex = 2.0
# )
# dev.off()


# magplot(
#   df$zQ50,
#   df$zsurvey
# )
# magplot(
#   df$AGNlumQ50_lit/10, 
#   df$StellarAGN_AGNlumQ50,
#   log = "xy"
# )
# abline(0,1)
# 
# median(
#   log10(df$AGNlumQ50_lit/10/df$StellarAGN_AGNlumQ50)
# )




check_if_in_data = function(data_xy, test_xy, rad = 1.0){
  
  test_xy[is.na(test_xy)] = -9999999
  nn = nn2(data  = data_xy,
           query = cbind(test_xy[,1], test_xy[,2]),         # no rounding needed
           k = 1, searchtype = "radius", radius = rad)
  
  foo        = nn$nn.idx[,1]                      # row of `pix`; 0 if none within k
  foo[foo == 0] = NA
  near_valid = !is.na(foo)
  
  # round to the nearest pixel (pixel centres are at integers in FITS convention)
  key_obj  = paste(round(test_xy[,1]), round(test_xy[,2]))
  key_grid = paste(data_xy[,1], data_xy[,2])
  
  foo = match(key_obj, key_grid)   # row of smiles_data_idx, NA if not in the valid footprint
  in_footprint = !is.na(foo)
  
  return(
    list(
      "rann" = near_valid,
      "exact" = in_footprint
    )
  )
}

smiles_data_idx = which(!is.na(smiles_footprint$image[,]$imDat), arr.ind = TRUE)
foobar = check_if_in_data(
  data_xy = smiles_data_idx,
  test_xy = Rwcs_s2p(RA = vadim_own$ra, Dec = vadim_own$dec, keyvalues = smiles_footprint$image$keyvalues)
)
sum(foobar$rann)
which(foobar$rann)
dim(vadim_own[foobar$rann, c("redshift", "ra", "dec", "flag", "class_own", "class")])
dim(vadim_match[, c("ra", "dec")])

in_footprint_only_coords <- anti_join(vadim_own[foobar$rann, c("ra", "dec")], vadim_match[, c("ra", "dec")], by = c("ra", "dec"))

plot_cutout = function(frame, RA, Dec, box = 200){
  plot(frame$image[RA, Dec, box = box, type = "coord"], coord.type = "deg", sparse = 1, locut = 0.5, hicut = 0.999)
  
  pix_coords = Rwcs_s2p(
    RA = RA, 
    Dec = Dec,
    keyvalues = frame$image$keyvalues
  )
  
  if(!is.null(frame$segim)){
    ProFound::profoundSegimPlot(
      image = matrix(0, nrow(frame$image[RA, Dec, box = box, type = "coord"]), ncol(frame$image[RA, Dec, box = box, type = "coord"])),
      segim = frame$segim[pix_coords[1], pix_coords[2], box = box]$imDat,
      sparse = 1, 
      col = rep(adjustcolor("yellow", 0.7), 19755),
      add = TRUE
    )
  }
  
  points(box/2.0, box/2.0, pch = 1, col = "red", cex = 1.5)
}
par(mfrow = c(2,2), mar = rep(0,4), oma = rep(0,4))
plot_cutout(
  frame = smiles_footprint,
  RA = in_footprint_only_coords$ra[1],
  Dec = in_footprint_only_coords$dec[1],
  box = 200
)
plot_cutout(
  frame = nircam_footprint,
  RA = in_footprint_only_coords$ra[1],
  Dec = in_footprint_only_coords$dec[1],
  box = 200
)
plot_cutout(
  frame = smiles_footprint,
  RA = in_footprint_only_coords$ra[2],
  Dec = in_footprint_only_coords$dec[2],
  box = 200
)
plot_cutout(
  frame = nircam_footprint,
  RA = in_footprint_only_coords$ra[2],
  Dec = in_footprint_only_coords$dec[2],
  box = 200
)

for(id in prospect_vadim_match$CATID){
  file.copy(
    paste0("~/Documents/JWST-MIRI-AGN/plots/fitted_candidates/gal_", id, ".png"),
    paste0("~/Documents/JWST-MIRI-AGN/plots/vadim_match/cutouts/gal_", id, ".png")
  )
}

# smiles_footprint = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")
# png("~/Desktop/smiles_footprint_vadim.png", width = nrow(smiles_footprint$image[,]), height= ncol(smiles_footprint$image[,]), res = 240)
# plot(smiles_footprint$image[,], coord.type = "deg", sparse = 1, locut = 0.5, hicut = 0.999)
# ProFound::profoundSegimPlot(
#   image = matrix(0, nrow(smiles_footprint$image[,]), ncol(smiles_footprint$image[,])),
#   segim = smiles_footprint$segim[,]$imDat,
#   sparse = 1, 
#   col = rep(adjustcolor("yellow", 0.7), 19755),
#   add = TRUE
# )
# points(
#   Rwcs_s2p(RA = vadim_own$ra, Dec = vadim_own$dec, keyvalues = smiles_footprint$image$keyvalues), 
#   pch = 16, 
#   col = adjustcolor("red", 0.8)
# )
# dev.off()

sum(vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "LBD")/sum(vadim_match_combine$class == "LBD")
sum(vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "LRD")/sum(vadim_match_combine$class == "LRD")
sum(vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "AGN")/sum(vadim_match_combine$class == "AGN")
sum(vadim_match_combine$FitType == "StellarAGN")/dim(vadim_match_combine)[1]

vadim_match_combine$CATID[vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class_own == "LBD"]
vadim_match_combine$CATID[vadim_match_combine$FitType == "Stellar" & vadim_match_combine$class_own == "LBD"]

table(vadim_match_combine$class_own)
table(vadim_match_combine$FitType)


log10(vadim_match_combine$AGNlumQ50[vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "LBD"])
log10(vadim_match_combine$AGNlumQ50_AGN[vadim_match_combine$FitType == "Stellar" & vadim_match_combine$class == "LBD"])

log10(vadim_match_combine$AGNlumQ50[vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "LRD"])
log10(vadim_match_combine$AGNlumQ50_AGN[vadim_match_combine$FitType == "Stellar" & vadim_match_combine$class == "LRD"])

log10(vadim_match_combine$AGNlumQ50[vadim_match_combine$FitType == "StellarAGN" & vadim_match_combine$class == "AGN"])
log10(vadim_match_combine$AGNlumQ50_AGN[vadim_match_combine$FitType == "Stellar" & vadim_match_combine$class == "AGN"])

magplot(
  prospect_cat$zQ50[prospect_cat$FitType == "StellarAGN"],
  log10(prospect_cat$AGNlumQ50[prospect_cat$FitType == "StellarAGN"]/prospect_cat$StellarMassQ50[prospect_cat$FitType == "StellarAGN"]),
  col = adjustcolor("red", 0.1),
  pch = 16,
  ylim = c(25, 48)
)
points(
  prospect_cat$zQ50[prospect_cat$FitType == "Stellar"],
  log10(AGN_cat$AGNlumQ50[prospect_cat$FitType == "Stellar"]/AGN_cat$StellarMassQ50[prospect_cat$FitType == "Stellar"]),
  col = adjustcolor("blue", 0.1),
  pch = 16,
)
