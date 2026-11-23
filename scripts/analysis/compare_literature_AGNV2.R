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

harikane23_muv = c(-19.50, -18.94, -19.20, -17.82, -21.66, -16.71, -19.50, -21.23, -21.49, -19.61, -19.72, -22.4)

harikane23_muv = c(-19.50, -18.94, -19.20, -17.82, -21.66, -16.71, -19.50, -21.23, -21.49, -19.61, -19.72, -22.4)
library(tibble)

harikane23 <- tribble(
  ~name,          ~L_Ha_broad, ~L_Ha_broad_lo, ~L_Ha_broad_hi, ~L_ratio, ~L_ratio_lo, ~L_ratio_hi, ~FWHM, ~FWHM_lo, ~FWHM_hi, ~M_BH,  ~M_BH_lo, ~M_BH_hi, ~L_bol, ~L_bol_lo, ~L_bol_hi, ~lambda_Edd, ~lambda_Edd_lo, ~lambda_Edd_hi,
  "CEERS_01244",  7.7e42, 0.2e42, 0.2e42,  2.51, 0.22, 0.17,  2228,  52, 75,    3.2e7, 0.2e7, 0.3e7,  1.9e45, 0.6e45, 3.8e45,  0.46, 0.16, 1.04,
  "GLASS_160133", 1.2e42, 0.1e42, 0.1e42,  0.26, 0.01, 0.02,  1028,  13, 19,    2.3e6, 0.1e6, 0.1e6,  1.1e45, 0.9e45, 6.1e45,  3.68, 2.84, 21.21,
  "GLASS_150029", 8.2e41, 0.5e41, 0.7e41,  0.22, 0.01, 0.02,  1429,  67, 10,    3.7e6, 0.3e6, 0.2e6,  9.1e44, 7.1e44, 13.8e44, 1.89, 1.49, 3.31,
  "CEERS_00746",  6.8e43, 0.4e43, 0.4e43,  5.06, 1.26, 1.10,  1660, 162, 157,   5.8e7, 1.3e7, 1.5e7,  1.1e46, 0.2e46, 0.1e46,  1.47, 0.52, 0.60,
  "CEERS_01665",  6.8e42, 0.6e42, 0.9e42,  0.28, 0.05, 0.05,  1794, 171, 282,   1.9e7, 0.5e7, 0.8e7,  4.8e45, 3.6e45, 7.2e45,  1.93, 1.59, 4.61,
  "CEERS_00672",  1.8e43, 0.2e43, 0.2e43,  1.95, 0.40, 0.45,  2208, 241, 277,   5.0e7, 1.3e7, 1.7e7,  4.3e45, 1.4e45, 0.4e45,  0.65, 0.33, 0.31,
  "CEERS_02782",  7.6e42, 0.7e42, 0.8e42,  0.85, 0.10, 0.14,  2534, 266, 260,   4.2e7, 1.0e7, 1.2e7,  2.8e45, 1.4e45, 2.3e45,  0.51, 0.32, 0.71,
  "CEERS_00397",  2.6e42, 0.4e42, 0.8e42,  0.19, 0.05, 0.07,  1739, 317, 359,   1.0e7, 0.5e7, 0.8e7,  2.8e45, 2.2e45, 14.2e45, 2.02, 1.80, 20.68,
  "CEERS_00717",  1.2e42, 0.2e42, 0.3e42,  0.57, 0.11, 0.17,  6279, 881, 805,   9.8e7, 3.2e7, 4.4e7,  7.1e44, 4.3e44, 51.8e44, 0.06, 0.04, 0.63,
  "CEERS_01236",  4.8e41, 0.8e41, 1.1e41,  2.39, 0.59, 1.12,  3521, 485, 649,   1.8e7, 0.6e7, 1.0e7,  1.8e44, 0.5e44, 2.6e44,  0.08, 0.04, 0.20
)

maghist(
  juodzbalis26$log_Lbol, breaks = seq(41, 46, 0.5), freq = TRUE
)
maghist(
  log10(harikane23$L_bol), breaks = seq(41, 46, 0.5), add = TRUE, col = adjustcolor("blue", 0.3), freq = TRUE
)
maghist(
  log10(prospect_cat$AGNlumQ50[prospect_cat$zQ50 > 3]), breaks = seq(41, 46, 0.5), add = TRUE, col = adjustcolor("red", 0.3), freq = TRUE
)

maghist(
  harikane23_muv, add = TRUE, col = adjustcolor("blue", 0.3)
)
maghist(
  juodzbalis26_combine$M_UV
)

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
vadim_catalogue = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/literature/vadim26_prismAGN/catalog_agn_classification.csv"))
vadim_own = vadim_catalogue
# vadim_own = vadim_catalogue[vadim_catalogue$flag == "OWN", ]

match_cat_vadim_idx = coordmatch(
  coordref = vadim_own[, c("ra", "dec")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 1.0
)

vadim_match = vadim_own[match_cat_vadim_idx$bestmatch$refID, ] 
prospect_vadim_match = prospect_cat[match_cat_vadim_idx$bestmatch$compareID, ]
vadim_match_combine = cbind(vadim_match[, c("ra", "dec", "class", "redshift")], prospect_vadim_match[, c("ID", "CATID", "FitType", "SFR10Q50", "StellarMassQ50", "AGNlumQ50", "zQ50", "DIC")], "DIC_Stellar" = Stellar_cat$DIC[match_cat_vadim_idx$bestmatch$compareID], "DIC_StellarAGN" = AGN_cat$DIC[match_cat_vadim_idx$bestmatch$compareID], "AGNlumQ50_AGN" = AGN_cat$AGNlumQ50[match_cat_vadim_idx$bestmatch$compareID])

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
vadim_match_literature_idx = coordmatch(
  coordref = vadim_own[, c("ra", "dec")],
  coordcompare = literature_AGN[, c("RA", "DEC")],
  rad = 1.0
)

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

dim(df)
sum(df$FitType == "StellarAGN")

sum(df$FitType == "StellarAGN" & df$class_from_vadim == "LRD", na.rm = TRUE)/sum(df$class_from_vadim == "LRD", na.rm = TRUE)
sum(df$FitType == "StellarAGN" & df$class_from_vadim == "LBD", na.rm = TRUE)/sum(df$class_from_vadim == "LBD", na.rm = TRUE)
sum(df$FitType == "StellarAGN" & df$class_from_vadim == "AGN", na.rm = TRUE)/sum(df$class_from_vadim == "AGN", na.rm = TRUE)

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

# write.table(
#   df$CATID[df$FitType == "StellarAGN"],
#   "~/Documents/JWST-MIRI-AGN/data/prospect/match_better_fit_AGN.txt"
# )