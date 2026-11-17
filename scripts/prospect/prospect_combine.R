library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)

phot_stats = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", "miri", "_profound_multiband_cat.fits"), 
  ext = "PhotStats"
))

raw_cat = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", "miri", "_profound_multiband_cat.fits"), 
  ext = "RawTot"
))
raw_cat16 = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", "miri", "_profound_multiband_cat.fits"), 
  ext = "Raw16"
))
raw_cat32 = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", "miri", "_profound_multiband_cat.fits"), 
  ext = "Raw32"
))

df_flux = raw_cat[, paste0(filter_info$filters, "_fluxt")]
df_flux_err = raw_cat[, paste0(filter_info$filters, "_scaled_fluxt_err")]

photoZ_bc03_snorm_fritz_norm_eazy = data.frame(read.csv("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy = data.frame(read.csv("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))
common_names = common_names[common_names %notin% c("ID", "CATID", "RA", "Dec")]
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy # start with Stellar+AGN

# AGN_UV_ratio = 10^(-0.4 * (photoZ_bc03_snorm_fritz_norm_eazy$UV1500AGNQ50 - photoZ_bc03_snorm_fritz_norm_eazy$UV1500Q50) )
AGN_UV_ratio = photoZ_bc03_snorm_fritz_norm_eazy$UV1500_AGN_ratioQ50
DIC_thresh = 0.0
deltaDIC = photoZ_bc03_snorm_norm_eazy$DIC - photoZ_bc03_snorm_fritz_norm_eazy$DIC

# sum(
#   (AGN_UV_ratio >= 0.1) | (photoZ_bc03_snorm_fritz_norm_eazy$AGN_20_frac_Q50 >= 0.1), na.rm = TRUE
# )
# sum(
#   (AGN_UV_ratio >= 0.1) | (photoZ_bc03_snorm_fritz_norm_eazy$AGN_1_20_frac_Q50 >= 0.1), na.rm = TRUE
# )
# sum(
#   (AGN_UV_ratio >= 0.1) | (photoZ_bc03_snorm_fritz_norm_eazy$AGN_3_20_frac_Q50 >= 0.1), na.rm = TRUE
# )
# sum(
#   (AGN_UV_ratio >= 0.1) | (photoZ_bc03_snorm_fritz_norm_eazy$AGN_Dale_frac_Q50 >= 0.1), na.rm = TRUE
# )

## Isolate a photometrically selected sample
replace_AGN_idx = !((AGN_UV_ratio >= 0.1) | (photoZ_bc03_snorm_fritz_norm_eazy$AGN_3_20_fracQ50 >= 0.1))
# replace_AGN_idx = !(deltaDIC > DIC_thresh)

combined_catalogue[which(replace_AGN_idx), common_names] = photoZ_bc03_snorm_norm_eazy[which(replace_AGN_idx), common_names] # now set better fit to Stellar
combined_catalogue[which(replace_AGN_idx), uncommon_names] = NA # set the AGN crap to NA

combined_catalogue$DeltaDIC = deltaDIC
# combined_catalogue$AGNFitLowerLimit =  replace_AGN_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(replace_AGN_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## remove the photometric flags
robust_catalogue = combined_catalogue

# hainline+26 brown dwarfs
hainline26 <- data.frame(
  object_name = c(
    "JADES-GS-BD-1", "JADES-GS-BD-2", "JADES-GS-BD-3",
    "JADES-GS-BD-4", "JADES-GS-BD-5", "JADES-GS-BD-6",
    "JADES-GS-BD-7", "JADES-GS-BD-8", "JADES-GS-BD-9",
    "JADES-GS-BD-10", "JADES-GS-BD-11", "JADES-GS-BD-12",
    "JADES-GS-BD-13", "JADES-GS-BD-14", "JADES-GS-BD-15",
    "JADES-GS-BD-16", "JADES-GS-BD-17", "PANO-GS-BD-18",
    "PANO-GS-BD-19", "JADES-GS-BD-20", "JADES-GS-BD-21",
    "JADES-GS-BD-22", "JADES-GS-BD-23", "JADES-GS-BD-24",
    "PANO-GS-BD-25"
  ),
  ID = c(
    44120, 20541, 42611, 13985, 190413, 2950, 101709, 92415,
    106126, 192282, 140057, 486694, 469315, 433690, 502828,
    495589, 2301, 441205, 501052, 318023, 7646, 128632,
    153326, 233407, 484735
  ),
  RA = c(
    53.026024, 53.036290, 53.056305, 53.076733, 53.084038,
    53.103905, 53.105540, 53.133760, 53.161140, 53.161800,
    53.200314, 52.938843, 52.985350, 52.989790, 53.005913,
    53.030766, 53.069588, 53.078167, 53.079990, 53.095350,
    53.096580, 53.104652, 53.112896, 53.128980, 53.162940
  ),
  DEC = c(
    -27.867167, -27.883734, -27.868404, -27.889942, -27.839346,
    -27.908504, -27.815155, -27.825521, -27.809166, -27.831615,
    -27.764128, -27.787360, -27.755032, -27.864613, -27.821442,
    -27.841581, -27.910807, -27.927135, -27.942880, -27.921080,
    -27.898750, -27.776674, -27.742144, -27.753494, -27.937887
  )
)

brown_dwarf_match_idx = coordmatch(
  coordref = hainline26[, c("RA", "DEC")],
  coordcompare = combined_catalogue[, c("RA", "Dec")],
  rad = 1.0
)
brown_dwarf_idx = (1:dim(combined_catalogue)[1]) %in% brown_dwarf_match_idx$bestmatch$compareID

message(
  "We have matched: ", sum(brown_dwarf_idx), " BD candidates with Hainline+26"
)


not_robust_idx =  which(phot_cat$ProFound$flag == 1 | combined_catalogue$StellarMassQ50 > 10^(12) | combined_catalogue$UV1500Q50 < -22 | brown_dwarf_idx)

robust_catalogue[not_robust_idx, c(common_names, uncommon_names, "DeltaDIC", "FitType")] = NA

message(
  "We have removed: ", length(which(phot_cat$ProFound$flag == 0 & (combined_catalogue$StellarMassQ50 > 10^(12) | combined_catalogue$UV1500Q50 < -22) & (combined_catalogue$zQ50 <= 10 & combined_catalogue$zQ50 >= 1))), " sources"
)

# compare how many we set as NAs from large Mstar/MUV
combined_catalogue_z = combined_catalogue[which(combined_catalogue$zQ50 >= 1.07467287 & combined_catalogue$zQ50 < 9.98276084), ]
robust_catalogue_z = robust_catalogue[which(robust_catalogue$zQ50 >= 1.07467287 & robust_catalogue$zQ50 < 9.98276084), ]

## Isolate a DIC selected sample
combined_DIC_catalogue = photoZ_bc03_snorm_fritz_norm_eazy # start with Stellar+AGN
DIC_idx = deltaDIC >= DIC_thresh
replace_DIC_idx = !(!replace_AGN_idx & DIC_idx)
combined_DIC_catalogue[which(replace_DIC_idx), common_names] = photoZ_bc03_snorm_norm_eazy[which(replace_DIC_idx), common_names] # now set better fit to Stellar
combined_DIC_catalogue[which(replace_DIC_idx), uncommon_names] = NA # set the AGN crap to NA

combined_DIC_catalogue$DeltaDIC = deltaDIC
# combined_DIC_catalogue$AGNFitLowerLimit =  replace_DIC_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_DIC_catalogue$FitType = ifelse(replace_DIC_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## remove the photometric flags
robust_DIC_catalogue = combined_DIC_catalogue

not_robust_DIC_idx =  which(phot_cat$ProFound$flag == 1 | combined_DIC_catalogue$StellarMassQ50 > 10^(12) | combined_DIC_catalogue$UV1500Q50 < -22 | brown_dwarf_idx)

robust_DIC_catalogue[not_robust_DIC_idx, c(common_names, uncommon_names, "DeltaDIC", "FitType")] = NA

fwrite(
  robust_catalogue, 
  file = "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)

fwrite(
  robust_DIC_catalogue, 
  file = "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robustDIC_catalogue.csv"
)

py_code <- '
import pandas as pd

prospect_with_agn = pd.read_csv(
    "/Users/jordan/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"
)
prospect_without_agn = pd.read_csv(
    "/Users/jordan/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"
)
robust_catalogue = pd.read_csv(
    "/Users/jordan/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)
robust_DIC_catalogue = pd.read_csv(
    "/Users/jordan/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robustDIC_catalogue.csv"
)

h5_path = "/Users/jordan/Documents/JWST-MIRI-AGN/data/analysis_data.h5"

with pd.HDFStore(h5_path, mode="a") as store:
    store.put("ProSpect/StellarAGN", prospect_with_agn, format="table")
    store.put("ProSpect/Stellar", prospect_without_agn, format="table")
    store.put("ProSpect/Robust", robust_catalogue, format="table")
    store.put("ProSpect/RobustDIC", robust_DIC_catalogue, format="table")
'

writeLines(py_code, "/tmp/write_h5.py")
system('source ~/miniconda3/etc/profile.d/conda.sh && conda activate space && python /tmp/write_h5.py')
unlink('/tmp/write_h5.py')


# test quantities from the fits compared with and without AGN component
pdf("~/Documents/JWST-MIRI-AGN/rough-plots/test_stellar_agn_stuff.pdf")
# magplot(
#   photoZ_bc03_snorm_norm_eazy$CATID[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
#   photoZ_bc03_snorm_fritz_norm_eazy$CATID[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]]
# )
magplot(
  photoZ_bc03_snorm_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
  photoZ_bc03_snorm_fritz_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
  # z = photoZ_bc03_snorm_fritz_norm_eazy$AGN_3_20_fracQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
  log = "xy",
  xlim = c(1e7, 1e12),
  ylim = c(1e7, 1e12),
  pch = 16,
  col = adjustcolor("black", 0.3),
  xlab = "Stellarmass from Stellar only fit but not significant AGN",
  ylab = "Stellarmass from Stellar+AGN only fit but not significant AGN"
)
abline(0,1, col = "green")

magplot(
  photoZ_bc03_snorm_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
  photoZ_bc03_snorm_fritz_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]],
  log = "xy",
  xlim = c(1e-5, 1e3),
  ylim = c(1e-5, 1e3),
  pch = 16,
  col = adjustcolor("black", 0.3),
  xlab = "SFR from Stellar only fit but not significant AGN",
  ylab = "SFR from Stellar+AGN only fit but not significant AGN"
)
abline(0,1, col = "green")

maghist(
  log10(
    photoZ_bc03_snorm_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]]/
      photoZ_bc03_snorm_fritz_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]]
  ),
  xlab = "Delta Mstar",
  log = "y"
)
maghist(
  log10(
    photoZ_bc03_snorm_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]]/
      photoZ_bc03_snorm_fritz_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "Stellar")]]
  ),
  xlab = "Delta SFR",
  log = "y"
)

# magplot(
#   photoZ_bc03_snorm_norm_eazy$CATID[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
#   photoZ_bc03_snorm_fritz_norm_eazy$CATID[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]]
# )
magplot(
  photoZ_bc03_snorm_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  photoZ_bc03_snorm_fritz_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  z = photoZ_bc03_snorm_fritz_norm_eazy$AGN_3_20_fracQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  log = "xy",
  xlim = c(1e7, 1e12),
  ylim = c(1e7, 1e12),
  xlab = "Stellarmass from Stellar only fit but significant AGN",
  ylab = "Stellarmass from Stellar+AGN only fit but significant AGN"
)
magplot(
  photoZ_bc03_snorm_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  photoZ_bc03_snorm_fritz_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  z = photoZ_bc03_snorm_fritz_norm_eazy$AGN_3_20_fracQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]],
  log = "xy",
  xlim = c(1e-5, 1e3),
  ylim = c(1e-5, 1e3),
  xlab = "SFR from Stellar only fit but significant AGN",
  ylab = "SFR from Stellar+AGN only fit but significant AGN"
)
maghist(
  log10(
    photoZ_bc03_snorm_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]]/
      photoZ_bc03_snorm_fritz_norm_eazy$StellarMassQ50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]]
  ),
  xlab = "Delta Mstar",
  log = "y"
)
maghist(
  log10(
    photoZ_bc03_snorm_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]]/
      photoZ_bc03_snorm_fritz_norm_eazy$SFR10Q50[robust_catalogue_z$CATID[which(robust_catalogue_z$FitType == "StellarAGN")]]
  ),
  xlab = "Delta SFR",
  log = "y"
)
dev.off()