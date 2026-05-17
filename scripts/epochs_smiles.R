library(Rfits)
library(magicaxis)
library(ProSpect)
library(celestial)
library(scales)
library(dftools)
library(foreach)
library(hyper.fit)
library(matrixStats)
library(data.table)

## match MIRI smiles catalogue in GOODS

filt_names = read.table(
  file = "~/Documents/Epochs/data/throughputs/wavelengths.txt"
)[,2]
wavelengths = read.table(
  file = "~/Documents/Epochs/data/throughputs/wavelengths.txt"
)[,1]

epochs = data.frame(Rfits_read_table("~/Documents/Epochs/data/EPOCHS_MASTER_public_final.fits"))
smiles = data.frame(
  Rfits_read_table("~/Documents/Epochs/data/hlsp_smiles_jwst_miri_goodss_multi_v1.0_catalog.fits", ext = 4),
  Rfits_read_table("~/Documents/Epochs/data/hlsp_smiles_jwst_miri_goodss_multi_v1.0_catalog.fits", ext = 5)
  # Rfits_read_table("~/Downloads/hlsp_smiles_jwst_miri_goodss_multi_v1.0_catalog.fits", ext = 3)
)
epochs$mag_aper_f444w = -2.5*log10(epochs$flux_aper_f444W)+31.4 ## nJy
smiles$mag_F560W_KRON = -2.5*log10(smiles$F560W_KRON)+31.4 ## nJy

smiles_nirspec = Rfits_read_table(
  "~/Documents/Epochs/data/hlsp_smiles_jwst_nirspec_goodss_multi_v1.0_catalog.fits"
)

match_idx = coordmatch(
  coordref = epochs[, c("RA", "DEC")],
  coordcompare = smiles[, c("RA", "DEC")], 
  rad = 0.3
)

fluxes = cbind( 
  epochs[match_idx$bestmatch$refID, grep("flux_aper", names(epochs), value = TRUE)],
  smiles[match_idx$bestmatch$compareID, grep("^F.+KRON$", names(smiles), value = TRUE)]
)
flux_errs = cbind( 
  epochs[match_idx$bestmatch$refID, grep("fluxerr_aper", names(epochs), value = TRUE)],
  smiles[match_idx$bestmatch$compareID, grep("^F.+KRON_en", names(smiles), value = TRUE)]
)

dim(fluxes)

df_combined = data.frame(
  "RA_epochs" = epochs$RA[match_idx$bestmatch$refID],
  "DEC_epochs" = epochs$DEC[match_idx$bestmatch$refID],
  "zbest_epochs" = epochs$zbest[match_idx$bestmatch$refID],
  "RA_smiles" = smiles$RA[match_idx$bestmatch$compareID],
  "DEC_smiles" = smiles$DEC[match_idx$bestmatch$compareID],
  "ID_smiles" = smiles$ID[match_idx$bestmatch$compareID],
  "epochs_sfr10_zgauss" = epochs$sfr_10Myr_pipes_zgauss[match_idx$bestmatch$refID],
  "epochs_sfr100_zgauss" = epochs$sfr_100Myr_pipes_zgauss[match_idx$bestmatch$refID],
  "epochs_stellar_mass_zgauss" = epochs$stellar_mass_pipes_zgauss[match_idx$bestmatch$refID],
  "epochs_sfr10_zfix" = epochs$sfr_10Myr_pipes_zfix[match_idx$bestmatch$refID],
  "epochs_sfr100_zfix" = epochs$sfr_100Myr_pipes_zfix[match_idx$bestmatch$refID],
  "epochs_stellar_mass_zfix" = epochs$stellar_mass_pipes_zfix[match_idx$bestmatch$refID]
)
# 
magplot(
  (df_combined$RA_epochs - df_combined$RA_smiles)/cos(median(df_combined$DEC_epochs) * pi / 180 * 3600),
  df_combined$DEC_epochs - df_combined$DEC_smiles
)
points(
  (df_combined$RA_epochs[1] - df_combined$RA_smiles[1])/cos(median(df_combined$DEC_epochs) * pi / 180 * 3600),
  df_combined$DEC_epochs[1] - df_combined$DEC_smiles[1], 
  pch =16
)

smiles_nirspec[smiles_nirspec$MIRI_ID %in% df_combined$ID_smiles, c("RA", "DEC", "Z_SPEC")]
df_combined[df_combined$ID_smiles %in% smiles_nirspec$MIRI_ID, c("RA_epochs", "DEC_epochs", "zbest_epochs")]

convert_freq2wave(flux_freq = Jansky2CGS(unlist(fluxes[kk,])*1e-9), wave = wavelengths)/1e-20

kk = 1
convert_freq2wave(flux_freq = Jansky2CGS(unlist(fluxes[kk,])*1e-9), wave = wavelengths)/1e-20

magplot(
  log10(wavelengths), 
  convert_freq2wave(flux_freq = Jansky2CGS(unlist(fluxes[kk,])*1e-9), wave = wavelengths)/1e-20,
  ylim = c(0, 2.0),
  xlim = c(3.8, 5.3)
)
magerr(
  wavelengths/1e4, 
  unlist(fluxes[kk,]),
  ylo = unlist(flux_errs[kk,])
)

# magplot(
#   wavelengths/1e4, 
#   unlist(fluxes[kk,]),
#   ylim = c(0, 2.0)
# )
# magerr(
#   wavelengths/1e4, 
#   unlist(fluxes[kk,]),
#   ylo = unlist(flux_errs[kk,])
# )

names(fluxes) = unlist(filt_names)
names(flux_errs) = paste0(unlist(filt_names), "_ERR")

df_phot = data.frame(
  "ID" = 1:dim(fluxes)[1],
  fluxes, 
  flux_errs
)

fwrite(
  df_phot, 
  file = "~/Documents/Epochs/data/smiles_epochs_match_cat.csv"
)
fwrite(
  df_combined, 
  file = "~/Documents/Epochs/data/smiles_epochs_match_astro_cat.csv"
)


## load in Epochsv2
epochsV2 = Rfits_read(
  filename = "~/Documents/Epochs/EPOCHS-v2_MASTER_Sel-F277W+F356W+F444W_v12_psfmatch_F444W_empirical+v13_psfmatch_F444W_empirical+v14_psfmatch_F444W_empirical.fits"
)

goods_south_idx = epochsV2$OBJECTS$SURVEY %in% c("JADES-DR3-GS-North", "JADES-DR3-GS-South", "JADES-DR3-GS-East", "JADES-DR3-GS-West")

GS_epochsV2 = data.frame(
  epochsV2$OBJECTS[goods_south_idx, ],
  epochsV2$EAZY_SFHZ[goods_south_idx, ],
  epochsV2$PROPERTIES_EAZY_SFHZ[goods_south_idx, ],
  epochsV2$SELECTION[goods_south_idx, ]
)
fwrite(
  GS_epochsV2, 
  file = "~/Documents/Epochs/data/EPOCHS-v2_GS_Sel-F277W+F356W+F444W_v14_psfmatch_F444W_empirical+v14_psfmatch_F444W_empirical.csv"
)
