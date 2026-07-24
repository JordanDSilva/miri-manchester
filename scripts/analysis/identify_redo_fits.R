library(data.table)
library(Rfits)
library(Rwcs)
library(rhdf5)
library(magicaxis)
library(celestial)

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))
photoZ_bc03_snorm_fritz_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
photoZ_bc03_snorm_skirtor_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_skirtor_norm_eazy_astro_cat.csv"))

filters_info = data.frame(fread(paste0("~/Documents/Epochs/data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0("~/Documents/Epochs/data/throughputs/filt_funcs.rds"))

poor_fits_idx = which(
  photoZ_bc03_snorm_norm_eazy$LPBest <= quantile(photoZ_bc03_snorm_norm_eazy$LPBest, probs = 0.25, na.rm = T) & phot_cat$ProFound$flag == 0L & photoZ_bc03_snorm_norm_eazy$NfinalMCMC != 2000
)

# fwrite(
#   list(which(poor_fits_idx)),
#   file = "~/Documents/Epochs/data/prospect/poor_fits_IDs.txt",
# )

## I had already refitted with more iterations the ones from above, these are gonna be those with Pz integral < 0.6
# foo = fread("~/Documents/Epochs/data/prospect/poor_fits_IDs.txt")
# fwrite(
#   list(poor_fits_idx[(poor_fits_idx %notin% foo$V1)]),
#   file = "~/Documents/Epochs/data/prospect/poor_fits_IDs_V2.txt",
# )
