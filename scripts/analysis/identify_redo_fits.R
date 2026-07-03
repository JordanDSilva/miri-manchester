library(data.table)
library(Rfits)
library(Rwcs)
library(rhdf5)
library(magicaxis)
library(celestial)

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)
combined_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

phot_cat32 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Phot32"
))

filters_info = data.frame(fread(paste0("~/Documents/Epochs/data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0("~/Documents/Epochs/data/throughputs/filt_funcs.rds"))

df_flux = phot_cat32[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = phot_cat32[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

poor_fits_idx = which(
  combined_cat$LPBest <= quantile(combined_cat$LPBest, probs = 0.25, na.rm = T) & phot_cat$ProFound$flag == 0L 
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
