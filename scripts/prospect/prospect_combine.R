library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_fritz_norm_eazy = data.frame(read.csv("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy = data.frame(read.csv("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))
common_names = common_names[common_names %notin% c("ID", "CATID", "RA", "Dec")]
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy # start with Stellar+AGN

AGN_UV_ratio = 10^(-0.4 * (photoZ_bc03_snorm_fritz_norm_eazy$UV1500AGNQ50 - photoZ_bc03_snorm_fritz_norm_eazy$UV1500Q50) )
DIC_thresh = 0.0
deltaDIC = photoZ_bc03_snorm_norm_eazy$DIC - photoZ_bc03_snorm_fritz_norm_eazy$DIC

## Isolate a photometrically selected sample
replace_AGN_idx = !((AGN_UV_ratio >= 0.1 | photoZ_bc03_snorm_fritz_norm_eazy$AGN_Dale_frac_Q50 >= 0.1))

combined_catalogue[which(replace_AGN_idx), common_names] = photoZ_bc03_snorm_norm_eazy[which(replace_AGN_idx), common_names] # now set better fit to Stellar
combined_catalogue[which(replace_AGN_idx), uncommon_names] = NA # set the AGN crap to NA

combined_catalogue$DeltaDIC = deltaDIC
# combined_catalogue$AGNFitLowerLimit =  replace_AGN_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(replace_AGN_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## remove the photometric flags
robust_catalogue = combined_catalogue
robust_catalogue[which(phot_cat$ProFound$flag == 1 | combined_catalogue$StellarMassQ50 > 10^(12) | combined_catalogue$UV1500Q50 < -22), c(common_names, uncommon_names, "DeltaDIC", "FitType")] = NA


## Isolate a DIC selected sample
combined_DIC_catalogue = photoZ_bc03_snorm_fritz_norm_eazy # start with Stellar+AGN
replace_DIC_idx = !(deltaDIC >= DIC_thresh & (AGN_UV_ratio >= 0.1 | photoZ_bc03_snorm_fritz_norm_eazy$AGN_Dale_frac_Q50 >= 0.1))
combined_DIC_catalogue[which(replace_DIC_idx), common_names] = photoZ_bc03_snorm_norm_eazy[which(replace_DIC_idx), common_names] # now set better fit to Stellar
combined_DIC_catalogue[which(replace_DIC_idx), uncommon_names] = NA # set the AGN crap to NA

combined_DIC_catalogue$DeltaDIC = deltaDIC
# combined_DIC_catalogue$AGNFitLowerLimit =  replace_DIC_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_DIC_catalogue$FitType = ifelse(replace_DIC_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## remove the photometric flags
robust_DIC_catalogue = combined_DIC_catalogue
robust_DIC_catalogue[which(phot_cat$ProFound$flag == 1 | combined_DIC_catalogue$StellarMassQ50 > 10^(12) | combined_DIC_catalogue$UV1500Q50 < -22), c(common_names, uncommon_names, "DeltaDIC", "FitType")] = NA

fwrite(
  robust_catalogue, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)

fwrite(
  robust_DIC_catalogue, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robustDIC_catalogue.csv"
)

py_code <- '
import pandas as pd

prospect_with_agn = pd.read_csv(
    "/Users/jordan/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"
)
prospect_without_agn = pd.read_csv(
    "/Users/jordan/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"
)
robust_catalogue = pd.read_csv(
    "/Users/jordan/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)
robust_DIC_catalogue = pd.read_csv(
    "/Users/jordan/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robustDIC_catalogue.csv"
)

h5_path = "/Users/jordan/Documents/Epochs/data/analysis_data.h5"

with pd.HDFStore(h5_path, mode="w") as store:
    store.put("ProSpect/StellarAGN", prospect_with_agn, format="table")
    store.put("ProSpect/Stellar", prospect_without_agn, format="table")
    store.put("ProSpect/Robust", robust_catalogue, format="table")
    store.put("ProSpect/RobustDIC", robust_DIC_catalogue, format="table")
'

writeLines(py_code, "/tmp/write_h5.py")
system('source ~/miniconda3/etc/profile.d/conda.sh && conda activate space && python /tmp/write_h5.py')
unlink('/tmp/write_h5.py')
