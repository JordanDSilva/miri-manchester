library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_fritz_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_astro_cat.csv"))

photoZ_bc03_snorm_fritz_norm_eazy_photom = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_fit_photom.csv"))
photoZ_bc03_snorm_norm_eazy_no_photom = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_fit_photom.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))
common_names = common_names[common_names %notin% c("ID", "RA", "Dec")]
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth # start with Stellar+AGN

AGN_UV_ratio = 10^(-0.4 * (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$UV1500AGNQ50 - photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$UV1500Q50) )
DIC_thresh = 2.0
deltaDIC = photoZ_bc03_snorm_norm_eazy_no_tau_birth$DIC - photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$DIC

replace_AGN_idx = !(deltaDIC >= DIC_thresh & (AGN_UV_ratio >= 0.1 | photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_photom$AGN_Dale_frac_Q50 >= 0.1))
combined_catalogue[which(replace_AGN_idx), common_names] = photoZ_bc03_snorm_norm_eazy_no_tau_birth[which(replace_AGN_idx), common_names] # now set better fit to Stellar
combined_catalogue[which(replace_AGN_idx), uncommon_names] = NA # set the AGN crap to NA

combined_catalogue$DeltaDIC = deltaDIC
combined_catalogue$AGNFitLowerLimit =  replace_AGN_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(replace_AGN_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## Now do photom
combined_photom = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_photom
combined_photom[which(replace_AGN_idx), ] = photoZ_bc03_snorm_norm_eazy_no_tau_birth_photom[which(replace_AGN_idx), ]

## remove the photometric flags
robust_catalogue = combined_catalogue
robust_catalogue[phot_cat$ProFound$flag == 1 | combined_catalogue$StellarMassQ50 > 10^(12), c(common_names, uncommon_names, "DeltaDIC", "AGNFitLowerLimit", "FitType")] = NA
robust_photom = combined_photom
robust_photom[phot_cat$ProFound$flag == 1 ,] = NA

fwrite(
  robust_catalogue, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)
fwrite(
  robust_photom, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_fit_photom.csv"
)

h5file = "~/Documents/Epochs/data/all_data.h5"