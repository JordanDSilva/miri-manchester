library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_astro_cat.csv"))

photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_photom = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_fit_photom.csv"))
photoZ_bc03_snorm_norm_eazy_no_tau_birth_photom = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_fit_photom.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))
common_names = common_names[common_names %notin% c("ID", "RA", "Dec")]
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth # start with Stellar+AGN

DIC_thresh = 2.0
deltaDIC = photoZ_bc03_snorm_norm_eazy_no_tau_birth$DIC - photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$DIC
# deltaDIC[is.na(deltaDIC)] = 999
combined_catalogue[which(deltaDIC <= DIC_thresh), common_names] = photoZ_bc03_snorm_norm_eazy_no_tau_birth[which(deltaDIC <= DIC_thresh), common_names] # now set better fit to Stellar
combined_catalogue[which(deltaDIC <= DIC_thresh), uncommon_names] = NA # set the AGN crap to NA

combined_catalogue$DeltaDIC = deltaDIC
combined_catalogue$AGNFitLowerLimit = deltaDIC <= DIC_thresh # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(deltaDIC <= DIC_thresh, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## Now do photom
combined_photom = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_photom
combined_photom[which(deltaDIC <= DIC_thresh), ] = photoZ_bc03_snorm_norm_eazy_no_tau_birth_photom[which(deltaDIC <= DIC_thresh), ]

## remove the photometric flags
robust_catalogue = combined_catalogue
robust_catalogue[phot_cat$ProFound$flag == 1 , c(common_names, uncommon_names, "DeltaDIC", "AGNFitLowerLimit", "FitType")] = NA
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