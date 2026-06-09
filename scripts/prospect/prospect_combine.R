library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_astro_cat.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth), names(photoZ_bc03_snorm_norm_eazy_no_tau_birth))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth # start with Stellar+AGN

deltaDIC = photoZ_bc03_snorm_norm_eazy_no_tau_birth$DIC - photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$DIC
# deltaDIC[is.na(deltaDIC)] = 999
combined_catalogue[which(deltaDIC <= 2), common_names] = photoZ_bc03_snorm_norm_eazy_no_tau_birth[which(deltaDIC <= 2), common_names] # now set better fit to Stellar
combined_catalogue$DeltaDIC = deltaDIC
combined_catalogue$AGNFitLowerLimit = deltaDIC <= 2 # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(deltaDIC <= 2, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

fwrite(
  combined_catalogue, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_combined_catalogue.csv"
)
