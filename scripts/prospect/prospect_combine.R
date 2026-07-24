library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))

phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)

photoZ_bc03_snorm_fritz_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

common_names = intersect(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))
common_names = common_names[common_names %notin% c("ID", "RA", "Dec")]
uncommon_names = setdiff(names(photoZ_bc03_snorm_fritz_norm_eazy), names(photoZ_bc03_snorm_norm_eazy))

combined_catalogue = photoZ_bc03_snorm_fritz_norm_eazy # start with Stellar+AGN

AGN_UV_ratio = 10^(-0.4 * (photoZ_bc03_snorm_fritz_norm_eazy$UV1500AGNQ50 - photoZ_bc03_snorm_fritz_norm_eazy$UV1500Q50) )
DIC_thresh = 0.0
deltaDIC = photoZ_bc03_snorm_norm_eazy$DIC - photoZ_bc03_snorm_fritz_norm_eazy$DIC

replace_AGN_idx = !(deltaDIC >= DIC_thresh & (AGN_UV_ratio >= 0.1 | photoZ_bc03_snorm_fritz_norm_eazy$AGN_Dale_frac_Q50 >= 0.1))
# replace_AGN_idx = !(deltaDIC >= DIC_thresh)
combined_catalogue[which(replace_AGN_idx), common_names] = photoZ_bc03_snorm_norm_eazy[which(replace_AGN_idx), common_names] # now set better fit to Stellar
combined_catalogue[which(replace_AGN_idx), uncommon_names] = NA # set the AGN crap to NA

combined_catalogue$DeltaDIC = deltaDIC
combined_catalogue$AGNFitLowerLimit =  replace_AGN_idx # DIC cant favour AGN, hence lower limit on AGN contribution to the fit
combined_catalogue$FitType = ifelse(replace_AGN_idx, "Stellar", "StellarAGN") # DIC cant favour AGN, hence lower limit on AGN contribution to the fit

## remove the photometric flags
robust_catalogue = combined_catalogue
robust_catalogue[which(phot_cat$ProFound$flag == 1 | combined_catalogue$StellarMassQ50 > 10^(12)), c(common_names, uncommon_names, "DeltaDIC", "AGNFitLowerLimit", "FitType")] = NA

fwrite(
  robust_catalogue, 
  file = "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
)

h5file = "~/Documents/Epochs/data/all_data.h5"

# robust_photom = combined_photom
# robust_photom[phot_cat$ProFound$flag == 1 ,] = NA

# aperture_correction_factor_temp = phot_cat$RawTot[, c("F560W_fluxt", "F770W_fluxt")]/phot_cat$Raw32[, c("F560W_flux_app_0.32", "F770W_flux_app_0.32")]
# aperture_correction_factor_temp$F560W_fluxt[is.infinite(aperture_correction_factor_temp$F560W_fluxt) | aperture_correction_factor_temp$F560W_fluxt < 0] = NA
# aperture_correction_factor_temp$F770W_fluxt[is.infinite(aperture_correction_factor_temp$F770W_fluxt) | aperture_correction_factor_temp$F770W_fluxt < 0] = NA
# aperture_correction_factor = pmax(1, rowMins(
#   as.matrix(aperture_correction_factor_temp)
# ))
# 
# Ngas_col_names <- grep("NgasMass", names(robust_catalogue))
# robust_catalogue[,Ngas_col_names] <- lapply(robust_catalogue[,Ngas_col_names], function(x) {
#   as.numeric(as.character(x))
# })
# 
# robust_catalogue[, grep("StellarMass", names(robust_catalogue))] = robust_catalogue[, grep("StellarMass", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("SFR", names(robust_catalogue))] = robust_catalogue[, grep("SFR", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("AGNlum", names(robust_catalogue))] = robust_catalogue[, grep("AGNlum", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("DustMass", names(robust_catalogue))] = robust_catalogue[, grep("DustMass", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("DustLum", names(robust_catalogue))] = robust_catalogue[, grep("DustLum", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("NgasMass", names(robust_catalogue))] = robust_catalogue[, grep("NgasMass", names(robust_catalogue))] * aperture_correction_factor
# robust_catalogue[, grep("mSFR", names(robust_catalogue))] = robust_catalogue[, grep("mSFR", names(robust_catalogue))] * aperture_correction_factor
# 
# # have to be more careful with logged quantities
# robust_catalogue[, grep("UV1500Best", names(robust_catalogue))] = robust_catalogue[, grep("UV1500Best", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500Q50", names(robust_catalogue))] = robust_catalogue[, grep("UV1500Q50", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500Q16", names(robust_catalogue))] = robust_catalogue[, grep("UV1500Q16", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500Q84", names(robust_catalogue))] = robust_catalogue[, grep("UV1500Q84", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500AGNBest", names(robust_catalogue))] = robust_catalogue[, grep("UV1500AGNBest", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500AGNQ50", names(robust_catalogue))] = robust_catalogue[, grep("UV1500AGNQ50", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500AGNQ16", names(robust_catalogue))] = robust_catalogue[, grep("UV1500AGNQ16", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# robust_catalogue[, grep("UV1500AGNQ84", names(robust_catalogue))] = robust_catalogue[, grep("UV1500AGNQ84", names(robust_catalogue))] - 2.5*log10(aperture_correction_factor)
# 
# robust_catalogue$apperture_correction = aperture_correction_factor