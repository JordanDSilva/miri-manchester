library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)
library(png)
library(dplyr)

robust_catalogue = data.frame(fread(file = "~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
robust_catalogue_z = robust_catalogue[which(robust_catalogue$zQ50 >= 1.07467287 & robust_catalogue$zQ50 < 9.98276084), ]

low_stellarmass_ids = robust_catalogue_z$CATID[
  robust_catalogue_z$StellarMassQ50 <= 10^8 & 
    robust_catalogue_z$FitType == "StellarAGN"
]

for(catid in low_stellarmass_ids){
  file.copy(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy/plot_", catid, ".pdf"),
    paste0("/Users/jordan/Documents/JWST-MIRI-AGN/rough-plots/low_mass_agns/stellarAGN/plot_", catid, ".pdf"),
    overwrite = TRUE
  )
  
  file.copy(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_norm_eazy/plot_", catid, ".pdf"),
    paste0("/Users/jordan/Documents/JWST-MIRI-AGN/rough-plots/low_mass_agns/stellar/plot_", catid, ".pdf"),
    overwrite = TRUE
  )
}

low_stellarmass_ids[10]
low_stellarmass_ids[16]
low_stellarmass_ids[66]
low_stellarmass_ids[107]
low_stellarmass_ids[109]
low_stellarmass_ids[112]


h5file = "~/Documents/JWST-MIRI-AGN/data/analysis_data.h5"

if (!file.exists(h5file)) {
  h5createFile(h5file)
}

h5createGroup(
  h5file, 
  "SEDs"
)
h5createGroup(
  h5file, 
  "SEDs/Stellar"
)
h5createGroup(
  h5file, 
  "SEDs/StellarAGN"
)
h5createGroup(
  h5file, 
  "SEDs/photom"
)

h5createGroup(
  h5file, 
  "SFHs"
)
h5createGroup(
  h5file, 
  "SFHs/Stellar"
)
h5createGroup(
  h5file, 
  "SFHs/StellarAGN"
)

h5write_replace = function(obj, file, name) {
  # create parent groups if needed (no-op/warning if they already exist)
  grp = dirname(name)
  suppressWarnings(try(h5createGroup(file, grp), silent = TRUE))
  
  # delete the old entry if present
  suppressWarnings(try(h5delete(file, name), silent = TRUE))
  
  h5write(obj, file = file, name = name)
}

save_sed = function(catid, h5file){
  
  fit_sed_data_stellarAGN = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", catid, ".rds"))
  fit_sed_stellarAGN = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy/galID_", catid, ".rds"))
  
  fit_sed_data_stellar = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy/galID_", catid, ".rds"))
  fit_sed_stellar = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_norm_eazy/galID_", catid, ".rds"))
  
  # save the SED to a h5 for plotting
  
  #observer frame
  stellarAGN_sed_df = data.frame(
    "obswave" = fit_sed_stellarAGN$obswave,
    
    "fluxQ50" = fit_sed_stellarAGN$fluxQ50,
    "fluxQ16" = fit_sed_stellarAGN$fluxQ16,
    "fluxQ84" = fit_sed_stellarAGN$fluxQ84,
    
    "AGNfluxQ50" = fit_sed_stellarAGN$AGNfluxQ50,
    "AGNfluxQ16" = fit_sed_stellarAGN$AGNfluxQ16,
    "AGNfluxQ84" = fit_sed_stellarAGN$AGNfluxQ84,
    
    "StarsUnAttenfluxQ50" = fit_sed_stellarAGN$StarsUnAttenfluxQ50,
    "StarsUnAttenfluxQ16" = fit_sed_stellarAGN$StarsUnAttenfluxQ16,
    "StarsUnAttenfluxQ84" = fit_sed_stellarAGN$StarsUnAttenfluxQ84,
    
    "StarsAttenfluxQ50" = fit_sed_stellarAGN$StarsAttenfluxQ50,
    "StarsAttenfluxQ16" = fit_sed_stellarAGN$StarsAttenfluxQ16,
    "StarsAttenfluxQ84" = fit_sed_stellarAGN$StarsAttenfluxQ84,
    
    "DustEmitfluxQ50" = fit_sed_stellarAGN$DustEmitfluxQ50,
    "DustEmitfluxQ16" = fit_sed_stellarAGN$DustEmitfluxQ16,
    "DustEmitfluxQ84" = fit_sed_stellarAGN$DustEmitfluxQ84
  )
  
  stellar_sed_df = data.frame(
    "obswave" = fit_sed_stellar$obswave,
    
    "fluxQ50" = fit_sed_stellar$fluxQ50,
    "fluxQ16" = fit_sed_stellar$fluxQ16,
    "fluxQ84" = fit_sed_stellar$fluxQ84,
    
    "AGNfluxQ50" = fit_sed_stellar$AGNfluxQ50,
    "AGNfluxQ16" = fit_sed_stellar$AGNfluxQ16,
    "AGNfluxQ84" = fit_sed_stellar$AGNfluxQ84,
    
    "StarsUnAttenfluxQ50" = fit_sed_stellar$StarsUnAttenfluxQ50,
    "StarsUnAttenfluxQ16" = fit_sed_stellar$StarsUnAttenfluxQ16,
    "StarsUnAttenfluxQ84" = fit_sed_stellar$StarsUnAttenfluxQ84,
    
    "StarsAttenfluxQ50" = fit_sed_stellar$StarsAttenfluxQ50,
    "StarsAttenfluxQ16" = fit_sed_stellar$StarsAttenfluxQ16,
    "StarsAttenfluxQ84" = fit_sed_stellar$StarsAttenfluxQ84,
    
    "DustEmitfluxQ50" = fit_sed_stellar$DustEmitfluxQ50,
    "DustEmitfluxQ16" = fit_sed_stellar$DustEmitfluxQ16,
    "DustEmitfluxQ84" = fit_sed_stellar$DustEmitfluxQ84
  )
  
  flux_data = fit_sed_data_stellar$bestfit$Data$flux
  
  stellarAGN_sfh_df = data.frame(
    "age"    = fit_sed_stellarAGN$ageQ50,
    
    "SFHQ50" = fit_sed_stellarAGN$SFHQ50,
    "SFHQ16" = fit_sed_stellarAGN$SFHQ16,
    "SFHQ84" = fit_sed_stellarAGN$SFHQ84,
    
    "ZHQ50"  = fit_sed_stellarAGN$ZHQ50,
    "ZHQ16"  = fit_sed_stellarAGN$ZHQ16,
    "ZHQ84"  = fit_sed_stellarAGN$ZHQ84
  )
  stellar_sfh_df = data.frame(
    "age"    = fit_sed_stellar$ageQ50,
    
    "SFHQ50" = fit_sed_stellar$SFHQ50,
    "SFHQ16" = fit_sed_stellar$SFHQ16,
    "SFHQ84" = fit_sed_stellar$SFHQ84,
    
    "ZHQ50"  = fit_sed_stellar$ZHQ50,
    "ZHQ16"  = fit_sed_stellar$ZHQ16,
    "ZHQ84"  = fit_sed_stellar$ZHQ84
  )
  
  h5write_replace(stellarAGN_sed_df, h5file, paste0("SEDs/StellarAGN/gal_", catid))
  h5write_replace(stellar_sed_df,    h5file, paste0("SEDs/Stellar/gal_",    catid))
  h5write_replace(flux_data,    h5file, paste0("SEDs/photom/gal_",    catid))
  
  h5write_replace(stellarAGN_sfh_df, h5file, paste0("SFHs/StellarAGN/gal_", catid))
  h5write_replace(stellar_sfh_df,    h5file, paste0("SFHs/Stellar/gal_",    catid))
}

save_sed(
  898, h5file
)

