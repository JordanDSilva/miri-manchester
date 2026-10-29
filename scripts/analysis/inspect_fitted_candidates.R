library(data.table)
library(Rfits)
library(Rwcs)
library(rhdf5)
library(magicaxis)
library(celestial)

phot_cat = Rfits_read("~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)
combined_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
eazy = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
flux_cat = data.frame(Rfits_read_table(
  "~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", 
  ext = "PhotTot"
))

filters_info = data.frame(fread(paste0("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/throughputs/filt_funcs.rds"))

df_flux = flux_cat[, paste0(filters_info$filters, "_fluxt")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_fluxt_err")]

smiles_Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F444W_smiles_warp_smiles.fits")
smiles_Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F200W_smiles_warp_smiles.fits")
smiles_Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F150W_smiles_warp_smiles.fits")
smiles_stack = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")

midis_Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/midis/midis/JADES-DR3-GS_F444W_midis_warp_midis.fits")
midis_Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/midis/midis/JADES-DR3-GS_F200W_midis_warp_midis.fits")
midis_Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/midis/midis/JADES-DR3-GS_F150W_midis_warp_midis.fits")
midis_stack = Rfits_read("~/Documents/JWST-MIRI-AGN/data/images/midis/midis_midis_midis_profound_stack.fits")

idx = which(!is.na(combined_cat$zQ50))

for(kk in 1:length(idx)){ 
  
  ID = combined_cat$CATID[idx[kk]]
  field = phot_cat$ProFound$FIELD[idx[kk]]
  
  png(
    paste0(
      "~/Documents/JWST-MIRI-AGN/plots/fitted_candidates/gal_", ID, ".png"
    ), width = 18, height = 6, units = "in", res = 240
  )
  par(mfrow = c(1,3), oma = rep(1.5, 4), mar = rep(1.5,4))
  
  if(field == "smiles"){
    Rimage = smiles_Rimage
    Gimage = smiles_Gimage
    Bimage = smiles_Bimage
    stack = smiles_stack
  }else if (field == "midis"){
    Rimage = midis_Rimage
    Gimage = midis_Gimage
    Bimage = midis_Bimage
    stack = midis_stack
  }
  
  R = Rimage[
    combined_cat$RA[idx[kk]],
    combined_cat$Dec[idx[kk]],
    box = 100,
    type = "coord"
  ]
  G = Gimage[
    combined_cat$RA[idx[kk]],
    combined_cat$Dec[idx[kk]],
    box = 100,
    type = "coord"
  ]
  B = Bimage[
    combined_cat$RA[idx[kk]],
    combined_cat$Dec[idx[kk]],
    box = 100,
    type = "coord"
  ]
  
  R$imDat[is.na(R$imDat)] = 0
  G$imDat[is.na(G$imDat)] = 0
  B$imDat[is.na(B$imDat)] = 0
  
  Rwcs_imageRGB(
    R = R,
    G = G,
    B = B
  )
  xy_coords = Rwcs_s2p(
    combined_cat$RA[idx[kk]],
    combined_cat$Dec[idx[kk]],
    keyvalues = Rimage$keyvalues
  )
  ProFound::profoundSegimPlot(
    image = matrix(0, 100, 100),
    stack$segim[
      xy_coords[1,1],
      xy_coords[1,2],
      box = 100
    ]$imDat, 
    add = TRUE
  )
  text(
    x = 30,
    y = 20,
    "F150W",
    col = "cyan", 
    cex = 1.0
  )
  text(
    x = 50,
    y = 20,
    "F200W",
    col = "green", 
    cex = 1.0
  )
  text(
    x = 70,
    y = 20,
    "F444W",
    col = "red", 
    cex = 1.0
  )
  
  ProFound::profoundSegimPlot(
    image = stack$image[
      xy_coords[1,1],
      xy_coords[1,2],
      box = 100
    ]$imDat,
    stack$segim[
      xy_coords[1,1],
      xy_coords[1,2],
      box = 100
    ]$imDat
  )
  
  sed_agn = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", ID, ".rds")
  )
  sed_stars = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy/galID_", ID, ".rds")
  )
  
  fit_agn = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy/galID_", ID, ".rds")
  )
  fit_stars = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_norm_eazy/galID_", ID, ".rds")
  )
  
  ff = unlist(df_flux[idx[kk], ])
  ylims = range(ff[ff>0], na.rm = TRUE)
  magplot(
    NA,
    ylim = c(ylims[1]*0.5, ylims[2]*2),
    xlim = c(1000, 1e6),
    log = "xy"
  )
  lines(
    sed_stars$bestfit$SEDout$FinalFlux$wave,
    sed_stars$bestfit$SEDout$FinalFlux$flux,
    col = "blue"
  )
  lines(
    sed_agn$bestfit$SEDout$FinalFlux$wave,
    sed_agn$bestfit$SEDout$FinalFlux$flux,
    col = "red"
  )
  magerr(
    fit_agn$obswave,
    fit_agn$fluxQ50,
    ylo = fit_agn$fluxQ50-fit_agn$fluxQ16,
    yhi = fit_agn$fluxQ84-fit_agn$fluxQ50,
    col = adjustcolor("red",0.3),
    border = NA,
    poly = TRUE
  )
  magerr(
    fit_stars$obswave,
    fit_stars$fluxQ50,
    ylo = fit_stars$fluxQ50-fit_stars$fluxQ16,
    yhi = fit_stars$fluxQ84-fit_stars$fluxQ50,
    col = adjustcolor("blue",0.3),
    border = NA,
    poly = TRUE
  )
  points(
    filters_info$cenwave,
    ff,
    pch = 16,
    col = adjustcolor("black", 0.5)
  )
  magerr(
    filters_info$cenwave,
    unlist(df_flux[idx[kk], ]),
    ylo = unlist(df_flux_err[idx[kk], ])
  )
  legend(
    x = "topleft",
    legend = c(
      paste0("DIC Stellar = ", round(sed_stars$galInfo$DIC, 3)), 
      paste0("DIC Stellar+AGN= ", round(sed_agn$galInfo$DIC, 3)),
      
      paste0("LP Stellar = ", round(sed_stars$astro$LPQ50, 3)), 
      paste0("LP Stellar+AGN= ", round(sed_agn$astro$LPQ50, 3))
    )
  )
  legend(
    x = "bottomright",
    legend = c(
      paste0("z = ", round(combined_cat$zQ50[idx[kk]], 3)), 
      paste0("ID = ", round(combined_cat$ID[idx[kk]], 3)), 
      "Stellar+AGN SED fit",
      "Stellar SED fit", 
      paste0("Better fit = ", combined_cat$FitType[idx[kk]])
    ), 
    lty = c(NA, NA, 1, 1, 1),
    col = c(NA, NA, "blue", "red", NA)
  )
  dev.off()
}
