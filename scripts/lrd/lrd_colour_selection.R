library(data.table)
library(magicaxis)
library(celestial)
library(Rfits)
library(Rwcs)

fname = "miri" # all data, smiles, midis, etc

phot_cat = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", fname, "_profound_multiband_cat.fits"), 
  ext = 2
))
phot_stats = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", fname, "_profound_multiband_cat.fits"), 
  ext = "PhotStats"
))

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

raw_cat = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", fname, "_profound_multiband_cat.fits"), 
  ext = "RawTot"
))
raw_cat16 = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", fname, "_profound_multiband_cat.fits"), 
  ext = "Raw16"
))
raw_cat32 = data.frame(Rfits_read_table(
  paste0("~/Documents/JWST-MIRI-AGN/data/catalogues/", fname, "_profound_multiband_cat.fits"), 
  ext = "Raw32"
))

eazy = Rfits_read_table(
  "~/Documents/JWST-MIRI-AGN/data/eazy/miri_agn_blue_sfhz_col/miri_zphot_agn_blue_sfhz_col.out.zout.fits"
)

filters_info = data.frame(fread(paste0("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0("~/Documents/JWST-MIRI-AGN/data/throughputs/filt_funcs.rds"))

df_flux = raw_cat[, paste0(filters_info$filters, "_fluxt")]
df_flux_err = raw_cat[, paste0(filters_info$filters, "_scaled_fluxt_err")]

# epochs2_goodss_bd_cat_objects = Rfits_read_table(
#   "~/Documents/JWST-MIRI-AGN/data/catalogues/epochs2_goodss_bd_cat.fits",
#   ext = 2
# )
# epochs2_goodss_bd_cat = Rfits_read_table(
#   "~/Documents/JWST-MIRI-AGN/data/catalogues/epochs2_goodss_bd_cat.fits",
#   ext = "SELECTION"
# )
# 
# epochs2_match_idx = coordmatch(
#   coordref = epochs2_goodss_bd_cat_objects[, c("ALPHA_J2000", "DELTA_J2000")],
#   coordcompare = phot_cat[, c("RAmax", "Decmax")],
#   rad = 0.3
# )
# epochs2_goodss_bd_cat_match = epochs2_goodss_bd_cat[epochs2_match_idx$bestmatch$refID, ]
# beta_powerlaw_df_match = beta_powerlaw_df[epochs2_match_idx$bestmatch$compareID, ]

get_flux = function(df_flux, band){
  fflux = df_flux[grepl(band, names(df_flux))]
  fflux[fflux != -999 & fflux <= 0] = phot_stats[paste0(band, "_1sig_local_depth")][fflux != -999 & fflux <= 0]
  return(fflux)
}

F115W_flux = unname(unlist(get_flux(df_flux, "F115W")))
F150W_flux = unname(unlist(get_flux(df_flux, "F150W")))
F200W_flux = unname(unlist(get_flux(df_flux, "F200W")))
F277W_flux = unname(unlist(get_flux(df_flux, "F277W")))
F444W_flux = unname(unlist(get_flux(df_flux, "F444W")))
F560W_flux = unname(unlist(get_flux(df_flux, "F560W")))
F770W_flux = unname(unlist(get_flux(df_flux, "F770W")))

## similar criteria used by Rinaldi+26
lrd_idx = (
  -2.5*log10(F150W_flux/F200W_flux) <= 1.0 & 
  -2.5*log10(F277W_flux/F444W_flux) >= 0.5 & 
  
  -2.5*log10(F115W_flux/F200W_flux) > -0.5 & # brown dwarf remove
  -2.5*log10(F115W_flux/F277W_flux) > -0.5 & 
  -2.5*log10(F277W_flux/F444W_flux) < 0.4 * -2.5*log10(F115W_flux/F277W_flux) + 1.9 &
  
  raw_cat32$F444W_flux_app_0.32/raw_cat16$F444W_flux_app_0.16 <= 1.7 & 
  (raw_cat32$F444W_flux_app_0.32 != -999 | raw_cat16$F444W_flux_app_0.16 != -999) & 
    
  F150W_flux != -999 &
  F200W_flux != -999 &
  F277W_flux != -999 &
  F444W_flux != -999 &
    
  phot_cat$flag == 0L & 
  phot_cat$flag == 0 & 
  eazy$z_ml > 0 &
  prospect_cat$zQ50 >= 1 & 
  prospect_cat$zQ50 < 10
)

for(kk in which(lrd_idx)){
  png(
    paste0("/Users/jordan/Documents/JWST-MIRI-AGN/plots/lrd_colour_select/fit_gal_",kk, ".png"),
    width = 7, height = 7, units = "in", res = 240
  )
  par(mfrow = c(1,1), oma = rep(1.5, 4), mar = rep(1.5,4))
  
  stack = Rfits_read(
    paste0(
      "~/mnt/morgan/raid/data/JP/ProFound/Detects/", 
      phot_cat$FIELD[kk], "/", phot_cat$FIELD[kk], "/", 
      phot_cat$FIELD[kk], "_", phot_cat$FIELD[kk], "_", phot_cat$FIELD[kk], "_profound_stack.fits"
    )
  )
  
  dir_path <- paste0(
    "~/mnt/morgan/raid/data/JP/ProFound/Data/", phot_cat$FIELD[kk], "/", phot_cat$FIELD[kk], "/"
  )
  
  Rimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F444W*.fits")),
    full.names = TRUE
  )[1])
  Gimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F277W*.fits")),
    full.names = TRUE
  )[1])
  Bimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F150W*.fits")),
    full.names = TRUE
  )[1])
  
  R = Rimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
    box = 100,
    type = "coord"
  ]
  G = Gimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
    box = 100,
    type = "coord"
  ]
  B = Bimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
    box = 100,
    type = "coord"
  ]
  
  R$imDat[is.na(R$imDat)] = 0
  G$imDat[is.na(G$imDat)] = 0
  B$imDat[is.na(B$imDat)] = 0
  
  Rwcs_imageRGB(
    R = R,
    G = G,
    B = B,
    grid.lty = 0,
    coord.type = "deg",
    com.col = NULL
  )
  xy_coords = Rwcs_s2p(
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
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
    "F277W",
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
  dev.off()
}

fwrite(
  list(which(lrd_idx)),
  file = "~/Documents/JWST-MIRI-AGN/data/prospect/LRD_colour_selection_IDs.txt", 
)


for(i in which(lrd_idx)){
  file.copy(
    paste0("~/Documents/JWST-MIRI-AGN/plots/fitted_candidates/gal_", i, ".png"),
    paste0("~/Documents/JWST-MIRI-AGN/plots/lrd_colour_select//gal_", i, ".png"),
    overwrite = TRUE,
    recursive = TRUE
  )
}

for(kk in which(lrd_idx)){
  png(
    paste0(
      "~/Documents/JWST-MIRI-AGN/plots/lrd_colour_select/fit_gal_",kk, ".png"
    ), width = 18, height = 6, units = "in", res = 240
  )
  par(mfrow = c(1,3), oma = rep(1.5, 4), mar = rep(1.5,4))
  
  stack = Rfits_read(
    paste0(
      "~/mnt/morgan/raid/data/JP/ProFound/Detects/", 
      phot_cat$FIELD[kk], "/", phot_cat$FIELD[kk], "/", 
      phot_cat$FIELD[kk], "_", phot_cat$FIELD[kk], "_", phot_cat$FIELD[kk], "_profound_stack.fits"
    )
  )
  
  dir_path <- paste0(
    "~/mnt/morgan/raid/data/JP/ProFound/Data/", phot_cat$FIELD[kk], "/", phot_cat$FIELD[kk], "/"
  )
  
  Rimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F444W*.fits")),
    full.names = TRUE
  )[1])
  Gimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F277W*.fits")),
    full.names = TRUE
  )[1])
  Bimage <- Rfits_point(list.files(
    dir_path,
    pattern    = glob2rx(paste0("*F150W*.fits")),
    full.names = TRUE
  )[1])
  
  R = Rimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
    box = 100,
    type = "coord"
  ]
  G = Gimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
    box = 100,
    type = "coord"
  ]
  B = Bimage[
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
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
    phot_cat$RAmax[kk],
    phot_cat$Decmax[kk],
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
    "F277W",
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
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy/galID_", kk, ".rds")
  )
  sed_stars = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy/galID_", kk, ".rds")
  )
  
  sed_agn_samp = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy/galID_", kk, ".rds")
  )
  sed_stars_samp = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_norm_eazy/galID_", kk, ".rds")
  )
  sed_mbb_samp = readRDS(
    paste0("~/Documents/JWST-MIRI-AGN/data/prospect_seds/photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd/galID_", kk, ".rds")
  )
  
  ff = unlist(df_flux[kk, ])
  ylims = range(ff[ff>0], na.rm = TRUE)
  magplot(
    filters_info$cenwave,
    ff,
    ylim = c(ylims[1]*0.5, ylims[2]*2),
    xlim = c(1000, 1e6),
    log = "xy"
  )
  magerr(
    filters_info$cenwave,
    unlist(df_flux[kk, ]),
    ylo = unlist(df_flux_err[kk, ])
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
    sed_agn_samp$obswave,
    sed_agn_samp$fluxQ50,
    ylo = sed_agn_samp$fluxQ50-sed_agn_samp$fluxQ16,
    yhi = sed_agn_samp$fluxQ84-sed_agn_samp$fluxQ50,
    col = adjustcolor("red",0.3),
    poly = TRUE
  )
  magerr(
    sed_stars_samp$obswave,
    sed_stars_samp$fluxQ50,
    ylo = sed_stars_samp$fluxQ50-sed_stars_samp$fluxQ16,
    yhi = sed_stars_samp$fluxQ84-sed_stars_samp$fluxQ50,
    col = adjustcolor("blue",0.3),
    poly = TRUE
    
  )
  lines(
    sed_mbb_samp$obswave,
    sed_mbb_samp$fluxQ50,
    col = "purple"
  )
  points(
    filters_info$cenwave,
    sed_mbb_samp$phot[paste0(filters_info$filters, "_fit_fluxQ50")],
    pch = 16
  )
  
  legend(
    x = "topleft",
    legend = c(
      paste0("DIC Stellar = ", round(sed_stars$galInfo$DIC, 3)), 
      paste0("DIC Stellar+AGN= ", round(sed_agn$galInfo$DIC, 3)),
      # paste0("DIC Stellar+LRD LIU= ", round(sed_liu$galInfo$DIC, 3)),
      
      paste0("LP Stellar = ", round(sed_stars$astro$LPQ50, 3)), 
      paste0("LP Stellar+AGN= ", round(sed_agn$astro$LPQ50, 3))
      # paste0("LP Stellar+LRD LIU= ", round(sed_liu$astro$LPQ50, 3))
    ), 
    lty = c(NA, NA, 1, 1, 1),
    col = c(NA, NA, "green", "blue", "red")
  )
  legend(
    x = "bottomright",
    legend = c(
      paste0("z = ", round(combined_cat$zQ50[kk], 3)), 
      paste0("ID = ", round(combined_cat$ID[kk], 3)), 
      "Stellar+AGN SED fit",
      "Stellar SED fit", 
      paste0("Better fit = ", combined_cat$FitType[kk])
    ), 
    lty = c(NA, NA, 1, 1, 1),
    col = c(NA, NA, "green", "blue", "red")
  )
  dev.off()
}