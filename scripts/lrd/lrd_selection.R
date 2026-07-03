library(data.table)
library(magicaxis)
library(celestial)
library(Rfits)
library(Rwcs)

phot_cat = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = 2
))
phot_cat16 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Phot16"
))
phot_cat32 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Phot32"
))
phot_cat60 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Phot60"
))
beta_powerlaw_df = data.frame(fread(
  "~/Documents/Epochs/data/prospect/power_laws_astro_cat.csv"
))
phot_stats = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "PhotStats"
))
raw_cat16 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Raw16"
))
raw_cat32 = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", 
  ext = "Raw32"
))

filters_info = data.frame(fread(paste0("~/Documents/Epochs/data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0("~/Documents/Epochs/data/throughputs/filt_funcs.rds"))

df_flux = phot_cat32[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = phot_cat32[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

# epochs2_goodss_bd_cat_objects = Rfits_read_table(
#   "~/Documents/Epochs/data/catalogues/epochs2_goodss_bd_cat.fits", 
#   ext = 2
# )
# epochs2_goodss_bd_cat = Rfits_read_table(
#   "~/Documents/Epochs/data/catalogues/epochs2_goodss_bd_cat.fits", 
#   ext = "SELECTION"
# )

# epochs2_match_idx = coordmatch(
#   coordref = epochs2_goodss_bd_cat_objects[, c("ALPHA_J2000", "DELTA_J2000")],
#   coordcompare = beta_powerlaw_df[, c("RA", "Dec")], 
#   rad = 0.3
# )
# epochs2_goodss_bd_cat_match = epochs2_goodss_bd_cat[epochs2_match_idx$bestmatch$refID, ]
# beta_powerlaw_df_match = beta_powerlaw_df[epochs2_match_idx$bestmatch$compareID, ]

magplot(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaUVQ50,
  xlim = c(0, 20),
  # ylim = c(-5, 5)
)
magerr(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaUVQ50,
  ylo = beta_powerlaw_df$betaUV_lo,
  yhi = beta_powerlaw_df$betaUV_hi
)
magplot(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaOptQ50,
  xlim = c(0, 20),
  # ylim = c(-5, 5)
)
magerr(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaOptQ50,
  ylo = beta_powerlaw_df$betaOpt_lo,
  yhi = beta_powerlaw_df$betaOpt_hi
)
# 

lrd_idx = which(beta_powerlaw_df$betaUVQ50 <= -0.37 & beta_powerlaw_df$betaOptQ50 > 0 & beta_powerlaw_df$betaUVQ50 > -2.8)
# pdf("~/Desktop/beta_powerlaws.pdf", width = 8, height = 5)
magplot(
  beta_powerlaw_df$betaUVQ50,
  beta_powerlaw_df$betaOptQ50,
  xlim = c(-5, 2),
  ylim = c(-4, 4),
  xlab = "Beta UV slope",
  ylab = "Beta Optical slope"
)
points(
  beta_powerlaw_df$betaUVQ50[lrd_idx],
  beta_powerlaw_df$betaOptQ50[lrd_idx],
  col = "red",
  pch = 16
)
# legend(
#   x = "topright",
#   legend = c(
#     paste0("Better fit with AGN models"), 
#     paste0("Kocevski+25 LRD criterion")
#   ),
#   pch = c(1, 16),
#   col = c("black", "red")
# )
# dev.off()
# 

Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F444W_smiles_warp_smiles.fits")
Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F200W_smiles_warp_smiles.fits")
Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F150W_smiles_warp_smiles.fits")
stack = Rfits_read("~/Documents/Epochs/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")

for(kk in 1:length(lrd_idx)){
  pdf(
    paste0(
      "~/Documents/Epochs/plots/lrd_candidates/gal_",beta_powerlaw_df$ID[lrd_idx][kk], ".pdf"
    ), width = 18, height = 6
  )
  par(mfrow = c(1,3), oma = rep(1.5, 4), mar = rep(1.5,4))
  
  Rwcs_imageRGB(
    R = Rimage[
      beta_powerlaw_df$RA[lrd_idx][kk],
      beta_powerlaw_df$Dec[lrd_idx][kk],
      box = 100,
      type = "coord"
    ],
    G = Gimage[
      beta_powerlaw_df$RA[lrd_idx][kk],
      beta_powerlaw_df$Dec[lrd_idx][kk],
      box = 100,
      type = "coord"
    ],
    B = Bimage[
      beta_powerlaw_df$RA[lrd_idx][kk],
      beta_powerlaw_df$Dec[lrd_idx][kk],
      box = 100,
      type = "coord"
    ]
  )
  xy_coords = Rwcs_s2p(
    beta_powerlaw_df$RA[lrd_idx][kk],
    beta_powerlaw_df$Dec[lrd_idx][kk],
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
    stack$image[
      xy_coords[1,1],
      xy_coords[1,2],
      box = 100
    ],
    stack$segim[
      xy_coords[1,1],
      xy_coords[1,2],
      box = 100
    ]$imDat, 
    add = FALSE
  )
  sed_pl = readRDS(
    paste0("~/Documents/Epochs/data/prospect/power_laws/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
  )
  sed_agn = readRDS(
    paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
  )
  sed_stars = readRDS(
    paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
  )

  magplot(
    filters_info$cenwave,
    unlist(df_flux[beta_powerlaw_df$ID[lrd_idx][kk], ]),
    ylim = c(1e-10, 1e-5),
    xlim = c(1000, 1e6),
    log = "xy"
  )
  magerr(
    filters_info$cenwave,
    unlist(df_flux[beta_powerlaw_df$ID[lrd_idx][kk], ]),
    ylo = unlist(df_flux_err[beta_powerlaw_df$ID[lrd_idx][kk], ])
  )
  lines(
    sed_pl$bestfit$FinalFluxUV$wave,
    sed_pl$bestfit$FinalFluxUV$flux,
    col = "green"
  )
  lines(
    sed_pl$bestfit$FinalFluxOpt$wave,
    sed_pl$bestfit$FinalFluxOpt$flux,
    col = "green"
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
  # lines(
  #   sed_lrd$bestfit$SEDout$FinalFlux$wave,
  #   sed_lrd$bestfit$SEDout$FinalFlux$flux,
  #   col = "purple"
  # )
  legend(
    x = "bottomright",
    legend = c(
      paste0("z = ", round(beta_powerlaw_df$zQ50[lrd_idx][kk], 3)), 
      paste0("ID = ", round(beta_powerlaw_df$ID[lrd_idx][kk], 3)), 
      # paste0("ID = ", round(phot_cat$flag_hotpixel[lrd_idx][kk], 3)),
      "UV-Optical power laws",
      "Stellar+AGN SED fit",
      "Stellar SED fit"
    ), 
    lty = c(NA, NA, 1, 1, 1),
    col = c(NA, NA, "green", "blue", "red")
  )
  dev.off()
}


kocevski25 = read.table("~/Documents/Epochs/data/literature/AGN-Catalogues/kocevski25_LRDs.txt", skip = 35)
match_kocevski = coordmatch(
  coordref = kocevski25[,3:4],
  coordcompare = beta_powerlaw_df[, c("RA", "Dec")],
  rad = 0.3
)

magplot(
  kocevski25$V5[match_kocevski$bestmatch$refID],
  beta_powerlaw_df$zQ50[match_kocevski$bestmatch$compareID]
)
abline(0,1)
# 
# kk = 4
# RA = kocevski25$V3[matchbox_faggot$bestmatch$refID][kk]
# Dec = kocevski25$V4[matchbox_faggot$bestmatch$refID][kk]
# Rwcs_imageRGB(
#   R = Rimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ],
#   G = Gimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ],
#   B = Bimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ]
# )
# 
# foo = readRDS(
#   paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", combined_catalogue$ID[matchbox_faggot$bestmatch$compareID][kk], ".rds")
# )
# magplot(
#   filters_info$cenwave,
#   unlist(df_flux[combined_catalogue$ID[matchbox_faggot$bestmatch$compareID][kk], ]), 
#   ylim = c(1e-8, 1e-6), 
#   log = "xy"
# )
# lines(
#   foo$bestfit$FinalFluxUV$wave,
#   foo$bestfit$FinalFluxUV$flux
# )
# lines(
#   foo$bestfit$FinalFluxOpt$wave,
#   foo$bestfit$FinalFluxOpt$flux
# )
# lines(
#   foo$bestfit$SEDout$FinalFlux
# )
# plot(foo$bestfit, ylim = c(1e-9, 1e-5))


magplot(
  phot_cat$Phot16$F444W_flux_app_0.16*phot_cat$AperStats$F444W_ap_correction_app_0.16,
  phot_cat$Phot60$F444W_flux_app_0.6*phot_cat$AperStats$F444W_ap_correction_app_0.6, 
  xlim = c(1e-11, 2),
  ylim = c(1e-11, 2),
  log = 'xy'
)
