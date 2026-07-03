library(numDeriv)
library(Rfits)
library(data.table)
library(Highlander)
library(magicaxis)
library(foreach)
library(pracma)
library(rhdf5)

# likely photoz wrong objects 

eazy_dir = "~/Documents/Epochs/data/eazy/"
fname = "smiles"
fname_append = "_agn_blue_sfhz"

phot_cat = Rfits_read(
  paste0(main_dir, "/data/catalogues/", fname, "_profound_multiband_cat.fits")
)
filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))

phot_cat32 = data.frame(phot_cat$Phot32)
df_flux = phot_cat32[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = phot_cat32[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

eazy = Rfits_read_table(
  paste0(eazy_dir, fname, fname_append, "/", fname, "_zphot", fname_append, ".out.zout.fits")
)
eazy_stats = data.frame(fread(paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv")))

maghist(
  eazy$z_ml, log = 'y', breaks = 30
)

select_idx = eazy_stats$zbestPzIntegral >= 0.6 & phot_cat$ProFound$flag == 0L

maghist(
  eazy$z_ml[select_idx], log = 'y', breaks = 30
)

idx = which(eazy$z_ml >= 10 & select_idx == 1L)

for(kk in idx){
  
  png(
    paste0(
      "~/Documents/Epochs/plots/high_z10_photoz/gal_",kk, ".png"
    ), width = 12, height = 6, units = "in", res = 100
  )
  par(mfrow = c(1,2), oma = rep(1.5, 4), mar = rep(1.5,4))
  
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
    ff,
    ylo = unlist(df_flux_err[kk, ])
  )
  
  xy_coords = Rwcs_s2p(
    phot_cat$ProFound$RAmax[kk],
    phot_cat$ProFound$Decmax[kk],
    keyvalues = stack$image$keyvalues
  )
  
  Rwcs_imageRGB(
    R = Rimage[
      phot_cat$ProFound$RAmax[kk],
      phot_cat$ProFound$Decmax[kk],
      box = 100,
      type = "coord"
    ],
    G = Gimage[
      phot_cat$ProFound$RAmax[kk],
      phot_cat$ProFound$Decmax[kk],
      box = 100,
      type = "coord"
    ],
    B = Bimage[
      phot_cat$ProFound$RAmax[kk],
      phot_cat$ProFound$Decmax[kk],
      box = 100,
      type = "coord"
    ]
  )
  legend(
    x = "topleft",
    legend = c(
      paste0("ID = ", kk),
      paste0("z = ", round(eazy$z_ml[kk],3))
    )
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
    ]$imDat, 
    add = TRUE
  )
  dev.off()
}
