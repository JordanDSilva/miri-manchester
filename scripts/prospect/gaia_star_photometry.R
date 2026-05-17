library(Rwcs)
library(Rfits)
library(ProPane)
library(ProFound)
library(celestial)
library(data.table)
library(foreach)
library(magicaxis)
library(scales)
library(stringr)

main_dir = "~/Documents/Epochs/"

segim_stack = Rfits_read(
    filename = paste0(main_dir, "data/images/smiles_F560W_F770W_invar_stack_propane_segim.fits")
)

filt_names = read.table(
  file = paste0(main_dir, "data/throughputs/wavelengths.txt")
)[12:19,2]

F560W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f560w_v1.0_drz.fits")
)
F770W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f770w_v1.0_drz.fits")
)
F1000W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f1000w_v1.0_drz.fits")
)
F1280W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f1280w_v1.0_drz.fits")
)
F1500W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f1500w_v1.0_drz.fits")
)
F1800W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f1800w_v1.0_drz.fits")
)
F2100W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f2100w_v1.0_drz.fits")
)
F2550W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f2550w_v1.0_drz.fits")
)

im_list = list(
  F560W, F770W, F1000W, F1280W, F1500W, F1800W, F2100W, F2550W
)
names(im_list) = filt_names

magzero_list = unname(sapply(im_list, function(x){
  -2.5*log10(x$SCI$keyvalues$PIXAR_SR) - 6.10
}))

star_cat = foreach(i = 1:length(im_list), .combine = cbind) %do% {
  profound_aperture_phot = profoundAperPhot(
    image = im_list[[i]]$SCI[,]$imDat,
    segim = segim_stack$STARMASK[,]$imDat,
    app_diam = c(0.16,0.32)/pixscale(im_list[[i]]$SCI) * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Jansky", 
    depth = 5,
    correction = TRUE
  )
  df = data.frame(
    profound_aperture_phot$flux_app_1,
    profound_aperture_phot$flux_app_2
  )
  names(df) = c(
    paste0(filt_names[i], "_flux_aper16"),
    paste0(filt_names[i], "_flux_aper32")
  )
  return(df)
}

gaia_trim = fread(
  file = paste0(main_dir, "data/GAIA/GAIA_53.14977_-27.80586_trim.csv")
)
star_cat$RA = gaia_trim$ra_fix
star_cat$DEC = gaia_trim$dec_fix

star_cat = star_cat[, c("RA", "DEC", setdiff(names(star_cat), c("RA", "DEC")))]
fwrite(
  star_cat,
  paste0(main_dir, "/data/GAIA/GAIA_53.14977_-27.80586_MIRI_photometry.csv")
)
