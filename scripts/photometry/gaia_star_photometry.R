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

# main_dir = "~/Documents/Epochs/"
main_dir = "~/mnt/morgan/raid/"
fname = "smiles"

star_mask = readRDS(
  paste0(main_dir, "data/stars/", fname, "/", fname, "_GAIA_star_mask.rds")
)

segim_stack = Rfits_read(
  filename = paste0(main_dir, "data/JP/ProFound/Data/", fname, "/", fname, )
)

im_list_stubs = list.files(
  paste0(main_dir, "data/JP/ProFound/Data/", fname, "/", fname, "/"), 
  recursive = TRUE, 
  full.names = TRUE
)
filt_names = toupper(unique( regmatches(im_list_stubs, regexpr("F\\d+[WML]P?", im_list_stubs, ignore.case = TRUE)) ))

print(im_list_stub)

im_list = Rfits_make_list(
  filelist = im_list_stubs, 
  extlist = 1
)
names(im_list) = filt_names

magzero_list = unname(sapply(im_list, function(x){
  -2.5*log10(x$keyvalues$PIXAR_SR) - 6.10
}))

star_cat = foreach(i = 1:length(im_list), .combine = cbind) %do% {
  profound_aperture_phot = profoundAperPhot(
    image = im_list[[i]][,]$imDat,
    segim = star_mask$stars,
    app_diam = c(0.16,0.32)/pixscale(im_list[[i]]) * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]),
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

gaia_trim <- fread(
  file = list.files(
    path       = paste0(main_dir, "/data/stars/", fname, "/"),
    pattern    = paste0(fname, "_trim_.*\\.csv"),
    full.names = TRUE
  )
)
star_cat$RA = gaia_trim$ra_fix
star_cat$DEC = gaia_trim$dec_fix

star_cat = star_cat[, c("RA", "DEC", setdiff(names(star_cat), c("RA", "DEC")))]
fwrite(
  star_cat,
  paste0(main_dir, "/data/stars/", fname, "GAIA_MIRI_photometry.csv")
)
