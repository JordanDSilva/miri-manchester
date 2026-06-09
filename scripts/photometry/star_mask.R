library(Rfits)
library(Rwcs)
library(celestial)
library(ProFound)
library(data.table)
library(foreach)
library(magicaxis)

# star mask and star photometry

main_dir = "~/Documents/Epochs/"
# main_dir = "/"

fname = 'midis'
gaia_cat_stub = "GAIA_53.14977_-27.80586.csv"
frame_dir = paste0(main_dir, "data/images/", fname, "/")

gaia = data.frame(
  fread(
    file = paste0(main_dir, "data/stars/", fname, "/", gaia_cat_stub)
  )
) # created with JUMPROPE GAIA interface
im_frame_stubs = list.files(
  frame_dir, 
  recursive = TRUE,
  full.names = TRUE
)
F770W_stub = im_frame_stubs[grepl("F770W", im_frame_stubs, ignore.case = TRUE)]
im_frame = Rfits_read_image(
  F770W_stub
)
im_magzero = ifelse(
  is.null(im_frame$keyvalues$MAGZERO),
  -6.10 - 2.5*log10(im_frame$keyvalues$PIXAR_SR), 
  im_frame$keyvalues$MAGZERO
)

## Load in GAIA catalogue 
imdim = dim(im_frame)
im_ra_dec = Rwcs_p2s(
  x = imdim[1], 
  y = imdim[2], 
  keyvalues = im_frame$keyvalues
) #get extent of frame in RA and DEC coords

delta_time = 2023 - gaia$ref_epoch ## observed in December 2022?
gaia[, c("xpix", "ypix")] = Rwcs_s2p(
  RA = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time, ## Do proper motion correction
  Dec = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time, ## pmra[pmdec] = mas/yr
  keyvalues = im_frame$keyvalues
)
gaia[, c("xpix_nopm", "ypix_nopm")] = Rwcs_s2p(
  RA = gaia$ra, 
  Dec = gaia$dec, 
  keyvalues = im_frame$keyvalues
)
gaia$ra_fix = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time
gaia$dec_fix = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time
gaia_idx = gaia$xpix >= 0 & gaia$xpix <= imdim[1] & gaia$ypix >= 0 & gaia$ypix <= imdim[2] & gaia$classprob_dsc_combmod_star > 0.98 & gaia$phot_bp_rp_excess_factor < 2.5
gaia_idx[is.na(gaia_idx)] = FALSE
gaia_ra_dec = gaia[gaia_idx,]
gaia_trim = gaia_ra_dec[order(gaia_ra_dec$phot_g_mean_mag), ]
fwrite(
  gaia_trim,
  file = paste0(main_dir, "data/stars/", fname, "/", fname, "_trim_", gaia_cat_stub)
)

star_mask = profoundApplyMask(
  image = im_frame, 
  xcen = gaia_trim$xpix, 
  ycen = gaia_trim$ypix,
  xsize = 11,
  ysize = 11
)

profound_stars = profoundProFound(
  image = im_frame,
  segim = star_mask$mask,
  fluxtype = "Jansky",
  magzero = im_magzero,
  redosky = FALSE,
  redosegim = FALSE,
)

star_aperture_sizes = seq(0.1, 10.1, 0.1)
profound_stars_aper = profoundAperPhot(
  image = im_frame, 
  segim = star_mask$mask, 
  app_diam = star_aperture_sizes, 
  pixscale = pixscale(im_frame), 
  magzero = im_magzero,
  fluxtype = "Jansky", 
  correction = FALSE
)

star_radii = foreach(i = 1:dim(gaia_trim)[1], .combine = c) %do% {
  
  if(is.na(profound_stars$segstats$flux[i])){
    1L
  }else{
    # max_flux = profound_stars$segstats$flux[i]
    curve_of_growth = unname(unlist(c(profound_stars_aper[i, paste0("flux_app_", 1:length(star_aperture_sizes))])))
    max_flux = max(curve_of_growth)
    # rate_of_change = c(diff(curve_of_growth),0)
    
    magplot(
      star_aperture_sizes,
      curve_of_growth
    )
    smoothingSpline = smooth.spline(star_aperture_sizes, curve_of_growth, spar = 0.5)
    pred_lo = predict(
      smoothingSpline,
      (seq(0,10,0.01))
    )
    lines((seq(0,10,0.01)), pred_lo$y, col='red', lwd=2)
    
    rad = head(star_aperture_sizes[curve_of_growth/max_flux > 0.99],1) / pixscale(im_frame) * 2 ## make the thing larger
    
    rad_new = 2 * round(rad/2.0) + 1
  }
}

gc()
rm(star_mask)
rm(profound_stars)
rm(profound_stars_aper)
gc()

na_mask = profoundDilate(
  is.na(im_frame[,]$imDat), 
  15
)

star_mask_new = profoundApplyMask(
  image = im_frame, 
  xcen = gaia_trim$xpix, 
  ycen = gaia_trim$ypix,
  xsize = star_radii*6 + 1,
  ysize = star_radii*6 + 1
)

profoundSegimPlot(
  image = im_frame,
  segim = star_mask_new$mask
)
star_mask_new$stars = star_mask_new$mask
star_mask_new$mask = star_mask_new$mask | na_mask

saveRDS(
  star_mask_new,
  paste0(main_dir, "data/stars/", fname, "/", fname, "_GAIA_star_mask.rds")
)

im_list = Rfits_make_list(im_frame_stubs, extlist = 2)
filt_names = toupper(sub(".*_(f\\w+w)_.*", "\\1", im_frame_stubs))

magzero_list = unname(sapply(im_list, function(x){
  -2.5*log10(x$keyvalues$PIXAR_SR) - 6.10
}))

star_cat = foreach(i = 1:length(im_list), .combine = cbind) %do% {
  profound_aperture_phot = profoundAperPhot(
    image = im_list[[i]],
    segim = star_mask_new$stars,
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

star_cat$RA = gaia_trim$ra_fix
star_cat$DEC = gaia_trim$dec_fix

star_cat = star_cat[, c("RA", "DEC", setdiff(names(star_cat), c("RA", "DEC")))]
fwrite(
  star_cat,
  paste0(main_dir, "/data/stars/", fname, "/GAIA_", fname, "_MIRI_photometry.csv")
)
