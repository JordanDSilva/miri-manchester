library(Rwcs)
library(Rfits)
library(ProPane)
library(ProFound)
library(data.table)
library(foreach)
library(magicaxis)
library(scales)

main_dir = "~/Documents/Epochs/"

## load in images
F560W = Rfits_read(
  filename = paste0(main_dir, "data/images/hlsp_smiles_jwst_miri_goodss_f560w_v1.0_drz.fits")
)
F770W = Rfits_read(
  filename = paste0(main_dir, "/data/images/hlsp_smiles_jwst_miri_goodss_f770w_v1.0_drz.fits")
)

F560W_F770W_stack = propaneStackWarpInVar(
  image_list = list(F560W$SCI, F770W$SCI), 
  inVar_list = list(F560W$ERR[,]$imDat^-2, F770W$ERR[,]$imDat^-2), 
  magzero_in = c(
    -2.5*log10(1e6 * pixarea(F560W$SCI, unit = "str")) + 8.9, 
    -2.5*log10(1e6 * pixarea(F770W$SCI, unit = "str")) + 8.9
  ), 
  magzero_out = -2.5*log10(1e6 * pixarea(F560W$SCI, unit = "str")) + 8.9,
  keyvalues_out = F560W$SCI$keyvalues,
  cores = 1,
  cores_warp = 1,
  direction = "forward"
)

## Load in GAIA catalogue 
imdim = dim(F560W_F770W_stack$image)
im_ra_dec = Rwcs_p2s(x = imdim[1], y = imdim[2], keyvalues = F560W_F770W_stack$image$keyvalues) #get extent of frame in RA and DEC coords
gaia = data.frame(
  fread(
    file = paste0(main_dir, "data/stars/GAIA_53.14977_-27.80586.csv")
  )
)
delta_time = 2023 - gaia$ref_epoch ## observed in December 2022?
gaia[, c("xpix", "ypix")] = Rwcs_s2p(RA = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time, ## Do proper motion correction
                                     Dec = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time, ## pmra[pmdec] = mas/yr
                                     keyvalues = F560W$SCI$keyvalues)
gaia[, c("xpix_nopm", "ypix_nopm")] = Rwcs_s2p(RA = gaia$ra, 
                                               Dec = gaia$dec, 
                                               keyvalues = F560W$SCI$keyvalues)
gaia$ra_fix = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time
gaia$dec_fix = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time
gaia_idx = gaia$xpix >= 0 & gaia$xpix <= imdim[1] & gaia$ypix >= 0 & gaia$ypix <= imdim[2] & gaia$classprob_dsc_combmod_star > 0.98 & gaia$phot_bp_rp_excess_factor < 2.5
gaia_idx[is.na(gaia_idx)] = FALSE
gaia_ra_dec = gaia[gaia_idx,]
gaia_trim = gaia_ra_dec[order(gaia_ra_dec$phot_g_mean_mag), ]
fwrite(
  gaia_trim,
  file = paste0(main_dir, "data/stars/GAIA_53.14977_-27.80586_trim.csv")
)

star_mask = profoundApplyMask(
  image = F560W_F770W_stack$image, 
  xcen = gaia_trim$xpix, 
  ycen = gaia_trim$ypix,
  xsize = 11,
  ysize = 11
)

profound_stars = profoundProFound(
  image = F560W_F770W_stack$image,
  segim = star_mask$mask,
  fluxtype = "Jansky",
  magzero = F560W_F770W_stack$image$keyvalues$MAGZERO,
  redosky = FALSE,
  redosegim = FALSE,
)

star_aperture_sizes = seq(0.1, 10.1, 0.1)
profound_stars_aper = profoundAperPhot(
  image = F560W_F770W_stack$image, 
  segim = star_mask$mask, 
  app_diam = star_aperture_sizes, 
  pixscale = pixscale(F560W_F770W_stack$image), 
  magzero = F560W_F770W_stack$image$keyvalues$MAGZERO,
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

    rad = head(star_aperture_sizes[curve_of_growth/max_flux > 0.99],1) / pixscale(F560W$SCI) * 2 ## make the thing larger
    
    rad_new = 2 * round(rad/2.0) + 1
  }
}

gc()
rm(star_mask)
rm(profound_stars)
rm(profound_stars_aper)
gc()

na_mask = profoundDilate(
  is.na(F560W_F770W_stack$image[,]$imDat), 
  15
)

star_mask_new = profoundApplyMask(
  image = F560W_F770W_stack$image, 
  xcen = gaia_trim$xpix, 
  ycen = gaia_trim$ypix,
  xsize = star_radii*6 + 1,
  ysize = star_radii*6 + 1
)

pro_detect = profoundProFound(
  image = F560W_F770W_stack$image,
  sky = 0,
  skyRMS = F560W_F770W_stack$inVar$imDat^-0.5,
  mask = star_mask_new$mask > 0 | na_mask,
  magzero = F560W_F770W_stack$image$keyvalues$MAGZERO, 
  pixscale = pixscale(F560W_F770W_stack$image),
    
  pixcut = 7.0,
  
  tolerance = 4.0,
  reltol = 0,
  cliptol = 30,
  
  iters = 4,
  
  sigma = 1.0,
  box = 100,
  grid = 100,
  boxiters = 0,
  
  fluxtype = "Jansky",
  
  boundstats = TRUE,
  groupstats = TRUE,
  
  redosky = TRUE,
  rem_mask = TRUE
)

segim_orig_fits = Rfits_create_image(
  data = pro_detect$segim,
  keyvalues = F560W_F770W_stack$image$keyvalues
)
segim_orig_fits$keyvalues$EXTNAME = "SEGIM_ORIG"
F560W_F770W_stack$SEGIM_ORIG = segim_orig_fits

segim_fits = Rfits_create_image(
  data = pro_detect$segim,
  keyvalues = F560W_F770W_stack$image$keyvalues
)
segim_fits$keyvalues$EXTNAME = "SEGIM"
F560W_F770W_stack$SEGIM = segim_fits

mask_fits = Rfits_create_image(
  data = pro_detect$mask,
  keyvalues = F560W_F770W_stack$image$keyvalues
)
mask_fits$keyvalues$EXTNAME = "MASK"
F560W_F770W_stack$MASK = mask_fits

star_mask_fits = Rfits_create_image(
  data = star_mask_new$mask,
  keyvalues = F560W_F770W_stack$image$keyvalues
)
star_mask_fits$keyvalues$EXTNAME = "STARMASK"
F560W_F770W_stack$STARMASK = star_mask_fits

profoundSegimPlot(
  image = F560W_F770W_stack$image[53.14416, -27.81356, box = 200, type = "coord"]$imDat,
  segim = F560W_F770W_stack$SEGIM[53.14416, -27.81356, box = 200, type = "coord"]$imDat,
  locut = 0.5,
  hicut = 0.9999
)

Rfits_write(
  data = F560W_F770W_stack,
  filename = paste0(main_dir, "data/images/smiles_F560W_F770W_invar_stack_propane_segim.fits")
)

fwrite(
   pro_detect$segstats, 
   paste0(main_dir, "data/catalogues/smiles_miri_F560W_F770W_stack_profound_segstats.csv")
)

saveRDS(
  pro_detect,
  paste0(main_dir, "data/profound/smiles_miri_F560W_F770W_stack_pro_detect.rds")
)

pdf(paste0(main_dir, "data/plots/smiles_miri_F560W_F770W_profound.pdf"))
plot(pro_detect)
dev.off()

pdf(paste0(main_dir, "data/plots/smiles_miri_F560W_F770W_segim.pdf"))
profoundSegimPlot(pro_detect, sparse = 1, flip = TRUE, locut = 0.5, hicut = 0.995)
dev.off()

save.image(
  file = paste0(main_dir, "data/profound/source_detection.rda")
)
