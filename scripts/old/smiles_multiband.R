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
library(matrixStats)

main_dir = "~/Documents/Epochs/"
## Do aperture photometry from my segim mask
segim_stack = Rfits_read(
  paste0(main_dir, "data/images/smiles/smiles_F560W_F770W_invar_stack_propane_segim.fits")
)

segim_cat = fread(
  paste0(main_dir, "data/catalogues/smiles_miri_F560W_F770W_stack_profound_segstats.csv")
)

filt_names = read.table(
  file = paste0(main_dir, "data/throughputs/wavelengths.txt"),
)[25:32,2]

F560W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f560w_v1.0_drz.fits")
)
F770W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f770w_v1.0_drz.fits")
)
F1000W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f1000w_v1.0_drz.fits")
)
F1280W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f1280w_v1.0_drz.fits")
)
F1500W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f1500w_v1.0_drz.fits")
)
F1800W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f1800w_v1.0_drz.fits")
)
F2100W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f2100w_v1.0_drz.fits")
)
F2550W = Rfits_read(
  filename = paste0(main_dir, "data/images/smiles/hlsp_smiles_jwst_miri_goodss_f2550w_v1.0_drz.fits")
)

im_list = list(
  F560W, F770W, F1000W, F1280W, F1500W, F1800W, F2100W, F2550W
)
names(im_list) = filt_names

magzero_list = unname(sapply(im_list, function(x){
  -2.5*log10(x$SCI$keyvalues$PIXAR_SR) - 6.10
}))

# profound_cat = foreach(i = 1:length(im_list), .combine = cbind) %do% {
profound_cat = foreach(i = 1, .combine = cbind) %do% {

  message("Running ", filt_names[i])
  
  pro_measure_fstub = paste0(main_dir, "/data/profound/", filt_names[i], "_pro_measure.rds")
  if(file.exists(pro_measure_fstub)){
    pro_measure = readRDS(pro_measure_fstub)
  }else{
    pro_measure = profoundProFound(
      image = im_list[[i]]$SCI, 
      segim = segim_stack$SEGIM[,]$imDat,
      sky = 0,
      skyRMS = im_list[[i]]$ERR[,]$imDat,
      magzero = magzero_list[i],
      fluxtype = "Jansky",
      
      iters = 2,
      size = 5,
      
      box = 100,
      grid = 100,
      boxiters = 2,
      
      # measurement mode
      # static_photom = TRUE,
      redosegim = FALSE,
      redosky = FALSE,
      # stats
      
      boundstats = TRUE,
      segstats = TRUE,
      
      # misc
      pixelcov = FALSE,
      rem_mask = TRUE,
      lowmemory = FALSE
    )
    
    message("Saving profound...")
    saveRDS(
      pro_measure,
      paste0(main_dir, "/data/profound/", filt_names[i], "_pro_measure.rds")
    )
  }
  
  message("Aperture photometry...")
  r_ap = c(0.16, 0.32, 0.5, 0.6)
  profound_aperture_phot = profoundAperPhot(
    image = im_list[[i]]$SCI[,]$imDat,
    segim = segim_stack$SEGIM[,]$imDat,
    app_diam = r_ap/pixscale(im_list[[i]]$SCI) * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Jansky", 
    depth = 5,
    correction = TRUE
  )
  profound_aperture_phot_uncorr = profoundAperPhot(
    image = im_list[[i]]$SCI[,]$imDat,
    segim = segim_stack$SEGIM[,]$imDat,
    app_diam = r_ap/pixscale(im_list[[i]]$SCI) * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Jansky", 
    depth = 5,
    correction = FALSE
  )
  
  profound_aperture_err = profoundAperPhot(
    image = im_list[[i]]$ERR[,]$imDat^2,
    segim = segim_stack$SEGIM[,]$imDat,
    app_diam = r_ap/pixscale(im_list[[i]]$SCI) * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Raw", 
    depth = 5,
    correction = TRUE
  )
  profound_aperture_err[,grep("^flux_app", names(profound_aperture_err), value = TRUE)] = sqrt(profound_aperture_err[,grep("^flux_app", names(profound_aperture_err), value = TRUE)]) * 1e6 * pixarea(im_list[[i]]$SCI, unit = "str") ## put into Jansky
  
  random_aperture_radii = seq(0.05, 5.05, 0.05) * 0.5
  random_apertures = profoundAperRan(
    image = im_list[[i]]$SCI[,]$imDat,
    segim = pro_measure$segim,
    app_diam = random_aperture_radii * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Jansky", 
    Nran = 100000,
    depth = 4,
    correction = TRUE
  )
  apertures_sizes_in_pixels = (random_aperture_radii/pixscale(im_list[[i]]$SCI))^2 * pi
  random_apertures_model = lm(log10(random_apertures$errors[random_apertures$errors>0]) ~ log10(apertures_sizes_in_pixels[random_apertures$errors>0]))
  random_aperture_model_alpha <- 10^(coef(random_apertures_model)[1])
  random_aperture_model_beta  <- coef(random_apertures_model)[2]

  nearest_apertures_idx = RANN::nn2(
    data = random_apertures$AperPhot[, c("xcen", "ycen")],
    query = pro_measure$segstats[, c("xmax", "ymax")], 
    k = 200
  )
  ## calc local depths
  Nphot_app = profound_aperture_phot[,grep("^N_app", names(profound_aperture_err), value = TRUE)]
  local_depths_app = lapply(1:dim(nearest_apertures_idx$nn.idx)[1], function(x){
    closest_apertures = random_apertures$AperPhot[, paste0("flux_app_", sapply(Nphot_app[x,], function(y) which.min(abs(apertures_sizes_in_pixels - y))))]
    colSds(
      as.matrix(closest_apertures[nearest_apertures_idx$nn.idx[x,],]), na.rm = TRUE
    )
  })
  local_depths_app = do.call(rbind, local_depths_app)
  
  ## Scale the errors with the Number of pixels that lie in the aperture
  # profound_aperture_err_scaled_temp = profound_aperture_phot[,grep("^N_app", names(profound_aperture_err), value = TRUE)]^random_aperture_model_beta * random_aperture_model_alpha
  # names(profound_aperture_err_scaled_temp) = str_replace_all(names(profound_aperture_err_scaled_temp), "N", "flux")
  
  # m1 <- as.matrix(profound_aperture_err_scaled_temp)
  # m2 <- as.matrix(profound_aperture_err[, c("flux_app_1", "flux_app_2", "flux_app_3", "flux_app_4")])
  # profound_aperture_err_scaled <- as.data.frame(pmax(m1, m2))
  
  m1 <- as.matrix(local_depths_app)
  m2 <- as.matrix(profound_aperture_err[, c("flux_app_1", "flux_app_2", "flux_app_3", "flux_app_4")])
  profound_aperture_err_scaled <- as.data.frame(pmax(m1, m2))
  
  ## add 5pc in quadrature 
  profound_aperture_err_scaled_5pc = sqrt( (profound_aperture_err_scaled)^2 + (0.05*profound_aperture_phot[, c("flux_app_1", "flux_app_2", "flux_app_3", "flux_app_4")])^2 )
  
  local_depths = sapply(1:dim(nearest_apertures_idx$nn.idx)[1], function(x){
    closest_apertures = random_apertures$AperPhot[, paste0("flux_app_", which.min(abs(apertures_sizes_in_pixels - pro_measure$segstats$N100[x])))]
    sd(closest_apertures, na.rm = TRUE)
  })
  # pro_measure$segstats$flux_err_scaled = pmax(pro_measure$segstats$N100^random_aperture_model_beta * random_aperture_model_alpha, pro_measure$segstats$flux_err)
  pro_measure$segstats$flux_err_scaled = pmax(local_depths, pro_measure$segstats$flux_err)
  pro_measure$segstats$flux_err_scaled_5pc = sqrt( (pro_measure$segstats$flux_err_scaled)^2 + (0.05*pro_measure$segstats$flux)^2 )
  
  phot_cat = data.frame(
    ifelse(is.na(pro_measure$segstats$flux), 0.0, pro_measure$segstats$flux),
    ifelse(is.na(pro_measure$segstats$flux_err_scaled), 0.0, pro_measure$segstats$flux_err_scaled),
    ifelse(is.na(pro_measure$segstats$flux_err_scaled_5pc), 0.0, pro_measure$segstats$flux_err_scaled_5pc),
    ifelse(is.na(pro_measure$segstats$R50), 0.0, pro_measure$segstats$R50),
    ifelse(is.na(pro_measure$segstats$R90), 0.0, pro_measure$segstats$R90),
    ifelse(is.na(pro_measure$segstats$cenfrac), 0.0, pro_measure$segstats$cenfrac),
    ifelse(is.na(pro_measure$segstats$edge_frac), 0.0, pro_measure$segstats$edge_frac),
    ifelse(is.na(pro_measure$segstats$edge_excess), 0.0, pro_measure$segstats$edge_excess),
    # ifelse(is.na(pro_measure$segstats$flag_border), 0.0, pro_measure$segstats$flag_border),
    
    ifelse(is.na(profound_aperture_phot$flux_app_1), 0.0, profound_aperture_phot$flux_app_1),
    ifelse(is.na(profound_aperture_phot$flux_app_2), 0.0, profound_aperture_phot$flux_app_2),
    ifelse(is.na(profound_aperture_phot$flux_app_3), 0.0, profound_aperture_phot$flux_app_3),
    ifelse(is.na(profound_aperture_phot$flux_app_4), 0.0, profound_aperture_phot$flux_app_4),
    
    ifelse(is.na(profound_aperture_err_scaled[,1]), 0.0, profound_aperture_err_scaled[,1]),
    ifelse(is.na(profound_aperture_err_scaled[,2]), 0.0, profound_aperture_err_scaled[,2]),
    ifelse(is.na(profound_aperture_err_scaled[,3]), 0.0, profound_aperture_err_scaled[,3]),
    ifelse(is.na(profound_aperture_err_scaled[,4]), 0.0, profound_aperture_err_scaled[,4]),
    
    ifelse(is.na(profound_aperture_err_scaled_5pc[,1]), 0.0, profound_aperture_err_scaled_5pc[,1]),
    ifelse(is.na(profound_aperture_err_scaled_5pc[,2]), 0.0, profound_aperture_err_scaled_5pc[,2]),
    ifelse(is.na(profound_aperture_err_scaled_5pc[,3]), 0.0, profound_aperture_err_scaled_5pc[,3]),
    ifelse(is.na(profound_aperture_err_scaled_5pc[,4]), 0.0, profound_aperture_err_scaled_5pc[,4]),
    
    ifelse(is.na(profound_aperture_phot$N_app_1), 0.0, profound_aperture_phot$N_app_1),
    ifelse(is.na(profound_aperture_phot$N_app_2), 0.0, profound_aperture_phot$N_app_2),
    ifelse(is.na(profound_aperture_phot$N_app_3), 0.0, profound_aperture_phot$N_app_3),
    ifelse(is.na(profound_aperture_phot$N_app_4), 0.0, profound_aperture_phot$N_app_4),
    
    ifelse(is.na(profound_aperture_phot$flux_app_1/profound_aperture_phot_uncorr$flux_app_1), 1.0, profound_aperture_phot$flux_app_1/profound_aperture_phot_uncorr$flux_app_1),
    ifelse(is.na(profound_aperture_phot$flux_app_2/profound_aperture_phot_uncorr$flux_app_2), 1.0, profound_aperture_phot$flux_app_2/profound_aperture_phot_uncorr$flux_app_2),
    ifelse(is.na(profound_aperture_phot$flux_app_3/profound_aperture_phot_uncorr$flux_app_3), 1.0, profound_aperture_phot$flux_app_3/profound_aperture_phot_uncorr$flux_app_3),
    ifelse(is.na(profound_aperture_phot$flux_app_4/profound_aperture_phot_uncorr$flux_app_4), 1.0, profound_aperture_phot$flux_app_4/profound_aperture_phot_uncorr$flux_app_4),
    
    local_depths_app[,1] * 1,
    local_depths_app[,2] * 1,
    local_depths_app[,3] * 1,
    local_depths_app[,4] * 1
  )
  names(phot_cat) = c(
    ## main profound stuff
    paste0(filt_names[i], "_flux_tot"),
    paste0(filt_names[i], "_flux_tot_err"),
    paste0(filt_names[i], "_flux_tot_err_5pc"),
    paste0(filt_names[i], "_R50"),
    paste0(filt_names[i], "_R90"),
    paste0(filt_names[i], "_cenfrac"),
    paste0(filt_names[i], "_edge_frac"),
    paste0(filt_names[i], "_edge_excess"),
    # paste0(filt_names[i], "_flag_border"),
    
    paste0(filt_names[i], "_flux_aper16"),
    paste0(filt_names[i], "_flux_aper32"),
    paste0(filt_names[i], "_flux_aper50"),
    paste0(filt_names[i], "_flux_aper60"),
  
    paste0(filt_names[i], "_flux_aper16_err"),
    paste0(filt_names[i], "_flux_aper32_err"),
    paste0(filt_names[i], "_flux_aper50_err"),
    paste0(filt_names[i], "_flux_aper60_err"),
    
    paste0(filt_names[i], "_flux_aper16_err_5pc"),
    paste0(filt_names[i], "_flux_aper32_err_5pc"),
    paste0(filt_names[i], "_flux_aper50_err_5pc"),
    paste0(filt_names[i], "_flux_aper60_err_5pc"),
    
    paste0(filt_names[i], "_N_aper16"),
    paste0(filt_names[i], "_N_aper32"),
    paste0(filt_names[i], "_N_aper50"),
    paste0(filt_names[i], "_N_aper60"),
    
    paste0(filt_names[i], "_app_corr_aper16"),
    paste0(filt_names[i], "_app_corr_aper32"),
    paste0(filt_names[i], "_app_corr_aper50"),
    paste0(filt_names[i], "_app_corr_aper60"),
    
    paste0(filt_names[i], "_1sigma_loc_depth_aper16"),
    paste0(filt_names[i], "_1sigma_loc_depth_aper32"),
    paste0(filt_names[i], "_1sigma_loc_depth_aper50"),
    paste0(filt_names[i], "_1sigma_loc_depth_aper60")
  )
  gc()
  rm(pro_measure)
  gc()
  message("Done \n")
  return(phot_cat)
}

## noise in detect band
random_aperture_radii = seq(0.05, 5.05, 0.05)
random_apertures = profoundAperRan(
  image = segim_stack$image[,]$imDat,
  segim = segim_stack$SEGIM[,]$imDat,
  app_diam = random_aperture_radii * 2, 
  magzero = magzero_list[i],
  pixscale = pixscale(segim_stack$image),
  fluxtype = "Jansky", 
  Nran = 1000,
  depth = 5,
  correction = TRUE
)
apertures_sizes_in_pixels = (random_aperture_radii/pixscale(segim_stack$image))^2 * pi
random_apertures$errors[random_apertures$errors == 0] = max(random_apertures$errors)
random_apertures_model = lm(log10(random_apertures$errors) ~ log10(apertures_sizes_in_pixels))
random_aperture_model_alpha <- 10^(coef(random_apertures_model)[1])
random_aperture_model_beta  <- coef(random_apertures_model)[2]

profound_cat$RA = segim_cat$RAmax
profound_cat$DEC = segim_cat$Decmax
profound_cat$segID = segim_cat$segID
profound_cat$flux = segim_cat$flux
profound_cat$flux_err = pmax(segim_cat$N100^random_aperture_model_beta * random_aperture_model_alpha, segim_cat$flux_err)
profound_cat$mag = segim_cat$mag
profound_cat$mag_err = (2.5 / log(10)) * abs(profound_cat$flux_err / profound_cat$flux)
profound_cat$R50 = segim_cat$R50
profound_cat$cenfrac = segim_cat$cenfrac
profound_cat$edge_frac = segim_cat$edge_frac
profound_cat$edge_excess = segim_cat$edge_excess
profound_cat$Nedge = segim_cat$Nedge
profound_cat$Nsky = segim_cat$Nsky
profound_cat$Nobject = segim_cat$Nobject
profound_cat$Nborder = segim_cat$Nborder
profound_cat$Nmask = segim_cat$Nmask
profound_cat$N100 = segim_cat$N100
# profound_cat$flag_border = segim_cat$flag_border
profound_cat = profound_cat[, c("RA", "DEC", setdiff(names(profound_cat), c("RA", "DEC")))]
fwrite(
  profound_cat,
  paste0(main_dir, "/data/catalogues/smiles_profound_MIRI_photometry_cat.csv")
)

#