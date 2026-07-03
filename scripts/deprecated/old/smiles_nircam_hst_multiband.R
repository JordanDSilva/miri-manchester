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

args = commandArgs(trailingOnly = TRUE)
# field_name = args[1]
field_name = "east"

main_dir = "~/Documents/Epochs/"
## Do aperture photometry from my segim mask
segim_stack = Rfits_read(
  paste0(main_dir, "data/images/smiles_F560W_F770W_invar_stack_propane_segim.fits")
)

segim_cat = fread(
  paste0(main_dir, "data/catalogues/smiles_miri_F560W_F770W_stack_profound_segstats.csv")
)

filt_names = read.table(
  file = paste0(main_dir, "data/throughputs/wavelengths.txt"),
)[1:23, 2]

# field_name = "east"


field_files = list.files(paste0(main_dir, "/data/images/combine/", field_name), full.names = TRUE)
im_list = lapply(
  field_files,
  function(x){
    Rfits_read(x)
  }
)

filt_names = toupper(regmatches(field_files, regexpr("F\\d+[WML]P?", field_files, ignore.case = TRUE)))

names(im_list) = filt_names

magzero_list = unname(sapply(im_list, function(x){
  if(is.null(x$SCI$keyvalues$ZEROPNT)){
    -2.5*log10(x$SCI$keyvalues$PIXAR_SR) - 6.10 ## using NIRCam pixel area
  }else{
    x$SCI$keyvalues$ZEROPNT
  }
}))

profound_cat = foreach(i = 1:length(im_list), .combine = cbind) %do% {
  
  message("Running ", filt_names[i])
  
  pro_measure_fstub = paste0(main_dir, "/data/profound/", filt_names[i], "_pro_measure.rds")
  if(file.exists(pro_measure_fstub)){
    pro_measure = readRDS(pro_measure_fstub)
  }else{
    pro_measure = profoundProFound(
      image = im_list[[i]]$SCI, 
      segim = segim_stack$SEGIM[,]$imDat,
      sky = 0,
      mask = im_list[[i]]$SCI[,]$imDat == 0 | is.infinite(im_list[[i]]$SCI[,]$imDat),
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
      paste0(main_dir, "/data/profound/", field_name, "/", filt_names[i], "_pro_measure.rds")
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
  err_corr_factor = 10^((magzero_list[i] - 8.9) * -0.4)
  profound_aperture_err[,grep("^flux_app", names(profound_aperture_err), value = TRUE)] = sqrt(profound_aperture_err[,grep("^flux_app", names(profound_aperture_err), value = TRUE)]) * err_corr_factor ## put into Jansky
  
  image_imdat = im_list[[i]]$SCI[,]$imDat
  image_imdat[image_imdat == 0] = NA
  random_aperture_radii = seq(0.05, 5.05, 0.05) * 0.5
  random_apertures = profoundAperRan(
    image = image_imdat,
    segim = pro_measure$segim,
    app_diam = random_aperture_radii * 2, 
    magzero = magzero_list[i],
    pixscale = pixscale(im_list[[i]]$SCI),
    fluxtype = "Jansky", 
    Nran = 1000,
    depth = 5,
    correction = TRUE
  )
  apertures_sizes_in_pixels = (random_aperture_radii/pixscale(im_list[[i]]$SCI))^2 * pi
  random_apertures_model = lm(log10(random_apertures$errors) ~ log10(apertures_sizes_in_pixels))
  random_aperture_model_alpha <- 10^(coef(random_apertures_model)[1])
  random_aperture_model_beta  <- coef(random_apertures_model)[2]
  
  nearest_apertures_idx = RANN::nn2(
    data = random_apertures$AperPhot[, c("xcen", "ycen")],
    query = pro_measure$segstats[, c("xmax", "ymax")], 
    k = 200
  )
  ## calc local depths
  Nphot_app = profound_aperture_phot[,grep("^N_app", names(profound_aperture_err), value = TRUE)]
  local_depths = lapply(1:dim(nearest_apertures_idx$nn.idx)[1], function(x){
    closest_apertures = random_apertures$AperPhot[, paste0("flux_app_", sapply(Nphot_app[x,], function(x) which.min(abs(apertures_sizes_in_pixels - x))))]
    colSds(
      as.matrix(closest_apertures[nearest_apertures_idx$nn.idx[x,],])
    )
  })
  local_depths = do.call(rbind, local_depths)
  
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
    
    local_depths_app[,1] * 5,
    local_depths_app[,2] * 5,
    local_depths_app[,3] * 5,
    local_depths_app[,4] * 5
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
    
    paste0(filt_names[i], "_5sigma_loc_depth_aper16"),
    paste0(filt_names[i], "_5sigma_loc_depth_aper32"),
    paste0(filt_names[i], "_5sigma_loc_depth_aper50"),
    paste0(filt_names[i], "_5sigma_loc_depth_aper60")
  )
  gc()
  rm(pro_measure)
  gc()
  message("Done \n")
  return(phot_cat)
}

fwrite(
  profound_cat,
  paste0(main_dir, "/data/profound/", filt_names[i], "_", field_name, "_nircam_hst.csv")
)

magplot(
  -2.5*log10(profound_cat$F090W_flux_tot)+8.9, 
  -2.5*log10(profound_cat$F814W_flux_tot)+8.9,
  xlim = c(16,40),
  ylim = c(16,40)
)
abline(0,1)


combine_cat = fread("~/Documents/Epochs/data/catalogues/smiles_profound_MIRI_matched_EPOCHSv2_photometry_cat.csv")
magplot(
  -2.5*log10(combine_cat$JWST_NIRCAM_F090W_flux_aper32)+8.9,
  -2.5*log10(profound_cat$F090W_flux_aper32)+8.9
)
abline(0,1)

maghist(
  combine_cat$JWST_NIRCAM_F090W_flux_aper32/profound_cat$F090W_flux_tot
)

maghist(
  phot_cat$F090W_flux_tot/phot_cat$F090W_flux_tot_err
)
maghist(
  phot_cat$F090W_flux_aper16/phot_cat$F090W_flux_aper16_err
)
