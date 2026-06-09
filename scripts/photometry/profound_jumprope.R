setwd("~/Documents/JUMPROPE/PROFOUND/")

source("all_codes_profound.R")

main_dir = "~/mnt/morgan/data/data/JP/"

convert_to_propane(
  input_args = list(
    in_file = "~/Documents/Epochs/data/images/smiles/hlsp_smiles_jwst_miri_goodss_f770w_v1.0_drz.fits",
    ref_dir = "~/Documents/Epochs/JP/",
    VID = "smiles",
    MODULE = "smiles",
    PIXSCALE = "miri"
  )
)

foo = Rfits_read(
  filename = paste0("~/Documents/Epochs/data/images/smiles/smiles_F560W_F770W_invar_stack_propane_segim.fits")
)
# saveRDS(
#   list("mask" = foo$STARMASK[,]$imDat),
#   "~/Documents/Epochs/JP/ProFound/Star_Masks/smiles/smiles/smiles_smiles_star_mask.rds"
# )
# 
mask = foo$MASK[,]$imDat

profound_detect_master_smiles = function(frame, skyRMS, star_mask, pix_mask = mask, segim = NULL){
  ## Inject this function into the detect runs
  ## frame is the stack with WCS
  ## stack is for the skyRMS matrix
  ## star_mask should be 1 or 0 matrix
  
  if(is.null(pix_mask)){
    mask = is.na(frame$imDat) | star_mask
  }else{
    mask = is.na(frame$imDat) | star_mask | pix_mask
  }
  pro = profoundProFound(
    image = frame,
    mask = mask | (frame$imDat == 0) | (is.na(frame$imDat) | (is.infinite(frame$imDat))),
    segim = segim,
    skyRMS = skyRMS,
    magzero = 23.9,
    
    sky = 0,
    skycut = 1.5,
    pixcut = 7.0,
    
    tolerance = 4.0,
    reltol = 0,
    cliptol = 30,
    
    iters = 4,
    
    sigma = 1.0,
    box = 100,
    grid = 100,
    boxiters = 0,
    
    roughpedestal = FALSE,
    pixelcov = FALSE,
    boundstats = TRUE,
    groupstats = TRUE,
    nearstats = TRUE, 
    redosky = TRUE,
    rem_mask = TRUE,
    # redoskysize = 11,
    
    fluxtype = "Jansky",
    pixscale = pixscale(frame, unit = "asec"),
    
    verbose = FALSE,
  )
  message("Finished profound")
  return(pro)
}

do_detect(
  input_args = list(
    ref_dir = "~/Documents/Epochs/JP/",
    VID = "smiles",
    MODULE = "smiles",
    PIXSCALE = "miri"
  ),
  detect_bands = "f560w|f770w",
  profound_function = profound_detect_master_smiles
)

profound_measure_master_smiles = function(super_img, inVar, segim, mask, redosegim = TRUE, sky = 0, redosky = TRUE){
  ## Do multiband measurements with ProFound
  ## set sky = 0 and redosky = FALSE if frames are already sky subtracted and you don't want another sky model
  
  if(inherits(inVar, 'Rfits_image')){
    skyRMS = inVar$imDat^-0.5
  }else if(is.null(inVar)){
    skyRMS = NULL
  }else{
    skyRMS = inVar^-0.5
  }
  
  if(redosegim){
    iters = 2
  }else{
    iters = 0
  }
  
  super_pro = profoundProFound(
    super_img,
    magzero = 23.9,
    mask = mask | super_img$imDat == 0 | is.na(super_img$imDat) | is.infinite(super_img$imDat),
    # detection, segmentation and dilation
    segim = segim,
    redosegim = redosegim,
    
    size = 5,
    iters = iters,
    
    # sky estimate
    sky = sky,        # Do we model the sky again?
    skyRMS = skyRMS,
    box = 100,
    grid = 100,
    boxiters = 2,
    # measurement mode
    dotot = TRUE,
    docol = TRUE,
    dogrp = TRUE,
    # stats
    boundstats = TRUE,
    segstats = TRUE,
    # misc
    pixelcov = FALSE,
    rem_mask = TRUE,
    redosky = redosky,

    # image units
    fluxtype = "Jansky",
    pixscale = pixscale(super_img, unit = "asec"),
    
    verbose = FALSE,
  )
  return(super_pro)
}
do_measure(
  input_args = list(
    ref_dir = main_dir,
    VID = "smiles",
    MODULE = "smiles",
    PIXSCALE = "smiles"
  ), 
  profound_function = profound_measure_master_smiles
)
