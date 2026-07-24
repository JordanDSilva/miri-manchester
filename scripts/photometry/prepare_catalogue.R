library(ProFound)
library(Rfits)
library(Rwcs)
library(celestial)
library(ProSpect)
library(data.table)
library(dplyr)
library(foreach)
library(stringr)
library(magicaxis)
library(scales)
library(matrixStats)
library(rhdf5)

do_stuff_goodss = function(fname){
  ## change as necessary
  star_cat = data.frame(fread(
    paste0("~/Documents/Epochs/data/stars/", fname, "/GAIA_", fname, "_MIRI_photometry.csv")
  ))
  jades_F277W_warp = Rfits_point("~/Documents/Epochs/data/images/smiles/warped/JADES-DR3-GS_F277W_smiles_warp.fits")
  ## needed to discern if point lies out of footprint
  ## load in Epochsv2
  ## change as necessary - GOODS-S default
  GS_epochsV2 = data.frame(fread(file = "~/Documents/Epochs/data/catalogues/EPOCHS-v2_GS_Sel-F277W+F356W+F444W_v14_psfmatch_F444W_empirical+v14_psfmatch_F444W_empirical.csv"))

  profound_cat = data.frame(fread(
    paste0("~/Documents/Epochs/data/catalogues/jp/", fname, "/", fname, "_", fname, "_", fname, "_segstats.csv")
  ))
  phot_cat = data.frame(fread(
    paste0("~/Documents/Epochs/data/catalogues/jp/", fname, "/", fname, "_", fname, "_", fname, "_photometry.csv")
  ))
  aper_cat = data.frame(fread(
    paste0("~/Documents/Epochs/data/catalogues/jp/", fname, "/", fname, "_", fname, "_", fname, "_aperture_photometry.csv")
  ))

  MW_foreground_extinction = data.frame(fread(
    paste0("~/Documents/Epochs/data/throughputs/MW_foreground_extinction.csv")
  ))
  stack = Rfits_read_image(
    paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"),
    ext = 1, header = TRUE
  )
  segim = Rfits_read_image(
    paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"),
    ext = "segim", header = FALSE
  )
  segim_wcs = Rfits_read_image(
    paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"),
    ext = "segim", header = TRUE
  )
  frame_info = data.frame(fread(paste0("~/Documents/Epochs/data/depths/", fname, "_profound_frame_stats.csv")))
  # 
  # # star flags
  # magplot(
  #   -2.5*log10(aper_cat$F560W_flux_app_0.32)+8.9,
  #   -2.5*log10(aper_cat$F1280W_flux_app_0.32/aper_cat$F560W_flux_app_0.32),
  #   pch = "."
  # )
  # points(
  #   -2.5*log10(star_cat$F560W_flux_aper32)+8.9,
  #   -2.5*log10(star_cat$F1280W_flux_aper32/star_cat$F560W_flux_aper32),
  #   pch = 16, col = "blue"
  # )
  # points(
  #   -2.5*log10(aper_cat$F560W_flux_app_0.32[174])+8.9,
  #   -2.5*log10(aper_cat$F1280W_flux_app_0.32[174]/aper_cat$F560W_flux_app_0.32[174]),
  #   pch = 16, col = "red"
  # )
  # star_colour_lo = quantile(
  #   -2.5*log10(star_cat$F1280W_flux_aper32/star_cat$F560W_flux_aper32),
  #   probs = 0.05,
  #   na.rm = TRUE
  # )
  # star_colour_hi = quantile(
  #   -2.5*log10(star_cat$F1280W_flux_aper32/star_cat$F560W_flux_aper32),
  #   probs = 0.95,
  #   na.rm = TRUE
  # )
  # abline(
  #   h = star_colour_lo
  # )
  # abline(
  #   h = star_colour_hi
  # )

  f560w_f770w_color = -2.5*log10(aper_cat$F560W_flux_app_0.16/aper_cat$F770W_flux_app_0.16)
  col_cut_median = quantile(f560w_f770w_color, 0.5, na.rm = TRUE) 
  flag_hot_pixel = profound_cat$R50 <= 0.1 & profound_cat$flux/profound_cat$det_1sig_local_depth >= 10 & (is.na(f560w_f770w_color) | f560w_f770w_color > col_cut_median+0.5 | f560w_f770w_color < col_cut_median-0.5)
  
  # for(i in which(flag_hot_pixel)){
  #   pdf(
  #     paste0("~/Documents/Epochs/data/plots/artefacts/hotpixels/", i, ".pdf"), width = 5, height = 5
  #   )
  #   plot(
  #     stack[
  #       profound_cat$RA[i], 
  #       profound_cat$DEC[i], 
  #       box = 100, 
  #       type = "coord"
  #     ], 
  #     flip = TRUE
  #   )
  #   dev.off()
  # }
  # 
  # which(profound_cat$mag <= 24 & profound_cat$R50 <= 0.1)
  # plot(
  #   jadesEast_F277W[
  #     profound_cat$RA[520], 
  #     profound_cat$DEC[520], 
  #     box = 100, 
  #     type = "coord", 
  #   ], flip = TRUE
  # )
  
  # flag_star_idx = (-2.5*log10(aper_cat$F1280W_flux_app_0.32/aper_cat$F560W_flux_app_0.32) > star_colour_lo & 
  #   -2.5*log10(aper_cat$F1280W_flux_app_0.32/aper_cat$F560W_flux_app_0.32) < star_colour_hi) & aper_cat$F560W_flux_app_0.6/aper_cat$F560W_flux_app_0.16 >= 1.5 & -2.5*log10(aper_cat$F560W_flux_app_0.32)+8.9 <= 23
  #   # profound_cat$R50 <= 0.1
  # flag_star_idx[is.na(flag_star_idx)] = FALSE
  
  ## mask flags
  flag_mask_idx = profound_cat$Nmask/profound_cat$N100 > 0
  
  ## detection significance flags
  flag_detection_idx = (profound_cat$flux/profound_cat$det_1sig_local_depth < 4)
  
  ## flag no nircam
  match_idx = coordmatch(
    coordref = GS_epochsV2[, c("ALPHA_J2000", "DELTA_J2000")], 
    coordcompare = profound_cat[, c("RAmax", "Decmax")],
    rad = 0.3
  )
  flag_epochs_unmatch_idx = !(1:dim(profound_cat)[1] %in% match_idx$bestmatch$compareID)
  
  flag_outside_nircam_idx = rep(0L, length(flag_epochs_unmatch_idx))
  for(i in which(flag_epochs_unmatch_idx)){
    cutout = jades_F277W_warp[
      profound_cat$RAmax[i],
      profound_cat$Decmax[i],
      box = 3, ## about 16asec ?
      type = "coord"
    ]$imDat
    
    if( sum(cutout) == 0 | all(is.na(cutout))){
      flag_outside_nircam_idx[i] = 1L
    }
  }
  
  ## combine it all
  flag_all_idx = flag_hot_pixel | flag_mask_idx | flag_epochs_unmatch_idx | flag_outside_nircam_idx
  
  profound_cat$flag_hotpixel = flag_hot_pixel
  # profound_cat$flag_star = flag_star_idx
  profound_cat$flag_mask = flag_mask_idx
  profound_cat$flag_det = flag_detection_idx
  profound_cat$flag_epochs_unmatch = flag_epochs_unmatch_idx
  profound_cat$flag_outside_nircam = flag_outside_nircam_idx
  profound_cat$flag = flag_all_idx
  
  message("MATCHING TO EPOCHSv2")
  if(fname == 'smiles'){
    ## Match NIRCAM photometry
    epochsv2_fluxes = data.frame(GS_epochsV2)[,grep(pattern = glob2rx("FLUX_APER*aper_corr_Jy"), x = names(GS_epochsV2), value = TRUE)]
    matched_photometry = data.frame(matrix(0.0, nrow = dim(profound_cat)[1], ncol = dim(epochsv2_fluxes)[2]))
    matched_photometry[match_idx$bestmatch$compareID,] = epochsv2_fluxes[match_idx$bestmatch$refID, ]
    names(matched_photometry) = names(epochsv2_fluxes)
    
    epochsv2_flux_errs = data.frame(GS_epochsV2)[ ,grep(pattern = glob2rx("FLUXERR_APER*_loc_depth_5pc_Jy"), x = names(GS_epochsV2), value = TRUE)]
    ## load in depths from GALFIND
    JADES_DR3_GS_EAST_depths = read.table(file = paste0("~/Documents/Epochs/data/depths/JADES-DR3-GS-East_depths.ecsv"))
    JADES_DR3_GS_EAST_depths_medians = as.numeric(head(JADES_DR3_GS_EAST_depths$V6[JADES_DR3_GS_EAST_depths$V5 == "all"],-1))
    JADES_DR3_GS_SOUTH_depths = read.table(file = paste0("~/Documents/Epochs/data/depths/JADES-DR3-GS-South_depths.ecsv"))
    JADES_DR3_GS_SOUTH_depths_medians = as.numeric(head(JADES_DR3_GS_SOUTH_depths$V6[JADES_DR3_GS_SOUTH_depths$V5 == "all"],-1))
    
    # epochsv2_flux_errs[epochsv2_flux_errs == -999] = 0.05*epochsv2_fluxes[epochsv2_flux_errs == -999]
    matched_photometry_errs = data.frame(matrix(0.0, nrow = dim(profound_cat)[1], ncol = dim(epochsv2_flux_errs)[2]))
    matched_photometry_errs[match_idx$bestmatch$compareID,] = epochsv2_flux_errs[match_idx$bestmatch$refID, ]
    names(matched_photometry_errs) = names(epochsv2_flux_errs)
    
    matched_survey_names = rep(0L, dim(profound_cat)[1])
    matched_survey_names[match_idx$bestmatch$compareID] = GS_epochsV2$SURVEY[match_idx$bestmatch$refID]
    
    ## Replace -999 with median of the local depths calculated with the GALFIND code
    ## Use maximum of 5% or the local depth median?
    for(i in 1:length(JADES_DR3_GS_EAST_depths_medians)){
      flux_col_name = paste0("FLUX_APER_", head(JADES_DR3_GS_EAST_depths$V1[JADES_DR3_GS_EAST_depths$V5 == "all"],-1)[i], "_aper_corr_Jy")
      flux_err_col_name = paste0("FLUXERR_APER_", head(JADES_DR3_GS_EAST_depths$V1[JADES_DR3_GS_EAST_depths$V5 == "all"],-1)[i], "_loc_depth_5pc_Jy")
      idx = matched_survey_names == "JADES-DR3-GS-East" & matched_photometry_errs[[flux_err_col_name]] == -999
      if(sum(idx) > 0){
        loc_depth_err = 10^(-0.4*(JADES_DR3_GS_EAST_depths_medians[i]-8.9))/5.0
        matched_photometry_errs[[flux_err_col_name]][idx] = pmax(0.05*matched_photometry[[flux_col_name]][idx],loc_depth_err)
      }
    }
    for(i in 1:length(JADES_DR3_GS_SOUTH_depths_medians)){
      flux_col_name = paste0("FLUX_APER_", head(JADES_DR3_GS_SOUTH_depths$V1[JADES_DR3_GS_SOUTH_depths$V5 == "all"],-1)[i], "_aper_corr_Jy")
      flux_err_col_name = paste0("FLUXERR_APER_", head(JADES_DR3_GS_SOUTH_depths$V1[JADES_DR3_GS_SOUTH_depths$V5 == "all"],-1)[i], "_loc_depth_5pc_Jy")
      idx = matched_survey_names == "JADES-DR3-GS-South" & matched_photometry_errs[[flux_err_col_name]] == -999
      if(sum(idx) > 0){
        loc_depth_err = 10^(-0.4*(JADES_DR3_GS_SOUTH_depths_medians[i]-8.9))/5.0
        matched_photometry_errs[[flux_err_col_name]][idx] = pmax(0.05*matched_photometry[[flux_col_name]][idx],loc_depth_err)
      }
    }
    
    epochsv2_redshifts = data.frame(GS_epochsV2)[ ,grep(pattern = glob2rx("zbest_*"), x = names(GS_epochsV2), value = TRUE)]
    matched_redshifts = data.frame(matrix(-999, nrow = dim(profound_cat)[1], ncol = dim(epochsv2_redshifts)[2]))
    matched_redshifts[match_idx$bestmatch$compareID,] = epochsv2_redshifts[match_idx$bestmatch$refID, ]
    
    names(matched_photometry) = c(
      "F435W_flux_aper32",
      "F606W_flux_aper32",
      "F070W_flux_aper32",
      "F775W_flux_aper32",
      "F814W_flux_aper32",
      "F850LP_flux_aper32",
      "F090W_flux_aper32",
      "F115W_flux_aper32",
      "F150W_flux_aper32",
      "F200W_flux_aper32",
      "F277W_flux_aper32",
      "F335M_flux_aper32",
      "F356W_flux_aper32",
      "F410M_flux_aper32",
      "F444W_flux_aper32",
      "F162M_flux_aper32",
      "F182M_flux_aper32",
      "F210M_flux_aper32",
      "F250M_flux_aper32",
      "F300M_flux_aper32",
      "F430M_flux_aper32",
      "F460M_flux_aper32",
      "F480M_flux_aper32"
    )
    names(matched_photometry_errs) = c(
      "F435W_flux_aper32_err_5pc",
      "F606W_flux_aper32_err_5pc",
      "F070W_flux_aper32_err_5pc",
      "F775W_flux_aper32_err_5pc",
      "F814W_flux_aper32_err_5pc",
      "F850LP_flux_aper32_err_5pc",
      "F090W_flux_aper32_err_5pc",
      "F115W_flux_aper32_err_5pc",
      "F150W_flux_aper32_err_5pc",
      "F200W_flux_aper32_err_5pc",
      "F277W_flux_aper32_err_5pc",
      "F335M_flux_aper32_err_5pc",
      "F356W_flux_aper32_err_5pc",
      "F410M_flux_aper32_err_5pc",
      "F444W_flux_aper32_err_5pc",
      "F162M_flux_aper32_err_5pc",
      "F182M_flux_aper32_err_5pc",
      "F210M_flux_aper32_err_5pc",
      "F250M_flux_aper32_err_5pc",
      "F300M_flux_aper32_err_5pc",
      "F430M_flux_aper32_err_5pc",
      "F460M_flux_aper32_err_5pc",
      "F480M_flux_aper32_err_5pc"
    )
    names(matched_redshifts) = c(
      "zbest_sfhz_zfree",
      "zbest_16_sfhz_zfree",
      "zbest_84_sfhz_zfree"
    )
  }
  
  ## Combine forced aperture photometry from my HST/NIRCam warps
  combined_cat = data.frame(
    cbind(
      profound_cat,
      phot_cat,
      aper_cat
    )
  )
  # names(combined_cat) = str_remove_all(names(combined_cat), "JWST_MIRI_")
  
  AperStats = data.frame(
    aper_cat[, grep(pattern = glob2rx("*corr*"), x = names(aper_cat), value = TRUE)],
    aper_cat[, grep(pattern = glob2rx("*N_app*"), x = names(aper_cat), value = TRUE)],
    aper_cat[, grep(pattern = glob2rx("*1sig_local_depth_app*"), x = names(aper_cat), value = TRUE)]
  )
  PhotStats = data.frame(
    phot_cat[, grep(pattern = glob2rx("*1sig_local_depth"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*R50$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*R90$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*R100$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*N50$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*N90$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*N100$"), x = names(phot_cat), value = TRUE)],
    phot_cat[, grep(pattern = glob2rx("*maskfrac$"), x = names(phot_cat), value = TRUE)]
  )
  
  combined_photometry_tot = combined_cat[, grep(pattern = glob2rx("*fluxt$"), x = names(combined_cat), value = TRUE)]
  combined_photometry_col = combined_cat[, grep(pattern = glob2rx("*fluxc$"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a16 = combined_cat[, grep(pattern = glob2rx("*flux_app_0.16$"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a32 = combined_cat[, grep(pattern = glob2rx("*flux_app_0.32$"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a50 = combined_cat[, grep(pattern = glob2rx("*flux_app_0.5$"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a60 = combined_cat[, grep(pattern = glob2rx("*flux_app_0.6$"), x = names(combined_cat), value = TRUE)]
  
  combined_photometry_tot_err = combined_cat[, grep(pattern = glob2rx("*scaled_fluxt_err"), x = names(combined_cat), value = TRUE)]
  combined_photometry_col_err = combined_cat[, grep(pattern = glob2rx("*scaled_fluxc_err"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a16_err = combined_cat[, grep(pattern = glob2rx("*scaled_flux_err_app_0.16"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a32_err = combined_cat[, grep(pattern = glob2rx("*scaled_flux_err_app_0.32"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a50_err = combined_cat[, grep(pattern = glob2rx("*scaled_flux_err_app_0.5"), x = names(combined_cat), value = TRUE)]
  combined_photometry_a60_err = combined_cat[, grep(pattern = glob2rx("*scaled_flux_err_app_0.6"), x = names(combined_cat), value = TRUE)]
  
  ## filter heavily masked fluxes
  maskfrac_filters = phot_cat[, grep(pattern = glob2rx("*maskfrac$"), x = names(phot_cat), value = TRUE)]
  maskfrac_filters_idx = maskfrac_filters >= 0.5
  # ## flag unrealistically negative total fluxes
  negative_fluxes_idx = combined_photometry_tot + 4*combined_photometry_tot_err <= 0 ## 4sigma still negative
  negative_fluxes_idx[is.na(negative_fluxes_idx)] = TRUE 
  
  MW_foreground_extinction_func = function(phot_cat, err_cat){
    ext_corr = sweep(phot_cat, 2, MW_foreground_extinction$A_flux_inv[1:dim(phot_cat)[2]], `*`)
    ext_err_corr = sweep(err_cat, 2, MW_foreground_extinction$A_flux_inv[1:dim(err_cat)[2]], `*`)
    return(cbind(ext_corr, ext_err_corr))
  }
  
  ## MW foreground extinction
  if(fname == "smiles"){
    matched_photometry_ext_corr = MW_foreground_extinction_func(matched_photometry,matched_photometry_errs)
  }else{
    matched_photometry_ext_corr = NULL
    matched_redshifts = NULL
  }
  
  combined_cat$fluxt
  
  # aperture_correction_factor_temp_16 = combined_cat[, c("F560W_fluxt", "F770W_fluxt")]/combined_cat[, c("F560W_flux_app_0.16", "F770W_flux_app_0.16")]
  # # aperture_correction_factor_temp = combined_cat$F560W_fluxt
  # aperture_correction_factor_temp_16$F560W_fluxt[is.infinite(aperture_correction_factor_temp_16$F560W_fluxt) | aperture_correction_factor_temp_16$F560W_fluxt < 0] = NA
  # aperture_correction_factor_temp_16$F770W_fluxt[is.infinite(aperture_correction_factor_temp_16$F770W_fluxt) | aperture_correction_factor_temp_16$F770W_fluxt < 0] = NA
  # aperture_correction_factor_16 = pmax(1, rowMins(
  #   as.matrix(aperture_correction_factor_temp_16), na.rm = TRUE
  # ))
  # 
  # aperture_correction_factor_temp_32 = combined_cat[, c("F560W_fluxt", "F770W_fluxt")]/combined_cat[, c("F560W_flux_app_0.32", "F770W_flux_app_0.32")]
  # # aperture_correction_factor_temp = combined_cat$F560W_fluxt
  # aperture_correction_factor_temp_32$F560W_fluxt[is.infinite(aperture_correction_factor_temp_32$F560W_fluxt) | aperture_correction_factor_temp_32$F560W_fluxt < 0] = NA
  # aperture_correction_factor_temp_32$F770W_fluxt[is.infinite(aperture_correction_factor_temp_32$F770W_fluxt) | aperture_correction_factor_temp_32$F770W_fluxt < 0] = NA
  # aperture_correction_factor_32 = pmax(1, rowMins(
  #   as.matrix(aperture_correction_factor_temp_32), na.rm = TRUE
  # ))
  # 
  # aperture_correction_factor_temp_50 = combined_cat[, c("F560W_fluxt", "F770W_fluxt")]/combined_cat[, c("F560W_flux_app_0.5", "F770W_flux_app_0.5")]
  # # aperture_correction_factor_temp = combined_cat$F560W_fluxt
  # aperture_correction_factor_temp_50$F560W_fluxt[is.infinite(aperture_correction_factor_temp_50$F560W_fluxt) | aperture_correction_factor_temp_50$F560W_fluxt < 0] = NA
  # aperture_correction_factor_temp_50$F770W_fluxt[is.infinite(aperture_correction_factor_temp_50$F770W_fluxt) | aperture_correction_factor_temp_50$F770W_fluxt < 0] = NA
  # aperture_correction_factor_50 = pmax(1, rowMins(
  #   as.matrix(aperture_correction_factor_temp_50), na.rm = TRUE
  # ))
  # 
  # aperture_correction_factor_temp_60 = combined_cat[, c("F560W_fluxt", "F770W_fluxt")]/combined_cat[, c("F560W_flux_app_0.6", "F770W_flux_app_0.6")]
  # # aperture_correction_factor_temp = combined_cat$F560W_fluxt
  # aperture_correction_factor_temp_60$F560W_fluxt[is.infinite(aperture_correction_factor_temp_60$F560W_fluxt) | aperture_correction_factor_temp_60$F560W_fluxt < 0] = NA
  # aperture_correction_factor_temp_60$F770W_fluxt[is.infinite(aperture_correction_factor_temp_60$F770W_fluxt) | aperture_correction_factor_temp_60$F770W_fluxt < 0] = NA
  # aperture_correction_factor_60 = pmax(1, rowMins(
  #   as.matrix(aperture_correction_factor_temp_60), na.rm = TRUE
  # ))
  # 
  # app_tot_flux_corr = data.frame(
  #   "app_corr_16" = aperture_correction_factor_16, 
  #   "app_corr_32" = aperture_correction_factor_32, 
  #   "app_corr_50" = aperture_correction_factor_50, 
  #   "app_corr_60" = aperture_correction_factor_60
  # )

  combined_photometry_tot_ext_corr = MW_foreground_extinction_func(combined_photometry_tot,combined_photometry_tot_err)
  combined_photometry_col_ext_corr = MW_foreground_extinction_func(combined_photometry_col,combined_photometry_col_err)
  combined_photometry_a16_ext_corr = MW_foreground_extinction_func(combined_photometry_a16,combined_photometry_a16_err)
  combined_photometry_a32_ext_corr = MW_foreground_extinction_func(combined_photometry_a32,combined_photometry_a32_err)
  combined_photometry_a50_ext_corr = MW_foreground_extinction_func(combined_photometry_a50,combined_photometry_a50_err)
  combined_photometry_a60_ext_corr = MW_foreground_extinction_func(combined_photometry_a60,combined_photometry_a60_err)
  
  col_idx <- 1:(ncol(combined_photometry_tot))
  combined_photometry_tot_ext_corr[, 1:(ncol(combined_photometry_tot))][maskfrac_filters_idx] = -999
  combined_photometry_col_ext_corr[, 1:(ncol(combined_photometry_col))][maskfrac_filters_idx] = -999
  combined_photometry_a16_ext_corr[, 1:(ncol(combined_photometry_a16))][maskfrac_filters_idx] = -999
  combined_photometry_a32_ext_corr[, 1:(ncol(combined_photometry_a32))][maskfrac_filters_idx] = -999
  combined_photometry_a50_ext_corr[, 1:(ncol(combined_photometry_a50))][maskfrac_filters_idx] = -999
  combined_photometry_a60_ext_corr[, 1:(ncol(combined_photometry_a60))][maskfrac_filters_idx] = -999
  
  combined_photometry_tot_ext_corr[, 1:(ncol(combined_photometry_tot))][negative_fluxes_idx] = -999
  combined_photometry_col_ext_corr[, 1:(ncol(combined_photometry_col))][negative_fluxes_idx] = -999
  combined_photometry_a16_ext_corr[, 1:(ncol(combined_photometry_a16))][negative_fluxes_idx] = -999
  combined_photometry_a32_ext_corr[, 1:(ncol(combined_photometry_a32))][negative_fluxes_idx] = -999
  combined_photometry_a50_ext_corr[, 1:(ncol(combined_photometry_a50))][negative_fluxes_idx] = -999
  combined_photometry_a60_ext_corr[, 1:(ncol(combined_photometry_a60))][negative_fluxes_idx] = -999
  
  profound_cat$ID = 1:dim(profound_cat)[1]
  
  the_catalogue = list(
    "ProFound" = cbind( data.frame(profound_cat)[, !grepl("flux|aper|F.+W", names(profound_cat))], profound_cat[, c("flux", "flux_err")], "FIELD" = fname),
    
    "PhotTot" = combined_photometry_tot_ext_corr,
    "PhotCol" = combined_photometry_col_ext_corr,
    "Phot16" = combined_photometry_a16_ext_corr,
    "Phot32" = combined_photometry_a32_ext_corr,
    "Phot50" = combined_photometry_a50_ext_corr,
    "Phot60" = combined_photometry_a60_ext_corr,
    
    "RawTot" = cbind(combined_photometry_tot,combined_photometry_tot_err),
    "RawCol" = cbind(combined_photometry_col,combined_photometry_col_err),
    "Raw16" = cbind(combined_photometry_a16,combined_photometry_a16_err),
    "Raw32" = cbind(combined_photometry_a32,combined_photometry_a32_err),
    "Raw50" = cbind(combined_photometry_a50,combined_photometry_a50_err),
    "Raw60" = cbind(combined_photometry_a60,combined_photometry_a60_err),
    
    "PhotStats" = PhotStats, 
    "AperStats" = AperStats,
    # "AperTotCorr" = app_tot_flux_corr,
    
    "EpochsV2Phot" = matched_photometry_ext_corr,
    "EpochsV2z" = matched_redshifts
  )
  
  ## NA will become -999 in the catalogue
  Rfits_write(
    the_catalogue, 
    filename = paste0("~/Documents/Epochs/data/catalogues/", fname, "_profound_multiband_cat.fits")
  )
  
  return(
    the_catalogue
  )
}

field_names = c(
  "smiles",
  "midis"
)

smiles = do_stuff_goodss(fname = field_names[1])
midis = do_stuff_goodss(fname = field_names[2])

match_idx = coordmatch(
  coordref = midis$ProFound[, c("RAmax", "Decmax")],
  coordcompare = smiles$ProFound[, c("RAmax", "Decmax")],
  rad = 0.3
)
idx_not_in_smiles = 1:dim(midis$ProFound)[1] %notin% match_idx$bestmatch$refID # ids of the midis sample

and_in_the_darkness_bind_them <- Map(
  function(a, b) {
    # common <- intersect(names(a), names(b))
    bind_rows(a, b[idx_not_in_smiles, ])
  },
  smiles, midis
)

col_names <- c(
  "ProFound",
  "PhotTot",
  "PhotCol",
  "Phot16",
  "Phot32",
  "Phot50",
  "Phot60",
  "RawTot",
  "RawCol",
  "Raw16",
  "Raw32",
  "Raw50",
  "Raw60",
  "PhotStats",
  "AperStats"
)
Rfits_write(
  and_in_the_darkness_bind_them[col_names], 
  filename = paste0("~/Documents/Epochs/data/catalogues/", "miri", "_profound_multiband_cat.fits")
)

h5file <- path.expand("~/Documents/Epochs/data/analysis_data.h5")
# Create the file if it doesn't exist yet; otherwise we append into it
if (!file.exists(h5file)) {
  h5createFile(h5file)
}

extname = "Phot"
# Create the group (returns FALSE + doesn't error if it already exists)
h5createGroup(h5file, extname)

for (nm in col_names) {
  df <- as.data.frame(and_in_the_darkness_bind_them[[nm]])
  path <- paste0(extname, "/", nm)
  
  fid <- H5Fopen(h5file)
  exists <- H5Lexists(fid, path)
  H5Fclose(fid)
  
  if (exists) {
    h5delete(h5file, path)
  }
  
  h5write(df, file = h5file, name = path)
}

h5createGroup(h5file, "frame_stats")
for(ff in field_names){
  frame_stats = data.frame(
    fread(
      paste0("~/Documents/Epochs/data/depths/", ff, "_profound_frame_stats.csv")
    )
  )
  
  path <- paste0("frame_stats", "/", ff)
  
  fid <- H5Fopen(h5file)
  exists <- H5Lexists(fid, path)
  H5Fclose(fid)
  
  if (exists) {
    h5delete(h5file, path)
  }
  
  h5write(frame_stats, file = h5file, name = path)
}
