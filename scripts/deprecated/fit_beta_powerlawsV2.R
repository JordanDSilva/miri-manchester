library(ProSpect)
library(Highlander)
library(celestial)
library(data.table)
library(foreach)
library(doParallel)
library(magicaxis)
library(matrixStats)
library(Rfits)
library(Rwcs)

set.seed(666)

# main_dir = "/scratch/"

CMANiters = 500
LDNiters = 1000
ncores = 8 ## will set to 1 for SLURM job array 
min_pc_err = 0.1

main_dir = "~/Documents/Epochs/"
fname = "smiles"
eazy = Rfits_read_table(
  paste0(main_dir, "/data/eazy/", fname, "_agn_blue_sfhz", "/", fname, "_zphot", "_agn_blue_sfhz", ".out.zout.fits")
)

combined_cat = data.frame(fread(
  file = paste0(main_dir, "data/prospect/photoZ_bc03_snorm_fritz_combined_catalogue.csv")
))

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)
## Frame info, filter info, throughputs loaded 
## prior and constraints loaded
## eazy zPDFs loaded

dir.create(
  paste0(main_dir, "/data/prospect/power_laws/"), 
  recursive = TRUE
)
dir.create(
  paste0(main_dir, "/data/prospect/power_laws/plots/"), 
  recursive = TRUE
)

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F444W_smiles_warp_smiles.fits")
Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F200W_smiles_warp_smiles.fits")
Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F150W_smiles_warp_smiles.fits")
stack = Rfits_read("~/Documents/Epochs/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")

# fit powerlaws with eazy normal on PDF, same as with prospect fitting
do_powerlaws_kocevski = function(kk, redshifts){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  if(combined_cat$FitType[kk] == "StellarAGN"){
    sed = readRDS(
      paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", combined_cat$ID[kk], ".rds")
    )
  }else{
    sed = readRDS(
      paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", combined_cat$ID[kk], ".rds")
    )
  }
  
  z = redshifts[kk]
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (min_pc_err * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  ## convert flux and err to CGS Fnu
  flux_data$fluxFnu = Jansky2CGS(flux_data$flux)
  flux_data$fluxerrFnu = Jansky2CGS(flux_data$fluxerr)
  ## convert to Flambda
  flux_data$fluxFlambda = convert_freq2wave(flux_freq = flux_data$fluxFnu, wave = flux_data$cenwave)
  flux_data$fluxerrFlamda = convert_freq2wave(flux_freq = flux_data$fluxerrFnu, wave = flux_data$cenwave)
  
  ## Table 2 selection in Kocevski+25
  if(z >= 2 & z < 3.25){
    ff_idxUV = filters_info$filters %in% c("F606W", "F814W", "F115W")
    ff_idxOpt = filters_info$filters %in% c("F150W", "F200W", "F277W")
  }else if(z >= 3.25 & z < 4.75){
    ff_idxUV = filters_info$filters %in% c("F814W", "F115W", "F150W")
    ff_idxOpt = filters_info$filters %in% c("F200W", "F277W", "F356W")
  }else if(z >= 4.75 & z < 8){
    ff_idxUV = filters_info$filters %in% c("F115W", "F150W", "F200W")
    ff_idxOpt = filters_info$filters %in% c("F277W", "F356W", "F444W")
  }else{
    ff_idxUV = filters_info$filters %in% c("F150W", "F200W", "F277W")
    ff_idxOpt = filters_info$filters %in% c("F356W", "F444W")
  }
  
  powerlaw_model_kocevski = function(parm, z){
    
    norm = 10**parm[1]
    beta = parm[2]
    
    x = 10**seq(2, 9.35, 0.01) ## rest wavelength
    
    flux = x**beta
    
    observed_wavelength = x
    # observed_frame = rest_frame * norm_BB/normRest
    observed_frame = flux * norm
    observed_frame_Fnu = CGS2Jansky(convert_wave2freq(
      observed_frame, 
      wave = observed_wavelength
    ))
    
    df = list(
      wave = observed_wavelength,
      flux = observed_frame,
      fluxJy = observed_frame_Fnu
    )
    
    return(df)
  }
  
  LL = function(parm, Data){
    
    yy_modelUV = Data$model(parm[1:2], Data$z)
    photUV = photom_flux(
      yy_modelUV$wave, 
      yy_modelUV$flux,
      outtype = "Jy",
      filters = Data$filter_transforms
    )
    
    yy_modelOpt = Data$model(parm[3:4], Data$z)
    photOpt = photom_flux(
      yy_modelOpt$wave, 
      yy_modelOpt$flux,
      outtype = "Jy",
      filters = Data$filter_transforms
    )
    
    idxUV = Data$ff_idxUV & Data$flux_data$flux > 0
    idxOpt = Data$ff_idxOpt & Data$flux_data$flux > 0
    
    # likeUV = sum(
    #   dnorm(
    #     x = log10(Data$flux_data$flux[idxUV]), 
    #     mean = log10(photUV[idxUV]),
    #     sd = Data$flux_data$fluxerr[idxUV]/(log(10) * Data$flux_data$flux[idxUV]),
    #     log = TRUE
    #   ), na.rm = TRUE
    # )
    # likeOpt = sum(
    #   dnorm(
    #     x = log10(Data$flux_data$flux[idxOpt]), 
    #     mean = log10(photOpt[idxOpt]),
    #     sd = Data$flux_data$fluxerr[idxOpt]/(log(10) * Data$flux_data$flux[idxOpt]),
    #     log = TRUE
    #   ), na.rm = TRUE
    # )
    
    likeUV = sum(
      dnorm(
        x = (Data$flux_data$flux[idxUV]), 
        mean = (photUV[idxUV]),
        sd = Data$flux_data$fluxerr[idxUV],
        log = TRUE
      ), na.rm = TRUE
    )
    likeOpt = sum(
      dnorm(
        x = (Data$flux_data$flux[idxOpt]), 
        mean = (photOpt[idxOpt]),
        sd = Data$flux_data$fluxerr[idxOpt],
        log = TRUE
      ), na.rm = TRUE
    )
    
    likelihood = likeUV + likeOpt + Data$prior(parm)
    return( likelihood )
  }
  
  highout = Highlander(
    parm = c(-21, 0, -21, 0),
    Data = list(
      flux_data = flux_data, 
      filter_transforms = filtout, 
      model = powerlaw_model_kocevski,
      prior = function(parm){0},
      z = z,
      ff_idxUV = ff_idxUV,
      ff_idxOpt = ff_idxOpt
    ),
    likefunc = LL,
    liketype = "max",
    lower = c(-35, -10, -35, -10),
    upper = c(10, 10, 10, 10), 
    Niters = c(500,1000),
    optim_iters = 2
  )
  
  png(paste0(main_dir, "/plots/power_laws/plot_gal_", kk, ".png"), width = 20, height = 10, units = "in", res = 240)
  par(mfrow = c(1,4), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))

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
  xy_coords = Rwcs_s2p(
    phot_cat$ProFound$RAmax[kk],
    phot_cat$ProFound$Decmax[kk],
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
  
  magplot(
    flux_data$cenwave / (1+z),
    flux_data$flux/1e-6,
    log = "xy",
    xlim = c(700, 20000),
    ylim = c(min(flux_data$flux[flux_data$flux>0], na.rm = TRUE)*0.5, max(flux_data$flux, na.rm = TRUE)*1.5)/1e-6,
    pch = 16,
    col = "red"
  )
  points(
    flux_data$cenwave[ff_idxUV] / (1+z),
    flux_data$flux[ff_idxUV]/1e-6,
    cex = 1.5,
    col = "blue"
  )
  points(
    flux_data$cenwave[ff_idxOpt] / (1+z),
    flux_data$flux[ff_idxOpt]/1e-6,
    cex = 1.5,
    col = "purple"
  )
  magerr(
    flux_data$cenwave / (1+z),
    flux_data$flux/1e-6,
    ylo = flux_data$fluxerr/1e-6,
    col = "red"
  )
  lines(
    powerlaw_model_kocevski(highout$parm[1:2], z = z)$wave/(1+z),
    powerlaw_model_kocevski(highout$parm[1:2], z = z)$fluxJy/1e-6,
    col = "black", lwd = 3
  )
  lines(
    powerlaw_model_kocevski(highout$parm[3:4], z = z)$wave/(1+z),
    powerlaw_model_kocevski(highout$parm[3:4], z = z)$fluxJy/1e-6,
    col = "black", lwd = 3
  )
  lines(
    sed$bestfit$SEDout$FinalFlux$wave / (1+z),
    sed$bestfit$SEDout$FinalFlux$flux/1e-6, 
    col = "grey", lwd = 2
  )
  abline(
    v = 3646, lty = 2, lwd = 2
  )
  legend(
    x = "bottomright",
    legend = c(
      paste0("z = ", round(z, 3)), 
      paste0("ID = ", kk),
      paste0("Fittype = ", combined_cat$FitType[kk])
    )
  )
  
  magplot(
    flux_data$cenwave / (1+z),
    flux_data$fluxFlambda/10**highout$parm[1],
    log = "x",
    xlim = c(700, 20000),
    ylim = c(0, max(flux_data$fluxFlambda/10**highout$parm[1], na.rm = TRUE)*1.5),
    pch = 16,
    col = "red"
  )
  points(
    flux_data$cenwave[ff_idxUV] / (1+z),
    flux_data$fluxFlambda[ff_idxUV]/10**highout$parm[1],
    cex = 1.5,
    col = "blue"
  )
  points(
    flux_data$cenwave[ff_idxOpt] / (1+z),
    flux_data$fluxFlambda[ff_idxOpt]/10**highout$parm[1],
    cex = 1.5,
    col = "purple"
  )
  magerr(
    flux_data$cenwave / (1+z),
    flux_data$fluxFlambda/10**highout$parm[1],
    ylo = flux_data$fluxerrFlamda/10**highout$parm[1],
    col = "red"
  )
  lines(
    powerlaw_model_kocevski(highout$parm[1:2], z)$wave / (1+z),
    powerlaw_model_kocevski(highout$parm[1:2], z)$flux/10**highout$parm[1],
    col = "black", lwd = 3
  )
  lines(
    powerlaw_model_kocevski(highout$parm[3:4], z)$wave / (1+z),
    powerlaw_model_kocevski(highout$parm[3:4], z)$flux/10**highout$parm[1],
    col = "black", lwd = 3
  )
  abline(
    v = 3646, lty = 2, lwd = 2
  )
  dev.off()
  
  out = list(
    "galInfo" = list(ID = kk, RAmax = phot_cat$ProFound$RAmax[kk], Decmax = phot_cat$ProFound$Decmax[kk], z = z),
    "bestfit" = list(
      "parm" = highout$parm, 
      "FinalFluxUV" = data.frame("wave" = powerlaw_model_kocevski(highout$parm[1:2], z)$wave, "flux" = powerlaw_model_kocevski(highout$parm[1:2], z)$fluxJy),
      "FinalFluxOpt" = data.frame("wave" = powerlaw_model_kocevski(highout$parm[3:4], z)$wave, "flux" = powerlaw_model_kocevski(highout$parm[3:4], z)$fluxJy)
    ),
    "highout" = highout,
    "save_stub" = save_stub
  )
  return(out)
}

fit_idx = which(
  combined_cat$selected == 1L & 
    phot_cat$ProFound$flag == 0 & 
    # phot_cat$Raw32$F444W_flux_app_0.32/phot_cat$Raw16$F444W_flux_app_0.16 <= 1.7 & 
    combined_cat$zQ50 >= 2.0
)

registerDoParallel(cores = 8)
beta_powerlaw_df = foreach(kk = fit_idx, .combine = rbind) %dopar% {
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  redshifts = combined_cat$zQ50
  
  if(!file.exists(save_stub)){
    fit = do_powerlaws_kocevski(kk, redshifts)
    
    astro_quantiles = colQuantiles(as.matrix(fit$highout$LD_last$Posterior1), probs = c(0.5, 0.16, 0.84))
    params = c("normUV", "betaUV", "normOpt", "betaOpt")
    rownames(astro_quantiles) = params
    
    best = fit$bestfit$parm
    names(best) = params
    
    summary_dt = list(
      best, 
      astro_quantiles[, "50%"], 
      astro_quantiles[, "16%"], 
      astro_quantiles[, "84%"], 
      astro_quantiles[, "50%"] - 
        astro_quantiles[, "16%"], 
      astro_quantiles[, "84%"] - 
        astro_quantiles[, "50%"]
    )
    Qnames = c("Best", "Q50", "Q16", "Q84", "_lo", "_hi")
    for (x in seq_along(summary_dt)) {
      names(summary_dt[[x]]) = paste0(
        names(summary_dt[[x]]), 
        Qnames[x]
      )
    }
    summary_dt_nice = as.data.frame(t(c(unlist(summary_dt))))
    
    fit$astro = summary_dt_nice
    
    saveRDS(
      fit,
      fit$save_stub
    )
  }
  return(NULL)
}
stopImplicitCluster()


fit_idx = which(
  combined_cat$selected == 1L & 
    phot_cat$ProFound$flag == 0 & 
    
    phot_cat$Raw32$F444W_flux_app_0.32/phot_cat$Raw16$F444W_flux_app_0.16 <= 1.7 & 
    combined_cat$zQ50 >= 2.0
)

for(kk in fit_idx){
  pdf(
    paste0(
      "~/Documents/Epochs/plots/lrd_colour_select//gal_", kk, ".pdf"
    ), width = 18, height = 6
  )
  par(mfrow = c(1,3), oma = rep(1.5, 4), mar = rep(1.5,4))
  
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
  xy_coords = Rwcs_s2p(
    phot_cat$ProFound$RAmax[kk],
    phot_cat$ProFound$Decmax[kk],
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
  
  
  magplot(
    filters_info$cenwave,
    unlist(df_flux[kk, ]),
    ylim = c(1e-10, 1e-5),
    xlim = c(1000, 1e6),
    log = "xy"
  )
  magerr(
    filters_info$cenwave,
    unlist(df_flux[kk, ]),
    ylo = unlist(df_flux_err[kk, ]),
  )
  
  legend(
    x = "bottomright",
    legend = c(
      paste0("z = ", round(eazy$z_ml[kk], 3)), 
      paste0("ID = ", kk)
    )
  )
  dev.off()
}
