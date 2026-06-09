library(ProSpect)
library(Highlander)
library(celestial)
library(data.table)
library(foreach)
library(doParallel)
library(magicaxis)
library(matrixStats)

set.seed(666)

# main_dir = "/scratch/"

CMANiters = 500
LDNiters = 1000
ncores = 1 ## will set to 1 for SLURM job array 

args <- commandArgs(trailingOnly = TRUE)
if (length(args[1]) == 0 | is.na(args[1])) {
  idx <- 1
} else if (grepl(":", args[1])) {
  bounds <- as.integer(strsplit(args[1], ":")[[1]])
  idx <- bounds[1]:bounds[2]
} else {
  idx <- as.integer(args[1])
}

if (length(args[2]) == 0 | is.na(args[2])) {
  main_dir = "~/Documents/Epochs/"
} else {
  main_dir <- args[2]
}

combined_catalogue = data.frame(fread(
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

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

# fit powerlaws with eazy normal on PDF, same as with prospect fitting
do_powerlaws = function(kk){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs_normals[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs_normals$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs_normals$zgrid[which.max(zPDF_load)]
  
  ## convert flux and err to CGS Fnu
  flux_data$fluxFnu = Jansky2CGS(flux_data$flux)
  flux_data$fluxerrFnu = Jansky2CGS(flux_data$fluxerr)
  ## convert to Flambda
  flux_data$fluxFlambda = convert_freq2wave(flux_freq = flux_data$fluxFnu, wave = flux_data$cenwave)
  flux_data$fluxerrFlamda = convert_freq2wave(flux_freq = flux_data$fluxerrFnu, wave = flux_data$cenwave)

  powerlaw_model = function(parm){
    
    z = 10**parm[1]
    norm_BB = 10**parm[2] ## two power laws should be equal at this point ?
    betaUV = parm[3]
    betaOpt = parm[4]
    x_BB = 10^parm[5] # wavelength at the break, angstrom
    
    x = 10**seq(2, 9.35, 0.01) ## rest wavelength
    
    # x_BB = (3645) ## angstrom
    smoothstep = 0.5 * (tanh((x-x_BB)) + 1)
    fluxUV = x**betaUV
    fluxOpt = x**betaOpt
    normUV = approx(x = x, y = fluxUV, xout = x_BB)$y
    normOpt = approx(x = x, y = fluxOpt, xout = x_BB)$y 
    
    rest_frame = fluxUV + smoothstep*(fluxOpt*normUV/normOpt - fluxUV) # per angstrom
    normRest = approx(x = x, y = rest_frame, xout = x_BB)$y
    
    observed_wavelength = x * (1+z)
    observed_frame = rest_frame * norm_BB/normRest
    observed_frame_Fnu = CGS2Jansky(convert_wave2freq(
      observed_frame, 
      wave = observed_wavelength
    ))
    
    ## observed IGM
    # igm = Inoue14_IGM(
    #   wave = observed_wavelength, 
    #   z = z, Inoue14_LAFcoef = Inoue14_LAFcoef, 
    #   Inoue14_DLAcoef = Inoue14_DLAcoef
    # )
    
    df = list(
      wave = observed_wavelength,
      flux = observed_frame,
      fluxJy = observed_frame_Fnu
    )
    
    return(df)
  }
  LL = function(parm, Data){
    yy_model = Data$model(parm)
    
    phot = photom_flux(
      yy_model$wave, 
      yy_model$flux,
      outtype = "Jy",
      filters = Data$filter_transforms
    )
    
    idx = Data$flux_data$cenwave >= 1216*(1+10**parm[1]) & Data$flux_data$cenwave <= 6562.79*(1+10**parm[1]) & Data$flux_data$flux > 0
    
    # print(Data$flux_data$cenwave[idx])
    
    likelihood = sum(
      dnorm(
        x = log10(Data$flux_data$flux[idx]), 
        mean = log10(phot[idx]),
        sd = Data$flux_data$fluxerr[idx]/(Data$flux_data$flux[idx] * log(10)),
        log = TRUE
      ),
      na.rm = TRUE
    ) + Data$prior(parm)
    return( likelihood )
  }
  
  highout = Highlander(
    parm = c(log10(zeazy), -21, -2, 2, log10(3646)),
    Data = list(
      flux_data = flux_data, 
      filter_transforms = filtout, 
      model = powerlaw_model, 
      prior = function(parm){
        zPDF_func(10**parm[1]) +
          dnorm(10^parm[5], mean = 3646, sd = 100, log = TRUE) +
          0
      }
    ),
    likefunc = LL,
    liketype = "max",
    lower = c(-2, -25, -10, -10, log10(2500)),
    upper = c(1.5, -16, 10, 10, log10(4500)), 
    Niters = c(500,1000),
    optim_iters = 2
  )

  # par(mfrow = c(1,2), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))
  magplot(
    flux_data$cenwave/1e4,
    flux_data$flux,
    log = "xy",
    xlim = c(0.0002, 30),
    ylim = c(1e-11, max(flux_data$flux, na.rm = TRUE)*1.5),
    pch = 16,
    col = "red"
  )
  magerr(
    flux_data$cenwave/1e4,
    flux_data$flux,
    ylo = flux_data$fluxerr,
    col = "red"
  )
  lines(
    powerlaw_model(highout$parm)$wave/1e4,
    powerlaw_model(highout$parm)$fluxJy,
    col = "black", lwd = 3
  )
  abline(
    v = 1216*(1 + 10**highout$parm[1])/1e4, lty = 2, lwd = 2
  )
  abline(
    v = 3646*(1 + 10**highout$parm[1])/1e4, lty = 3, lwd = 2
  )
  abline(
    v = 6562.79*(1 + 10**highout$parm[1])/1e4, lty = 4, lwd = 2
  )
  
  # magplot(
  #   flux_data$cenwave/1e4,
  #   flux_data$fluxFlambda/10**highout$parm[2],
  #   log = "x",
  #   xlim = c(0.1, 30),
  #   ylim = c(0, max(flux_data$fluxFlambda/10**highout$parm[2], na.rm = TRUE)*1.5),
  #   pch = 16,
  #   col = "red"
  # )
  # magerr(
  #   flux_data$cenwave/1e4,
  #   flux_data$fluxFlambda/10**highout$parm[2],
  #   ylo = flux_data$fluxerrFlamda/10**highout$parm[2],
  #   col = "red"
  # )
  # lines(
  #   powerlaw_model(highout$parm)$wave/1e4,
  #   powerlaw_model(highout$parm)$flux/10**highout$parm[2],
  #   col = "black", lwd = 3
  # )
  # abline(
  #   v = 3646*(1 + 10**highout$parm[1])/1e4, lty = 2, lwd = 2
  # )
  # abline(
  #   v = 6562.79*(1 + 10**highout$parm[1])/1e4, lty = 3, lwd = 2
  # )
  
  out = list(
    "bestfit" = list(
      "parm" = highout$parm, 
      "FinalFlux" = data.frame("wave" = powerlaw_model(highout$parm)$wave, "flux" = powerlaw_model(highout$parm)$fluxJy),
      "flux_data" = flux_data
      ),
    "highout" = highout,
    "save_stub" = save_stub
  )
  return(out)
}
do_powerlaws_zfix = function(kk, redshifts){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  z = redshifts[kk]
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  ## convert flux and err to CGS Fnu
  flux_data$fluxFnu = Jansky2CGS(flux_data$flux)
  flux_data$fluxerrFnu = Jansky2CGS(flux_data$fluxerr)
  ## convert to Flambda
  flux_data$fluxFlambda = convert_freq2wave(flux_freq = flux_data$fluxFnu, wave = flux_data$cenwave)
  flux_data$fluxerrFlamda = convert_freq2wave(flux_freq = flux_data$fluxerrFnu, wave = flux_data$cenwave)
  
  powerlaw_model_fixz = function(parm, z){

    norm_BB = 10**parm[1] ## two power laws should be equal at this point ?
    betaUV = parm[2]
    betaOpt = parm[3]
    x_BB = 10^parm[4] # wavelength at the break, angstrom
    
    x = 10**seq(2, 9.35, 0.01) ## rest wavelength
    
    # x_BB = (3645) ## angstrom
    smoothstep = 0.5 * (tanh((x-x_BB)) + 1)
    fluxUV = x**betaUV
    fluxOpt = x**betaOpt
    normUV = approx(x = x, y = fluxUV, xout = x_BB)$y
    normOpt = approx(x = x, y = fluxOpt, xout = x_BB)$y 
    
    rest_frame = fluxUV + smoothstep*(fluxOpt*normUV/normOpt - fluxUV) # per angstrom
    normRest = approx(x = x, y = rest_frame, xout = x_BB)$y
    
    observed_wavelength = x * (1+z)
    observed_frame = rest_frame * norm_BB/normRest
    observed_frame_Fnu = CGS2Jansky(convert_wave2freq(
      observed_frame, 
      wave = observed_wavelength
    ))
    
    ## observed IGM
    # igm = Inoue14_IGM(
    #   wave = observed_wavelength, 
    #   z = z, Inoue14_LAFcoef = Inoue14_LAFcoef, 
    #   Inoue14_DLAcoef = Inoue14_DLAcoef
    # )

    df = list(
      wave = observed_wavelength,
      flux = observed_frame,
      fluxJy = observed_frame_Fnu
    )
    
    return(df)
  }
  LL = function(parm, Data){
    yy_model = Data$model(parm, Data$z)
    
    phot = photom_flux(
      yy_model$wave, 
      yy_model$flux,
      outtype = "Jy",
      filters = Data$filter_transforms
    )
    
    idx = Data$flux_data$cenwave >= (1215.7)*(1+Data$z) & Data$flux_data$flux > 0
    
    likelihood = sum(
      dnorm(
        x = log10(Data$flux_data$flux[idx]), 
        mean = log10(phot[idx]),
        sd = Data$flux_data$fluxerr[idx]/(log(10) * Data$flux_data$flux[idx]),
        log = TRUE
      ),
      na.rm = TRUE
    ) + Data$prior(parm)
    return( likelihood )
  }
  
  highout = Highlander(
    parm = c(-21, -2, 2, log10(3646)),
    Data = list(
      flux_data = flux_data, 
      filter_transforms = filtout, 
      model = powerlaw_model_fixz, 
      prior = function(parm){
          dnorm(10^parm[4], mean = 3646, sd = 100, log = TRUE) +
          0
      },
      # constraints = beta_constraint_func,
      z = z
    ),
    
    likefunc = LL,
    liketype = "max",
    lower = c(-25, -10, -10, log10(2500)),
    upper = c(-16, 10, 10, log10(4500)), 
    Niters = c(500,1000),
    optim_iters = 2
  )
  
  # par(mfrow = c(1,2), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))
  # magplot(
  #   flux_data$cenwave / (1+z),
  #   flux_data$flux/1e-6,
  #   log = "xy",
  #   xlim = c(700, 20000),
  #   ylim = c(min(flux_data$flux[flux_data$flux>0], na.rm = TRUE)*0.5, max(flux_data$flux, na.rm = TRUE)*1.5)/1e-6,
  #   pch = 16,
  #   col = "red"
  # )
  # magerr(
  #   flux_data$cenwave / (1+z),
  #   flux_data$flux/1e-6,
  #   ylo = flux_data$fluxerr/1e-6,
  #   col = "red"
  # )
  # lines(
  #   powerlaw_model_fixz(highout$parm, z = z)$wave/(1+z),
  #   powerlaw_model_fixz(highout$parm, z = z)$fluxJy/1e-6,
  #   col = "black", lwd = 3
  # )
  # abline(
  #   v = 3646, lty = 2, lwd = 2
  # )
  
  # magplot(
  #   flux_data$cenwave/1e4,
  #   flux_data$fluxFlambda/10**highout$parm[1],
  #   log = "x",
  #   xlim = c(0.1, 30),
  #   ylim = c(0, max(flux_data$fluxFlambda/10**highout$parm[1], na.rm = TRUE)*1.5),
  #   pch = 16,
  #   col = "red"
  # )
  # magerr(
  #   flux_data$cenwave/1e4,
  #   flux_data$fluxFlambda/10**highout$parm[1],
  #   ylo = flux_data$fluxerrFlamda/10**highout$parm[1],
  #   col = "red"
  # )
  # lines(
  #   powerlaw_model_fixz(highout$parm, z)$wave/1e4,
  #   powerlaw_model_fixz(highout$parm, z)$flux/10**highout$parm[1],
  #   col = "black", lwd = 3
  # )
  # abline(
  #   v = 3646*(1 + z)/1e4, lty = 2, lwd = 2
  # )
  # abline(
  #   v = 6562.79*(1 + z)/1e4, lty = 3, lwd = 2
  # )
  
  out = list(
    "bestfit" = list(
      "parm" = highout$parm, 
      "FinalFlux" = data.frame("wave" = powerlaw_model_fixz(highout$parm, z)$wave, "flux" = powerlaw_model_fixz(highout$parm, z)$fluxJy)
    ),
    "highout" = highout,
    "save_stub" = save_stub
  )
  return(out)
}
do_powerlaws_kocevski = function(kk, redshifts){
  message("Running gal: ", kk)
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  z = redshifts[kk]
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  ## convert flux and err to CGS Fnu
  flux_data$fluxFnu = Jansky2CGS(flux_data$flux)
  flux_data$fluxerrFnu = Jansky2CGS(flux_data$fluxerr)
  ## convert to Flambda
  flux_data$fluxFlambda = convert_freq2wave(flux_freq = flux_data$fluxFnu, wave = flux_data$cenwave)
  flux_data$fluxerrFlamda = convert_freq2wave(flux_freq = flux_data$fluxerrFnu, wave = flux_data$cenwave)
  
  # Table 2 selection in Kocevski+25
  # if(z >= 2 & z < 3.25){
  #   ff_idxUV = filters_info$filters %in% c("F606W", "F814W", "F115W")
  # }else if(z >= 3.25 & z < 4.75){
  #   ff_idxUV = filters_info$filters %in% c("F814W", "F115W", "F150W")
  # }else if(z >= 4.75 & z < 8){
  #   ff_idxUV = filters_info$filters %in% c("F115W", "F150W", "F200W")
  # }else{
  #   ff_idxUV = filters_info$filters %in% c("F150W", "F200W", "F277W")
  # }
  # 
  # # Table 2 selection in Kocevski+25
  # if(z >= 2 & z < 3.25){
  #   ff_idxOpt = filters_info$filters %in% c("F150W", "F200W", "F277W")
  # }else if(z >= 3.25 & z < 4.75){
  #   ff_idxOpt = filters_info$filters %in% c("F200W", "F277W", "F356W")
  # }else if(z >= 4.75 & z < 8){
  #   ff_idxOpt = filters_info$filters %in% c("F277W", "F356W", "F444W")
  # }else{
  #   ff_idxOpt = filters_info$filters %in% c("F356W", "F444W")
  # }
  
  ff_idxUV = filters_info$cenwave/(1+z) > 1216 & filters_info$cenwave/(1+z) <= 3000 
  ff_idxOpt = filters_info$cenwave/(1+z) > 3000 & filters_info$cenwave/(1+z) <= 5000
  
  powerlaw_model_kocevski = function(parm, z){
    
    norm = 10**parm[1]
    beta = parm[2]
    
    x = 10**seq(2, 9.35, 0.01) ## rest wavelength
    
    flux = x**beta

    observed_wavelength = x * (1+z)
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
    
    likeUV = sum(
      dnorm(
        x = log10(Data$flux_data$flux[idxUV]), 
        mean = log10(photUV[idxUV]),
        sd = Data$flux_data$fluxerr[idxUV]/(log(10) * Data$flux_data$flux[idxUV]),
        log = TRUE
      ), na.rm = TRUE
    )
    likeOpt = sum(
      dnorm(
        x = log10(Data$flux_data$flux[idxOpt]), 
        mean = log10(photOpt[idxOpt]),
        sd = Data$flux_data$fluxerr[idxOpt]/(log(10) * Data$flux_data$flux[idxOpt]),
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
  
  # highout$parm
  # colQuantiles(as.matrix(highout$LD_last$Posterior1), probs =0.5)
  
  pdf(paste0(main_dir, "/data/prospect/power_laws/plots/plot_gal_", kk, ".pdf"), width = 12, height = 8)
  par(mfrow = c(1,2), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))
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
  abline(
    v = 3646, lty = 2, lwd = 2
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

# fit_idx = which(combined_catalogue$zQ50 >= 2 & combined_catalogue$FitType == "StellarAGN")
fit_idx = which(combined_catalogue$FitType != "")

registerDoParallel(cores = 8)
beta_powerlaw_df = foreach(kk = fit_idx, .combine = rbind) %dopar% {
  
  save_stub = paste0(
    main_dir, "/data/prospect/power_laws/galID_",kk, ".rds" 
  )
  
  if(!file.exists(save_stub)){
    fit = do_powerlaws_kocevski(kk, combined_catalogue$zQ50)
    
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

beta_powerlaw_df = data.frame(fread(
  "~/Documents/Epochs/data/prospect/power_laws_astro_cat.csv"
))
magplot(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaUVQ50,
  xlim = c(0, 20),
  # ylim = c(-5, 5)
)
magerr(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaUVQ50,
  ylo = beta_powerlaw_df$betaUV_lo,
  yhi = beta_powerlaw_df$betaUV_hi
)
magplot(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaOptQ50,
  xlim = c(0, 20),
  # ylim = c(-5, 5)
)
magerr(
  beta_powerlaw_df$zQ50,
  beta_powerlaw_df$betaOptQ50,
  ylo = beta_powerlaw_df$betaOpt_lo,
  yhi = beta_powerlaw_df$betaOpt_hi
)
# 
lrd_idx = which(beta_powerlaw_df$betaUVQ50 <= -0.37 & beta_powerlaw_df$betaOptQ50 > 0 & beta_powerlaw_df$betaUVQ50 > -2.8)
# pdf("~/Desktop/beta_powerlaws.pdf", width = 8, height = 5)
magplot(
  beta_powerlaw_df$betaUVQ50,
  beta_powerlaw_df$betaOptQ50,
  xlim = c(-5, 2),
  ylim = c(-4, 4),
  xlab = "Beta UV slope",
  ylab = "Beta Optical slope"
)
points(
  beta_powerlaw_df$betaUVQ50[lrd_idx],
  beta_powerlaw_df$betaOptQ50[lrd_idx],
  col = "red",
  pch = 16
)
# legend(
#   x = "topright",
#   legend = c(
#     paste0("Better fit with AGN models"), 
#     paste0("Kocevski+25 LRD criterion")
#   ),
#   pch = c(1, 16),
#   col = c("black", "red")
# )
# dev.off()
# 

Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F444W_smiles_warp_smiles.fits")
Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F200W_smiles_warp_smiles.fits")
Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F150W_smiles_warp_smiles.fits")

stack = Rfits_read("~/Documents/Epochs/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")

R = Rimage[
  beta_powerlaw_df$RA[lrd_idx][kk],
  beta_powerlaw_df$Dec[lrd_idx][kk],
  box = 100,
  type = "coord"
]
pix_coords = Rwcs_s2p(beta_powerlaw_df$RA[lrd_idx][kk], beta_powerlaw_df$Dec[lrd_idx][kk], keyvalues = Rimage$keyvalues)
segim = stack$segim[pix_coords[1],
                    pix_coords[2],
                    box = 100]$imDat


foo=profoundAperPhot(
  image = R, 
  segim = segim,
  app_diam = 0.16,
  fluxtype = "Jansky",
  pixscale = pixscale(R),
  magzero = R$keyvalues$MAGZERO, 
  correction = TRUE
)
foo2=profoundAperPhot(
  image = R, 
  segim = segim,
  app_diam = 0.60,
  fluxtype = "Jansky",
  pixscale = pixscale(R),
  magzero = R$keyvalues$MAGZERO, 
  correction = TRUE
)

par(mfrow = c(1,2), oma = rep(1.5, 4), mar = rep(1.5,4))
kk = 2
Rwcs_imageRGB(
  R = Rimage[
    beta_powerlaw_df$RA[lrd_idx][kk],
    beta_powerlaw_df$Dec[lrd_idx][kk],
    box = 100,
    type = "coord"
  ],
  G = Gimage[
    beta_powerlaw_df$RA[lrd_idx][kk],
    beta_powerlaw_df$Dec[lrd_idx][kk],
    box = 100,
    type = "coord"
  ],
  B = Bimage[
    beta_powerlaw_df$RA[lrd_idx][kk],
    beta_powerlaw_df$Dec[lrd_idx][kk],
    box = 100,
    type = "coord"
  ]
)
sed_pl = readRDS(
  paste0("~/Documents/Epochs/data/prospect/power_laws/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
)
sed_agn = readRDS(
  paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
)
sed_stars = readRDS(
  paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", beta_powerlaw_df$ID[lrd_idx][kk], ".rds")
)

magplot(
  filters_info$cenwave,
  unlist(df_flux[beta_powerlaw_df$ID[lrd_idx][kk], ]),
  ylim = c(1e-11, 1e-4),
  log = "xy"
)
magerr(
  filters_info$cenwave,
  unlist(df_flux[beta_powerlaw_df$ID[lrd_idx][kk], ]),
  ylo = unlist(df_flux_err[beta_powerlaw_df$ID[lrd_idx][kk], ])
)
lines(
  sed_pl$bestfit$FinalFluxUV$wave,
  sed_pl$bestfit$FinalFluxUV$flux
)
lines(
  sed_pl$bestfit$FinalFluxOpt$wave,
  sed_pl$bestfit$FinalFluxOpt$flux
)
lines(
  sed_stars$bestfit$SEDout$FinalFlux$wave,
  sed_stars$bestfit$SEDout$FinalFlux$flux, 
  col = "blue"
)
lines(
  sed_agn$bestfit$SEDout$FinalFlux$wave,
  sed_agn$bestfit$SEDout$FinalFlux$flux, 
  col = "red"
)
legend(
  x = "bottomright",
  legend = c(
    paste0("z = ", round(beta_powerlaw_df$zQ50[lrd_idx][kk], 3))
  )
)


kocevski25 = read.table("~/Documents/Epochs/data/literature/AGN-Catalogues/kocevski25_LRDs.txt", skip = 35)
match_kocevski = coordmatch(
  coordref = kocevski25[,3:4],
  coordcompare = beta_powerlaw_df[, c("RA", "Dec")],
  rad = 0.3
)

magplot(
  kocevski25$V5[match_kocevski$bestmatch$refID],
  beta_powerlaw_df$zQ50[match_kocevski$bestmatch$compareID]
)
abline(0,1)
# 
# kk = 4
# RA = kocevski25$V3[matchbox_faggot$bestmatch$refID][kk]
# Dec = kocevski25$V4[matchbox_faggot$bestmatch$refID][kk]
# Rwcs_imageRGB(
#   R = Rimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ],
#   G = Gimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ],
#   B = Bimage[
#     RA,
#     Dec, 
#     box = 100,
#     type = "coord"
#   ]
# )
# 
# foo = readRDS(
#   paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", combined_catalogue$ID[matchbox_faggot$bestmatch$compareID][kk], ".rds")
# )
# magplot(
#   filters_info$cenwave,
#   unlist(df_flux[combined_catalogue$ID[matchbox_faggot$bestmatch$compareID][kk], ]), 
#   ylim = c(1e-8, 1e-6), 
#   log = "xy"
# )
# lines(
#   foo$bestfit$FinalFluxUV$wave,
#   foo$bestfit$FinalFluxUV$flux
# )
# lines(
#   foo$bestfit$FinalFluxOpt$wave,
#   foo$bestfit$FinalFluxOpt$flux
# )
# lines(
#   foo$bestfit$SEDout$FinalFlux
# )
# plot(foo$bestfit, ylim = c(1e-9, 1e-5))


magplot(
  phot_cat$Phot16$F444W_flux_app_0.16*phot_cat$AperStats$F444W_ap_correction_app_0.16,
  phot_cat$Phot60$F444W_flux_app_0.6*phot_cat$AperStats$F444W_ap_correction_app_0.6, 
  xlim = c(1e-11, 2),
  ylim = c(1e-11, 2),
  log = 'xy'
)
