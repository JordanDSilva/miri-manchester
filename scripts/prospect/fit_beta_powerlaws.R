library(ProSpect)
library(Highlander)
library(celestial)
library(data.table)
library(foreach)
library(doParallel)
library(magicaxis)

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

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)
## Frame info, filter info, throughputs loaded 
## prior and constraints loaded
## eazy zPDFs loaded

dir.create(
  paste0(main_dir, "/prospect/power_laws/"), 
  recursive = TRUE
)

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_0.32")]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

do_powerlaws = function(kk){
  message("Running gal: ", kk)
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$cenwave,
    "flux" = as.numeric(df_flux[kk, ])
    # "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (0.05 * flux_data$flux)^2 ) ## put 5pc minimum error
  
  filtout = filter_transforms
  
  zPDF_load = eazy_zPDFs[[paste0('gal',kk)]]
  zPDF_load[zPDF_load <= 0] = 1e-100
  zPDF_func = approxfun(
    x = eazy_zPDFs$zgrid,
    y = log(zPDF_load),
    rule = 2
  )## check z logged
  
  zeazy = eazy_zPDFs$zgrid[which.max(zPDF_load)]
  
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
    
    x = 10**seq(2, 9.35, 0.01) ## rest wavelength

    x_BB = (3645) ## angstrom
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
    igm = Inoue14_IGM(
      wave = observed_wavelength, 
      z = z, Inoue14_LAFcoef = Inoue14_LAFcoef, 
      Inoue14_DLAcoef = Inoue14_DLAcoef
    )
    
    # magplot(
    #   x, fluxUV, log = "xy", ylim = c(1e30, 1e40), xlim = c(2e3, 2e5)
    # )
    # lines(
    #   x, fluxOpt * normUV/normOpt
    # )
    # lines(
    #   x, rest_frame, col = "red"
    # )
    
    # magplot(
    #   observed_wavelength,
    #   observed_frame, 
    #   type = "l", log = "xy", 
    #   xlim = c(10, 3e5), 
    #   ylim = c(1e-32, 1e-16)
    # )
    # lines(
    #   observed_wavelength,
    #   observed_frame * igm, 
    #   col = "blue"
    # )
    # abline(v = x_BB * (1+z))
    # 
    # 
    # magplot(
    #   observed_wavelength,
    #   observed_frame_Fnu, 
    #   type = "l", log = "xy", 
    #   xlim = c(10, 3e5), 
    #   ylim = c(1e-11, 1e-6)
    # )
    # lines(
    #   observed_wavelength,
    #   observed_frame_Fnu * igm, 
    #   col = "blue"
    # )
    # abline(v = x_BB * (1+z))

    df = list(
      wave = observed_wavelength,
      flux = observed_frame*igm,
      fluxJy = observed_frame_Fnu*igm
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
    
    idx = Data$flux_data$cenwave <= (65602.79)*(1+10**parm[1])
    
    # print(Data$flux_data$cenwave[idx])
    
    likelihood = sum(
      dnorm(
        x = Data$flux_data$flux[idx], 
        mean = phot[idx],
        sd = Data$flux_data$fluxerr[idx],
        log = TRUE
      ),
      na.rm = TRUE
    ) + Data$prior(parm)
    return( likelihood )
  }
  
  highout = Highlander(
    parm = c(log10(zeazy), -21, -2, 2),
    Data = list(
      flux_data = flux_data, 
      filter_transforms = filtout, 
      model = powerlaw_model, 
      prior = function(parm){
        zPDF_func(10**parm[1]) +
          0
      }
    ),
    likefunc = LL,
    liketype = "max",
    lower = c(-2, -25, -10, -10),
    upper = c(1.5, -16, 10, 10), 
    Niters = c(500,1000),
    optim_iters = 2
  )

  par(mfrow = c(1,2), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))
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
    v = 3646*(1 + 10**highout$parm[1])/1e4, lty = 2, lwd = 2
  )
  abline(
    v = 6562.79*(1 + 10**highout$parm[1])/1e4, lty = 3, lwd = 2
  )
  
  magplot(
    flux_data$cenwave/1e4,
    flux_data$fluxFlambda/10**highout$parm[2], 
    log = "x", 
    xlim = c(0.1, 30), 
    ylim = c(0, max(flux_data$fluxFlambda/10**highout$parm[2], na.rm = TRUE)*1.5),
    pch = 16, 
    col = "red"
  )
  magerr(
    flux_data$cenwave/1e4,
    flux_data$fluxFlambda/10**highout$parm[2],
    ylo = flux_data$fluxerrFlamda/10**highout$parm[2], 
    col = "red"
  )
  lines(
    powerlaw_model(highout$parm)$wave/1e4,
    powerlaw_model(highout$parm)$flux/10**highout$parm[2], 
    col = "black", lwd = 3
  )
  abline(
    v = 3646*(1 + 10**highout$parm[1])/1e4, lty = 2, lwd = 2
  )
  abline(
    v = 6562.79*(1 + 10**highout$parm[1])/1e4, lty = 3, lwd = 2
  )

}

red2 =
  -2.5*log10(df_flux$F115W_flux_app_0.32/df_flux$F150W_flux_app_0.32) < 0.8 & -2.5*log10(df_flux$F200W_flux_app_0.32/df_flux$F277W_flux_app_0.32) > 0.7 & -2.5*log10(df_flux$F200W_flux_app_0.32/df_flux$F356W_flux_app_0.32) > 1.0 & 
  phot_cat$Phot50$F444W_flux_app_0.5/phot_cat$Phot16$F444W_flux_app_0.16 < 1.7 & df_flux$F444W_flux_app_0.32/df_flux_err$F444W_scaled_flux_err_app_0.32 > 4

which(red2)

do_powerlaws(3)
