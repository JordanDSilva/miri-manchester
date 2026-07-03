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
library(segmented)

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
combined_cat = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_combined_catalogue.csv"))

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

# fit powerlaws with eazy normal on PDF, same as with prospect fitting
do_powerlaws_segments = function(kk, redshifts){
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
  flux_data$flux[flux_data$flux <= 0 | flux_data$flux == -999 | grepl("F.+M", filters_info$filters)] = NA
  # flux_data$flux[flux_data$flux <= 0 | flux_data$flux == -999] = NA
  flux_data$fluxerr = sqrt( as.numeric(df_flux_err[kk, ])^2 + (min_pc_err * flux_data$flux)^2 ) ## put 5pc minimum error
  
  ## convert flux and err to CGS Fnu
  flux_data$fluxFnu = Jansky2CGS(flux_data$flux)
  flux_data$fluxerrFnu = Jansky2CGS(flux_data$fluxerr)
  ## convert to Flambda
  flux_data$fluxFlambda = convert_freq2wave(flux_freq = flux_data$fluxFnu, wave = flux_data$cenwave)
  flux_data$fluxerrFlambda = convert_freq2wave(flux_freq = flux_data$fluxerrFnu, wave = flux_data$cenwave)
  
  # work in log-log space
  
  idx = !is.na(flux_data$fluxFlambda)
  x <- log10(flux_data$cenwave[idx])
  y <- log10(flux_data$fluxFlambda[idx])
  
  df_fit <- data.frame(x = x, y = y)
  
  # fit base linear model
  lm_fit <- lm(y ~ x, data = df_fit)
  
  # fit piecewise with 2 breakpoints (3 segments)
  seg_fit <- tryCatch({
    
    # make sure breakpoints are strictly inside data range
    x_min <- min(x) + 0.05
    x_max <- max(x) - 0.05
    
    bp_lo <- x_min + (x_max - x_min) * 0.33
    bp_hi <- x_min + (x_max - x_min) * 0.66
    
    segmented(
      lm_fit,
      seg.Z = ~x,
      npsi = 2
    )
    
  }, error = function(e) {
    cat("Galaxy", kk, "failed:", e$message, "\n")
    NULL
  }, warning = function(w) {
    cat("Galaxy", kk, "warning:", w$message, "\n")
    NULL   # treat warnings as failures too
  })
  
  if(is.null(seg_fit)){
    return(NULL)
  }else{
    waveout = seq(2, 9.35, length.out = 100000)
    yhat       <- predict(seg_fit, newdata = data.frame(x = waveout))
    
    bp           <- seg_fit$psi[, "Est."] ## break points
    slope_values = slope(seg_fit)$x ## slopes
    
    out = data.frame(
      galID     = kk,
      beta1     = slope_values[1, "Est."] - 2.0,
      beta2     = slope_values[2, "Est."] - 2.0,
      beta3     = slope_values[3, "Est."] - 2.0,
      beta1_err = slope_values[1, "St.Err."],
      beta2_err = slope_values[2, "St.Err."],
      beta3_err = slope_values[3, "St.Err."],
      break1    = 10^bp[1] / (1+z),
      break2    = 10^bp[2] / (1+z),
      RAmax = phot_cat$ProFound$RAmax[kk], 
      Decmax = phot_cat$ProFound$Decmax[kk], 
      z = z
    )
    
    pdf(paste0(main_dir, "/plots/powerlaw_segments/plot_gal_", kk, ".pdf"), width = 12, height = 8)
    par(mfrow = c(1,1), mar = c(2.5, 2.5, 0.5, 0.5), oma = rep(0, 4))
    magplot(
      flux_data$cenwave / (1+z),
      flux_data$fluxFlambda,
      log = "xy",
      xlim = c(100, 200000),
      ylim = c(min(flux_data$fluxFlambda[flux_data$fluxFlambda>0], na.rm = TRUE)*0.5, max(flux_data$fluxFlambda, na.rm = TRUE)*1.5),
      pch = 16,
      col = "red"
    )
    magerr(
      flux_data$cenwave / (1+z),
      flux_data$fluxFlambda,
      ylo = flux_data$fluxerrFlambda,
      col = "red"
    )
    lines(10^waveout/(1+z), 10^yhat, col = "blue", lwd = 2)
    dev.off()
    
    return(out)
  }
}

fit_idx = which(
  combined_cat$selected == 1L & 
    phot_cat$ProFound$flag == 0L & 
    # -2.5*log10(phot_cat$Phot32$F150W_flux_app_0.32/phot_cat$Phot32$F200W_flux_app_0.32) < 1.0 & 
    # -2.5*log10(phot_cat$Phot32$F277W_flux_app_0.32/phot_cat$Phot32$F444W_flux_app_0.32) > 0.5 &
    phot_cat$Raw32$F444W_flux_app_0.32/phot_cat$Raw16$F444W_flux_app_0.16 <= 1.7
)

Rimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F444W_smiles_warp_smiles.fits")
Gimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F200W_smiles_warp_smiles.fits")
Bimage = Rfits_point("~/mnt/morgan/raid/data/JP/ProFound/Data/smiles/smiles/JADES-DR3-GS_F150W_smiles_warp_smiles.fits")
stack = Rfits_read("~/Documents/Epochs/data/images/smiles/smiles_smiles_smiles_profound_stack.fits")

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

redshifts = combined_cat$zQ50
registerDoParallel(cores = ncores)
beta_powerlaw_df = foreach(kk = fit_idx, .combine = rbind) %do% {
  print(kk)
  fit = do_powerlaws_segments(kk, redshifts)
  return(fit)
}
stopImplicitCluster()

power_law_cat <- data.frame(
  galID     = integer(),
  beta1     = numeric(),
  beta2     = numeric(),
  beta3     = numeric(),
  beta1_err = numeric(),
  beta2_err = numeric(),
  beta3_err = numeric(),
  break1    = numeric(),
  break2    = numeric()
)
for (kk in fit_idx) {
  fit = do_powerlaws_segments(kk, redshifts)
  power_law_cat[kk, ] <- fit
}


for(i in 1:1000){
  dev.off()
}
