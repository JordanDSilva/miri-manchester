library(ProSpect)
library(celestial)
library(foreach)
library(matrixStats)
library(Rfits)
library(data.table)
library(magicaxis)
library(scales)
library(RColorBrewer)

set.seed(666)

main_dir = "~/Documents/Epochs/"

frame_info = data.frame(fread(paste0(main_dir, "data/depths/smiles_profound_frame_stats.csv")))
filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))

phot_cat = Rfits_read(paste0(main_dir, "/data/catalogues/smiles_profound_multiband_cat.fits"))

# extract column names from the header
lines <- readLines(paste0(main_dir, "/data/catalogues/hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1-1photom_cat.txt"))
header_lines <- lines[grepl("^#", lines)]
# parse column names from "# N  COLNAME ..." format
col_names <- sub("^#\\s+\\d+\\s+(\\S+).*", "\\1", header_lines)
goodss_wfc3_cat = data.frame(
  fread(
    paste0(main_dir, "/data/catalogues/hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1-1photom_cat.txt"), header = FALSE, skip = 73
  )
)
setnames(goodss_wfc3_cat, col_names)

match_idx = coordmatch(
  coordref = phot_cat$ProFound[, c("RAmax", "Decmax")],
  goodss_wfc3_cat[, c("RA", "DEC")],
  rad = 0.3
)

mag_mag_plot = function(mag1, mag2, main = ""){
  
  mag_diff = mag2[mag1 <= 100] - mag1[mag1 <= 100]
  
  med  <- median(mag_diff, na.rm = TRUE)
  p16  <- quantile(mag_diff, 0.16, na.rm = TRUE)
  p84  <- quantile(mag_diff, 0.84, na.rm = TRUE)
  
  layout(matrix(c(1,2), 1, 2), widths = c(3,1))
  par(mar = c(5, 4, 4, 0.5))
  magplot( 
    mag1,
    mag2,
    xlim = c(18, 32),
    ylim = c(18, 32),
    xlab = "ProFound+SMILES [mag]",
    ylab = "Guo+13 GOODSS [mag]",
    main = main
  )
  abline(0,1, col = "red")
  par(mar = c(5, 4, 4, 0.5))
  hh = maghist(
    mag_diff,
    breaks = seq(-2, 2, 0.1), 
    log = "y",
    xlim = c(-2,2),
    ylim = c(1, 1200),
    xlab = "Guo13 - ProFound"
  )

  abline(v = c(p16, med, p84), 
         col = c("blue", "red", "blue"), 
         lty = c(2, 1, 2), lwd = 2)
  legend(
    x = "topright",
    legend = paste0(c("Q16 ", "Q50 ", "Q84 "), round(c(p16, med, p84), 3))
  )
  layout(1)
}

pdf("~/Documents/Epochs/plots/mag-mag_compare.pdf", width = 12, height = 9)
mag_mag_plot(
  -2.5*log10(phot_cat$RawTot$F356W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH1_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F356W / IRAC CH1"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F444W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH2_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F444W / IRAC CH2"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F560W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH3_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F560W / IRAC CH3"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F770W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH4_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F770W / IRAC CH4"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F435W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$ACS_F435W_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F435W / F435W"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F150W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$WFC3_F160W_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F150W / F160W"
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F200W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$HAWKI_KS_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F200W / HWAK-K Ks"
)
dev.off()

# GS_epochsV2 = fread(
#   file = "~/Documents/Epochs/data/EPOCHS-v2_GS_Sel-F277W+F356W+F444W_v14_psfmatch_F444W_empirical+v14_psfmatch_F444W_empirical.csv"
# )
# match_idx_epochsV2 = coordmatch(
#   coordref = GS_epochsV2[, c("ALPHA_J2000", "DELTA_J2000")],
#   goodss_wfc3_cat[, c("RA", "DEC")],
#   rad = 0.3
# )
# mag_mag_plot(
#   # -2.5*log10(phot_cat$PhotTot$F150W_flux_tot[match_idx$bestmatch$refID])+8.9,
#   -2.5*log10(GS_epochsV2$FLUX_AUTO_F150W[match_idx_epochsV2$bestmatch$refID])+28.01473,
#   -2.5*log10(goodss_wfc3_cat$WFC3_F160W_FLUX[match_idx_epochsV2$bestmatch$compareID])+23.9, 
#   main = "F150W / F160W"
# )

## number counts 
stone25 = df3 <- data.frame(
  Flux_uJy = c(3.98e-2, 6.31e-2, 1.00e-1, 1.58e-1, 2.51e-1, 3.98e-1, 6.31e-1,
               1.00e+0, 1.58e+0, 2.51e+0, 3.98e+0, 6.31e+0, 1.00e+1, 1.58e+1,
               2.51e+1, 3.98e+1, 6.31e+1, 1.00e+2, 1.58e+2, 2.51e+2, 3.98e+2,
               6.31e+2, 1.00e+3),
  
  F560W       = c(NA, NA, NA, NA, 71.3, 52.2, 38.7, 30.0, 23.0, 18.0, 13.7,  9.73,  6.60,  4.37, 2.46, 1.27, 0.46, 0.12, 0.06,   NA,   NA,   NA,   NA),
  F560W_err   = c(NA, NA, NA, NA,  1.4,  1.2,  1.1,  0.9,  0.8,  0.7,  0.6,  0.53,  0.44,  0.35, 0.27, 0.19, 0.12, 0.06, 0.04,   NA,   NA,   NA,   NA),
  
  F770W_SMILES     = c(NA, NA, NA, NA, 55.0, 41.7, 32.6, 25.3, 20.6, 16.1, 11.7,  8.48,  5.41,  3.17, 1.79, 0.85, 0.35, 0.21, 0.18, 0.09, 0.06, 0.03,   NA),
  F770W_SMILES_err = c(NA, NA, NA, NA,  1.3,  1.1,  1.0,  0.9,  0.8,  0.7,  0.6,  0.50,  0.40,  0.30, 0.23, 0.16, 0.10, 0.08, 0.07, 0.05, 0.04, 0.03,   NA),
  
  F770W_JADES     = c(218, 167, 130, 104, 81.8, 64.5, 50.8, 41.2, 31.6, 23.6, 17.6, 12.5,  8.01,  4.62, 2.06, 1.13, 0.31, 0.10,   NA,   NA,   NA,   NA,   NA),
  F770W_JADES_err = c(  5,   4,   4,   3,  2.9,  2.6,  2.3,  2.0,  1.8,  1.6,  1.3,  1.1,  0.90,  0.69, 0.46, 0.34, 0.18, 0.10,   NA,   NA,   NA,   NA,   NA),
  
  F1000W     = c(NA, NA, NA, NA, NA, 41.3, 31.1, 25.0, 20.0, 15.8, 12.1,  8.81,  5.96,  4.16, 2.92, 1.86, 1.27, 0.83, 0.53, 0.24, 0.21, 0.12, 0.06),
  F1000W_err = c(NA, NA, NA, NA, NA,  1.1,  1.0,  0.9,  0.8,  0.7,  0.6,  0.51,  0.42,  0.35, 0.29, 0.23, 0.19, 0.16, 0.12, 0.08, 0.08, 0.06, 0.04),
  
  F1280W     = c(NA, NA, NA, NA, NA, NA, 35.6, 27.1, 22.0, 18.0, 14.5, 11.2,  8.67,  6.17, 4.79, 3.58, 2.53, 1.73, 0.91, 0.56, 0.26, 0.09, 0.06),
  F1280W_err = c(NA, NA, NA, NA, NA, NA,  1.0,  0.9,  0.8,  0.7,  0.6,  0.6,  0.50,  0.42, 0.37, 0.32, 0.27, 0.22, 0.16, 0.13, 0.09, 0.05, 0.04),
  
  F1500W     = c(NA, NA, NA, NA, NA, NA, 39.3, 30.0, 24.1, 19.6, 15.6, 12.7,  9.96,  7.81, 5.99, 4.41, 3.03, 1.91, 0.94, 0.53, 0.21, 0.09, 0.03),
  F1500W_err = c(NA, NA, NA, NA, NA, NA,  1.0,  0.9,  0.8,  0.8,  0.7,  0.6,  0.54,  0.48, 0.42, 0.36, 0.30, 0.24, 0.17, 0.12, 0.08, 0.05, 0.03),
  
  F1800W     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 22.7, 18.5, 15.0, 12.6,  9.58, 7.26, 5.58, 4.03, 2.58, 1.56, 0.68, 0.29, 0.12, 0.03),
  F1800W_err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA,  0.8,  0.7,  0.7,  0.6,  0.53, 0.46, 0.40, 0.34, 0.27, 0.21, 0.14, 0.09, 0.06, 0.03),
  
  F2100W     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 27.7, 21.0, 16.4, 13.5, 10.7,  8.29, 6.23, 4.12, 2.79, 1.47, 0.79, 0.35, 0.09, 0.03),
  F2100W_err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA,  0.9,  0.8,  0.7,  0.6,  0.6,  0.49, 0.43, 0.35, 0.28, 0.21, 0.15, 0.10, 0.05, 0.03),
  
  F2550W     = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 13.5,  9.13, 6.50, 4.41, 2.70, 1.53, 0.76, 0.32, 0.18, 0.06),
  F2550W_err = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  0.6,  0.51, 0.43, 0.36, 0.28, 0.21, 0.15, 0.10, 0.07, 0.04)
)

stone25_counts = data.frame(
  stone25$F560W,
  stone25$F770W_SMILES,
  stone25$F1000W, 
  stone25$F1280W, 
  stone25$F1500W, 
  stone25$F1800W, 
  stone25$F2100W, 
  stone25$F2550W
)
calc_number_counts = function(mag, area, mag_bins, mag_bw, limiting_mag = 100){
  mag[is.na(mag) | mag == 0 | is.infinite(mag)] = 100
  
  hh = maghist(
    mag, mag_bins, plot = FALSE
  )
  
  df = data.frame(
    "count" = hh$counts / (area * mag_bw),
    "err" = sqrt(hh$counts) / (area * mag_bw),
    "cumul" = cumsum(hh$counts) / (area),
    "cumulerr" = sqrt(cumsum(hh$counts)) / (area)
  )
  
  df$count[hh$mids >= limiting_mag] = 0
  df$err[hh$mids >= limiting_mag] = 0
  df$cumul[hh$mids >= limiting_mag] = 0
  df$cumulerr[hh$mids >= limiting_mag] = 0
  
  return(df)
}

mag_bw = 0.5
mag_bins = seq(1, 50, mag_bw)
mag_mids = head(mag_bins,-1) + mag_bw*0.5

# tompkins25_f160w = fread("~/Desktop/digitized_data.csv")
# magplot(
#   tompkins25_f160w$mag,
#   tompkins25_f160w$N,
#   type = "l", lwd = 2
# )
# n_F150W = calc_number_counts(
#   mag = -2.5*log10(phot_cat$PhotTot$F150W_flux_tot)+8.9,
#   area = frame_info$JWST_MIRI_F560W_effarea_amin2/3600,
#   mag_bins = mag_bins,
#   mag_bw = mag_bw, 
#   limiting_mag = 23
# )
# points(
#   mag_mids, 
#   log10(n_F150W$count)
# )

MIRI_area = c(frame_info$F560W_effarea_amin2, frame_info$F770W_effarea_amin2, frame_info$F1000W_effarea_amin2, frame_info$F1280W_effarea_amin2,
              frame_info$F1500W_effarea_amin2, frame_info$F1800W_effarea_amin2, frame_info$F2100W_effarea_amin2, frame_info$F2550W_effarea_amin2)
MIRI_limiting_mag = -2.5*log10(c(frame_info$F560W_5sigma_depth_32asec, frame_info$F770W_5sigma_depth_32asec, frame_info$F1000W_5sigma_depth_32asec, frame_info$F1280W_5sigma_depth_32asec,
              frame_info$F1500W_5sigma_depth_32asec, frame_info$F1800W_5sigma_depth_32asec, frame_info$F2100W_5sigma_depth_32asec, frame_info$F2550W_5sigma_depth_32asec))+8.9
kk = 1
MIRI_counts = lapply(
  list(phot_cat$PhotTot$F560W_fluxt, phot_cat$PhotTot$F770W_fluxt, phot_cat$PhotTot$F1000W_fluxt, phot_cat$PhotTot$F1280W_fluxt, 
       phot_cat$PhotTot$F1500W_fluxt, phot_cat$PhotTot$F1800W_fluxt, phot_cat$PhotTot$F2100W_fluxt, phot_cat$PhotTot$F2550W_fluxt),
  function(x){
    out = calc_number_counts(
      mag = -2.5*log10(x)+8.9,
      area = MIRI_area[kk],
      mag_bins = mag_bins,
      mag_bw = mag_bw, 
      limiting_mag = MIRI_limiting_mag[kk]
    )
    kk <<- kk + 1
    return(out)
  }
)

cols <- c(
  "#2166AC",  # deep blue
  "#4393C3",  # mid blue
  "#92C5DE",  # light blue
  "#FDDBC7",  # pale yellow-orange
  "#F4A736",  # yellow-orange
  "#E8601C",  # orange
  "#D6301F",  # red-orange
  "#A50026"   # deep red
)

MIRI_names = filters_info$filters[24:31]
stone25_counts = list(
  stone25$F560W, stone25$F770W_SMILES, stone25$F1000W, stone25$F1280W,
  stone25$F1500W, stone25$F1800W, stone25$F2100W, stone25$F2550W
)
stone25_errs = list(
  stone25$F560W_err, stone25$F770W_SMILES_err, stone25$F1000W_err, stone25$F1280W_err,
  stone25$F1500W_err, stone25$F1800W_err, stone25$F2100W_err, stone25$F2550W_err
)

pdf("~/Documents/Epochs/plots/miri_counts.pdf", width = 7, height = 7)
par(mfrow = c(4,2), mar = c(3.5, 3.5, 0.5, 0.5), bg = "grey10", 
    col.axis = "white",   # tick labels
    col.lab  = "white",   # axis titles
    fg       = "white" 
)
# tick marks and plot border)
for(i in 1:8){
  magplot(
    NA,
    log = "xy",
    xlim = 10^( -0.4 * (c(26.5, 16) - 23.9)),
    ylim = c(1e-3, 3e2),
    xlab = "Flux [microJy]",
    ylab = "N(>)",
    grid = FALSE
  )
  lines(
    stone25$Flux_uJy,
    stone25_counts[[i]],
    col = cols[i],
    lty = 1
  )
  lines(
    stone25$Flux_uJy,
    stone25_counts[[i]] - stone25_errs[[i]],
    col = cols[i],
    lty = 3
  )
  lines(
    stone25$Flux_uJy,
    stone25_counts[[i]] + stone25_errs[[i]],
    col = cols[i],
    lty = 3
  )
  
  points(
    10^(-0.4 * (mag_mids - 8.9)) * 1e6,
    MIRI_counts[[i]]$cumul, 
    col = cols[i]
  )
  magerr(
    10^(-0.4 * (mag_mids - 23.9)),
    MIRI_counts[[i]]$cumul,
    ylo = MIRI_counts[[i]]$cumulerr,
    col = cols[i]
  )
  
  abline(
    v = 10^(-0.4 * (MIRI_limiting_mag[i] - 8.9))/1e-6, 
    col = cols[i]
  )
  legend(
    x = "bottomleft", 
    lty = c(1, NA),
    pch = c(NA, 1),
    legend = c("Stone+25 (SMILES)", paste0("ProFound = ", MIRI_names[i])),
    col = cols[i],
    cex = 1.1
  )
}
dev.off()


df_flux = as.data.frame(phot_cat$Phot32)[, paste0(filters_info$filters, "_flux_app_0.32")] 
df_flux_err = as.data.frame(phot_cat$Phot32)[, paste0(filters_info$filters, "_scaled_flux_err_app_0.32")]

# df_flux = df_flux[phot_cat$ProFound$flag == 0L, ]
# df_flux_err = df_flux_err[phot_cat$ProFound$flag == 0L, ]

# df_flux_err[df_flux == -999] = NA
# df_flux[df_flux == -999] = NA

weighted_stack = foreach(i = 1:31, .combine = c) %do% {
  
  ff = df_flux[,i]
  ff_err = df_flux_err[,i]
  
  weighted.mean(
    x = ff[ff > 0], 
    w = ff_err[ff > 0]^-2,
    na.rm = TRUE
  )
  
}


filter_area = unlist(as.data.frame(frame_info[, paste0(filters_info$filters, "_effarea_amin2")]))
# stacked_sed = colSums(as.matrix(df_flux*df_flux_err^-2), na.rm = TRUE) / colSums(as.matrix(df_flux_err), na.rm = TRUE)
stacked_sed = colMedians(as.matrix(df_flux), na.rm = TRUE)
# stacked_sed = colMeans(as.matrix(df_flux), na.rm = TRUE)
magplot(
  filters_info$cenwave,
  stacked_sed,
  log = "xy", 
  ylim = c(1e-9, 1e-5)
)

kk = 3006
profound_cat$flag[kk]
magplot(
  filters_info$cenwave,
  unlist(df_flux[kk,]),
  log = "xy"
)
magerr(
  filters_info$cenwave,
  unlist(df_flux[kk,]),
  ylo = unlist(df_flux_err[kk,])
)

magbin(
  -2.5*log10(phot_cat$Phot32$F560W_flux_app_0.32)+8.9,
  -2.5*log10(phot_cat$PhotTot$F560W_fluxt/phot_cat$Phot32$F560W_flux_app_0.32), 
  xlim = c(15, 30),
  ylim = c(-1, 1)
)
magbin(
  -2.5*log10(phot_cat$Phot32$F770W_flux_app_0.32)+8.9,
  -2.5*log10(phot_cat$PhotTot$F770W_fluxt/phot_cat$Phot32$F770W_flux_app_0.32), 
  xlim = c(15, 30),
  ylim = c(-1, 1)
)

