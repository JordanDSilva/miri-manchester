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

phot_cat = Rfits_read(paste0(main_dir, "/data/catalogues/miri_profound_multiband_cat.fits"))

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

smiles_hlsp = Rfits_read(
  "~/Documents/Epochs/data/catalogues/hlsp_smiles_jwst_miri_goodss_multi_v1.0_catalog.fits"
)

test_f560w_ee_data = data.frame(fread("~/Documents/Epochs/data/depths/ee_f560w_gordon25.csv"))
test_f560w_ee_corr_fun = approxfun(
  x = pmax(0, test_f560w_ee_data$r),
  y = test_f560w_ee_data$ee,
  yleft = 0,
  yright = 1
)

match_idx = coordmatch(
  coordref = phot_cat$ProFound[, c("RAmax", "Decmax")],
  coordcompare = goodss_wfc3_cat[, c("RA", "DEC")],
  rad = 0.3
)

mag_mag_plot = function(mag1, mag2, maglim = 30, main = ""){
  
  mag_diff = mag2[mag1 <= maglim] - mag1[mag1 <= maglim]
  
  med  <- median(mag_diff, na.rm = TRUE)
  p16  <- quantile(mag_diff, 0.16, na.rm = TRUE)
  p84  <- quantile(mag_diff, 0.84, na.rm = TRUE)
  
  layout(matrix(c(1,2), 1, 2), widths = c(3,1))
  par(mar = c(5, 4, 4, 0.5))
  magplot( 
    mag1[mag1 <= maglim],
    mag2[mag1 <= maglim],
    xlim = c(18, 32),
    ylim = c(18, 32),
    xlab = "magx [mag]",
    ylab = "magy [mag]",
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
    xlab = "magy - magx"
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
  -2.5*log10(phot_cat$PhotTot$F356W_fluxt[match_idx$bestmatch$refID])+8.9,
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

pdf("~/Documents/Epochs/plots/mag-mag_compare_32.pdf", width = 12, height = 9)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F356W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH1_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F356W / IRAC CH1"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F444W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH2_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F444W / IRAC CH2"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F560W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH3_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F560W / IRAC CH3"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F770W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$IRAC_CH4_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F770W / IRAC CH4"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Raw32$F435W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$ACS_F435W_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F435W / F435W"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F150W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$WFC3_F160W_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F150W / F160W"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Phot32$F200W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(goodss_wfc3_cat$HAWKI_KS_FLUX[match_idx$bestmatch$compareID])+23.9, 
  main = "F200W / HWAK-K Ks"
)
dev.off()


match_idx = coordmatch(
  coordref = goodss_wfc3_cat[, c("RA", "DEC")],
  coordcompare = smiles_hlsp$CIRC[, c("RA", "DEC")],
  rad = 0.3
)
mag_mag_plot(
  -2.5*log10(goodss_wfc3_cat$IRAC_CH3_FLUX[match_idx$bestmatch$refID])+23.9,
  -2.5*log10(smiles_hlsp$KRON$F560W_KRON[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F560W / IRAC CH3"
)
mag_mag_plot(
  -2.5*log10(goodss_wfc3_cat$IRAC_CH4_FLUX[match_idx$bestmatch$refID])+23.9,
  -2.5*log10(smiles_hlsp$KRON$F770W_KRON[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F560W / IRAC CH3"
)


match_idx = coordmatch(
  coordref = phot_cat$ProFound[, c("RAmax", "Decmax")],
  smiles_hlsp$KRON[, c("RA", "DEC")],
  rad = 0.3
)
mag_mag_plot(
  -2.5*log10(phot_cat$RawTot$F560W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(smiles_hlsp$KRON$F560W_KRON[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F560W - Me vs. SMILES"
)
mag_mag_plot(
  -2.5*log10(phot_cat$Raw32$F560W_flux_app_0.32[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(smiles_hlsp$CIRC$F560W_CIRC4[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F560W - Me vs. SMILES"
)

mag_mag_plot(
  -2.5*log10(phot_cat$RawTot$F770W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(smiles_hlsp$KRON$F770W_KRON[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F770W / IRAC CH4"
)
mag_mag_plot(
  -2.5*log10(phot_cat$RawTot$F2550W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(smiles_hlsp$KRON$F2550W_KRON[match_idx$bestmatch$compareID] * 1e-9)+8.9, 
  main = "F770W / IRAC CH4"
)

GS_epochsV2 = data.frame(fread(file = "~/Documents/Epochs/data/catalogues/EPOCHS-v2_GS_Sel-F277W+F356W+F444W_v14_psfmatch_F444W_empirical+v14_psfmatch_F444W_empirical.csv"))

match_idx = coordmatch(
  coordref = phot_cat$ProFound[, c("RAmax", "Decmax")],
  coordcompare = GS_epochsV2[, c("ALPHA_J2000", "DELTA_J2000")],
  rad = 0.3
)
mag_mag_plot(
  -2.5*log10(phot_cat$PhotTot$F356W_fluxt[match_idx$bestmatch$refID])+8.9,
  -2.5*log10(GS_epochsV2$FLUX_AUTO_F356W[match_idx$bestmatch$compareID])+(-6.10 - 2.5*log10(2.29e-14)), 
  main = "F356W / IRAC CH1"
)

match_idx = coordmatch(
  coordref = goodss_wfc3_cat[, c("RA", "DEC")],
  coordcompare = GS_epochsV2[, c("ALPHA_J2000", "DELTA_J2000")],
  rad = 0.3
)
mag_mag_plot(
  -2.5*log10(goodss_wfc3_cat$IRAC_CH1_FLUX[match_idx$bestmatch$refID])+23.9,
  # -2.5*log10(GS_epochsV2$FLUX_AUTO_F356W[match_idx$bestmatch$compareID])+(-6.10 - 2.5*log10(2.29e-14)), 
  -2.5*log10(GS_epochsV2$FLUX_AUTO_F356W[match_idx$bestmatch$compareID])+(-6.10 - 2.5*log10(2.29e-14)), 
  main = "F356W / IRAC CH1"
)
