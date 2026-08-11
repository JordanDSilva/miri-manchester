library(magicaxis)
library(data.table)
library(Rfits)
library(celestial)
library(Rwcs)
library(rhdf5)

## number counts

main_dir = "~/Documents/Epochs/"
h5file <- path.expand("~/Documents/Epochs/data/analysis_data.h5")

field_names = c(
  "smiles",
  "midis"
)

mag_bw = 0.5
mag_bins = seq(1, 50, mag_bw)
mag_mids = head(mag_bins,-1) + mag_bw*0.5

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))

df_flux_name = "PhotTot"
flux_col_names = paste0(filters_info$filters, "_fluxt")


.calc_number_counts = function(mag, area, mag_bins, mag_bw, limiting_mag = 100){
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

number_counts = function(df_flux, limiting_mags, area){
  MIRI_counts = Map(
    function(x, a, lim) {
      .calc_number_counts(
        mag = -2.5*log10(x) + 8.9,
        area = a,
        mag_bins = mag_bins,
        mag_bw = mag_bw,
        limiting_mag = lim
      )
    },
    df_flux,
    area,
    limiting_mags
  )
  
  names(MIRI_counts) = filters_info$filters
  n_df = data.frame(do.call(cbind, MIRI_counts))
  n_df$mag = mag_mids
  n_df$uJy_mids = 10^(-0.4 * (mag_mids - 23.9))
  names(n_df) = gsub("\\.", "_", names(n_df))
  
  return(
    n_df
  )
}

safe_select <- function(df, cols) {
  df = data.frame(df)
  out <- data.frame(matrix(-999, nrow = nrow(df), ncol = length(cols)))
  names(out) <- cols
  existing <- cols[cols %in% names(df)]
  out[, existing] <- df[, existing]
  out
}

h5createGroup(h5file, "number_counts")
for(i in 1:length(field_names)){
  
  phot = Rfits_read(paste0(main_dir, "/data/catalogues/", field_names[i], "_profound_multiband_cat.fits"))

  frame_info = h5read(
    h5file, name = paste0("frame_stats/", field_names[i])
  )

  df = number_counts(
    df_flux = safe_select(phot[[df_flux_name]][phot$ProFound$flag == 0, ], flux_col_names),
    limiting_mags = -2.5*log10(
      safe_select(frame_info, paste0(filters_info$filters, "_5sigma_depth_32asec"))
    ) + 8.9,
    area = safe_select(frame_info, paste0(filters_info$filters, "_effarea_amin2"))
  )
  
  path <- paste0("number_counts", "/", field_names[i])
  
  fid <- H5Fopen(h5file)
  exists <- H5Lexists(fid, path)
  H5Fclose(fid)
  
  if (exists) {
    h5delete(h5file, path)
  }
  
  h5write(as.data.frame(df), file = h5file, name = path)
  
}


## number counts from literature
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
fwrite(
  stone25,
  "~/Documents/Epochs/data/literature/number_counts/stone25.csv"
)

phot = Rfits_read(paste0(main_dir, "/data/catalogues/", "miri", "_profound_multiband_cat.fits"))
df_flux = safe_select(phot[[df_flux_name]][phot$ProFound$flag == 0, ], flux_col_names)
per_band_sources = vector(length = 31)
for(i in 1:31){
  per_band_sources[i] = sum(df_flux[,i] != -999)
}
magplot(
  filters_info$cenwave,
  per_band_sources/nrow(df_flux),
  log = "xy"
)
text(
  filters_info$cenwave,
  per_band_sources/nrow(df_flux),
  labels = filters_info$filters
)
