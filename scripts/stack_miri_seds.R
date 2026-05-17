library(scales)
library(RColorBrewer)
library(ProSpect)
library(Highlander)
library(celestial)
library(data.table)
library(ParmOff)
library(foreach)
library(matrixStats)
library(magicaxis)
library(Rfits)

combined_cat = data.frame(fread(
  paste0("~/Documents/Epochs/data/catalogues/smiles_profound_MIRI_matched_EPOCHSv2_photometry_cat.csv")
))

filter_transforms = readRDS("~/Documents/Epochs/data/throughputs/filt_funcs.rds")
filters_info = fread("~/Documents/Epochs/data/throughputs/wavelengths.txt", col.names = c("wavelengths", "filters"))

colors = colorRampPalette(
  colors =rev(brewer.pal(10, name = "Spectral"))
)(length(filters_info$filters))

zcut_lo = 0.5
zcut_hi = 1.5

z_idx = combined_cat$zbest_sfhz_zfree <= zcut_hi & combined_cat$zbest_sfhz_zfree >= zcut_lo & combined_cat$zbest_sfhz_zfree != -999
med_z = median(combined_cat$zbest_sfhz_zfree[z_idx])

df_flux = combined_cat[z_idx, paste0(filters_info$filters, "_flux_aper32")]
df_flux_err = combined_cat[z_idx, paste0(filters_info$filters, "_flux_aper32_err_5pc")]

average_flux = unlist(colSums(as.matrix(df_flux), na.rm = TRUE)) / sapply(df_flux, function(x){sum(x>0)})
average_flux_err = sqrt( unlist(colSums(as.matrix(df_flux_err)^2)) ) / sapply(df_flux, function(x){sum(x>0)})

magplot(
  NA,
  log = "xy",
  ylim = c(1e-9, max(average_flux, na.rm = TRUE)*1.5),
  xlim = c(0.4, 30)
)
points(
  filters_info$wavelengths/1e4,
  average_flux
)
magerr(
  filters_info$wavelengths/1e4,
  average_flux,
  ylo = (average_flux_err)
)
for(i in 1:length(filters_info$filters)){
  lines(
    10^seq(1,8,0.001)/1e4,
    filter_transforms[[i]](10^seq(1,8,0.001)) * 1e-7,
    col = alpha(colors[i], 0.5)
  )
}