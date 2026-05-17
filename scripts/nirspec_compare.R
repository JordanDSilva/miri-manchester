library(ProSpect)
library(ProSpectUtils)
library(celestial)
library(ParmOff)
library(Highlander)
library(nloptr)
library(foreach)
library(matrixStats)
library(Rfits)
library(doParallel)
library(magicaxis)
library(data.table)
library(scales)
library(RColorBrewer)

set.seed(666)

main_dir = "~/Documents/Epochs/"
# main_dir = "/scratch/"

Niters = 3000
ncores = 40

frame_info = data.frame(fread(paste0(main_dir, "data/depths/smiles_profound_MIRI_frame_stats.csv")))
filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

phot_cat = Rfits_read(paste0(main_dir, "/data/catalogues/smiles_profound_multiband_cat.fits"))
combined_cat = data.frame(phot_cat$Phot32)

df_flux = combined_cat[, paste0(filters_info$filters, "_flux_aper32")]
df_flux_err = combined_cat[, paste0(filters_info$filters, "_flux_aper32_err_5pc")]

miri_nirspec = Rfits_read_table(paste0(main_dir, "data/hlsp_smiles_jwst_nirspec_goodss_multi_v1.0_catalog.fits"))

match_nirspec_idx = coordmatch(
  coordref = miri_nirspec[, c("RA", "DEC")],
  coordcompare = phot_cat$ProFound[, c("RA", "DEC")], 
  rad = 0.3
)
match_nirspec_idx$Nmatch

nirspec_z = miri_nirspec$Z_SPEC[match_nirspec_idx$bestmatch$refID]

photoZ_bc03_snorm_fritz_astro = fread("~/mnt/morgan/canopus/prospect/photoZ_bc03_snorm_fritz_astro_cat.csv")

calc_nmad_z = function(zspec, zphot){
  deltaz = zphot - zspec
  inner = (deltaz - median(deltaz)) / (1+zspec)
  sigmaNMAD = 1.48 * median(abs(inner))
  return(sigmaNMAD)
}

calc_nmad_z(
  nirspec_z, photoZ_bc03_snorm_fritz_astro$zQ50
)

magplot(
  nirspec_z,
  photoZ_bc03_snorm_fritz_astro$zQ50,
  xlim = c(0, 10),
  ylim = c(0, 25)
)
magerr(
  nirspec_z, 
  photoZ_bc03_snorm_fritz_astro$zQ50,
  ylo = photoZ_bc03_snorm_fritz_astro$z_lo,
  yhi = photoZ_bc03_snorm_fritz_astro$z_hi,
)
abline(0,1)

magplot(
  photoZ_bc03_snorm_fritz_astro$zQ50,
  log10(photoZ_bc03_snorm_fritz_astro$StellarMassQ50)
)
magplot(
  log10(photoZ_bc03_snorm_fritz_astro$StellarMassQ50),
  log10(photoZ_bc03_snorm_fritz_astro$SFR10Q50), 
  z = log10(photoZ_bc03_snorm_fritz_astro$AGNlumQ50), 
  xlim = c(5, 13),
  ylim = c(-5, 5)
)
