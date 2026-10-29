library(data.table)
library(Rfits)
library(Rwcs)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)
library(bit64)

filter_info = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)
h5file = "~/Documents/JWST-MIRI-AGN/data/analysis_data.h5"

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
AGN_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
Stellar_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

sample_idx = prospect_cat$zQ50 >= 1.07467287 & prospect_cat$zQ50 < 9.98276084 & !is.na(prospect_cat$zQ50)

prospect_cat = prospect_cat[sample_idx, ]

#covelo paz Halpha http://arxiv.org/abs/2409.17241
covelo_paz_fname = "~/Downloads/J_A+A_694_A178/catalog.dat"
lines <- readLines(covelo_paz_fname)
lines <- gsub(",\\s+", ",", lines)   # "12143, 12144, 12145" -> "12143,12144,12145"

covelo_paz25 <- read.table(
  text = lines,
  header = FALSE,
  stringsAsFactors = FALSE,
  na.strings = "-",                  # fresco-only rows have no member IDs
  col.names = c("field", "id", "members", "ra", "dec", "grade", "z",
                "flux", "flux_err", "muv", "muv_err", "ew", "ew_err")
)

covelo_paz25_match_idx = coordmatch(
  coordref = covelo_paz25[, c("ra", "dec")],
  coordcompare = prospect_cat[, c("RA", "Dec")],
  rad = 0.3
)

covelo_paz25_match = covelo_paz25[covelo_paz25_match_idx$bestmatch$refID, ]
prospect_covelo_paz25_match = prospect_cat[covelo_paz25_match_idx$bestmatch$compareID, ]

A_Ha = 0.47

covelo_paz25_match$LHaRaw = covelo_paz25_match$flux * 
  (cosdistLumDist(z = covelo_paz25_match$z, ref = "737")*3.086e+24)^2 * 4*pi
covelo_paz25_match$LHaCorr = covelo_paz25_match$LHaRaw * 10^(0.4 * A_Ha)
covelo_paz25_match$SFRHa = covelo_paz25_match$LHaCorr * 10^-41.36

covelo_paz25_match$LHaRawErr = covelo_paz25_match$flux_err * 
  (cosdistLumDist(z = covelo_paz25_match$z, ref = "737")*3.086e+24)^2 * 4*pi
covelo_paz25_match$LHaCorrErr = covelo_paz25_match$LHaRawErr * 10^(0.4 * A_Ha)
covelo_paz25_match$SFRHaErr = covelo_paz25_match$LHaCorrErr * 10^-41.36

median(
  log10(prospect_covelo_paz25_match$SFR10Q84) - log10(prospect_covelo_paz25_match$SFR10Q50)
)
median(
  log10(prospect_covelo_paz25_match$SFR10Q50) - log10(prospect_covelo_paz25_match$SFR10Q16)
)

magplot(
  covelo_paz25_match$SFRHa,
  prospect_covelo_paz25_match$SFR10Q50,
  log = "xy",
  col = ifelse(prospect_covelo_paz25_match$FitType == "Stellar", "blue", "red"),
  xlim = c(1e0, 2e2),
  ylim = c(1e0, 2e2)
)
magerr(
  covelo_paz25_match$SFRHa,
  prospect_covelo_paz25_match$SFR10Q50,
  ylo = prospect_covelo_paz25_match$SFR10_lo,
  yhi = prospect_covelo_paz25_match$SFR10_hi,
  xlo = covelo_paz25_match$SFRHaErr
)
abline(0,1)
maghist(
  log10(covelo_paz25_match$SFRHa/prospect_covelo_paz25_match$SFR10Q50)
)

magplot(
  prospect_covelo_paz25_match$UV1500Q50,
  covelo_paz25_match$muv,
  # log = "xy",
  col = ifelse(prospect_covelo_paz25_match$FitType == "Stellar", "blue", "red"),
  # xlim = c(1e-3, 2e2),
  # ylim = c(1e-3, 2e2)
)
abline(0,1)

