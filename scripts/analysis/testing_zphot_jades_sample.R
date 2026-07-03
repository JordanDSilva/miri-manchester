library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)

compute_sigmaNMAD <- function(zspec, zphot, zphot_lim = 0) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz   <- zphot[idx] - zspec[idx]
  
  sigmaNMAD <- 1.48 * median( abs( (deltaz - median(deltaz)) / (1 + zspec[idx]) ) )
  return(sigmaNMAD)
}
compute_outliers <- function(zspec, zphot, thresh = 0.15, zphot_lim = -2) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz     <- abs(zphot[idx] - zspec[idx])
  Noutlier   <- sum(deltaz/(1 + zspec[idx]) > thresh)
  outlier_frac <- Noutlier / length(zspec[idx])
  
  cat("Noutlier:", Noutlier, "\n")
  cat("N total: ", length(zspec[idx]), "\n")
  
  return(outlier_frac)
}
zspec_zphot_plot = function(zphot, zspec, zphot_lim = 0){
  magplot(
    zphot,
    zspec, 
    xlim = c(0, 15),
    ylim = c(0, 15),
    xlab = "zphot",
    ylab = "zspec"
  )
  legend(
    x = "topright",
    legend = c(
      paste0("sigmaNMAD = ", round(compute_sigmaNMAD(zspec, zphot, zphot_lim = zphot_lim), 3)),
      paste0("Outliers = ", round(compute_outliers(zspec, zphot, zphot_lim = zphot_lim), 3))
    )
  )
}

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)
h5file = "~/Documents/Epochs/data/analysis_data.h5"

# Eazy stuff
load_eazy = function(fname, fname_append){
  # h5ls(
  #   "~/Documents/Epochs/data/eazy/smiles/smiles_zphot.out.h5"
  # )
  eazy_zPDFs = data.frame(fread(
    paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF", fname_append, ".csv")
  ))
  eazy_ntemplates = h5read(
    paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_zphot", fname_append, ".out.h5"),
    'fit/coeffs_best'
  )
  eazy_ntemplates = colSums(
    eazy_ntemplates > 0
  )
  eazy = Rfits_read_table(
    paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_zphot", fname_append, ".out.zout.fits")
  )
  eazy$ntemplates = eazy_ntemplates
  
  if(!file.exists(paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv"))){
    eazy_stats = NULL
  }else{
    eazy_stats = data.frame(fread(
      paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv")
    ))
    eazy_stats$z_ml = eazy$z_ml
  }

  if(!file.exists(paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF_normal", fname_append, ".csv"))){
    eazy_zPDFs_normals = NULL
  }else{
    eazy_zPDFs_normals = data.frame(fread(
      paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF_normal", fname_append, ".csv")
    ))
  }
  return(list(eazy = eazy, eazy_stats = eazy_stats, eazy_zPDFs = eazy_zPDFs, eazy_zPDFs_normals = eazy_zPDFs_normals))
}

compute_peak_integrals = function(kk, eazy_zPDFs, eazy){
  zgrid = eazy_zPDFs$zgrid
  zPDF = eazy_zPDFs[[paste0('gal',kk)]]
  
  zbest = eazy$z_ml[kk]
  
  idx = zgrid >= 0.9*zbest & zgrid <= 1.1*zbest
  
  trapz(
    zgrid[idx],
    zPDF[idx]
  )
}

agn_blue_sfhz = load_eazy("smiles", "_agn_blue_sfhz")
agn_blue_sfhz_no_medium = load_eazy("smiles", "_agn_blue_sfhz_no_medium")

sfhz_carnall = load_eazy("smiles", "_sfhz_carnall")
larson_default = load_eazy("smiles", "_larson")

jades_prism_spec = data.frame(Rfits_read_table(
  "~/Documents/Epochs/data/spec/hlsp_jades_jwst_nirspec_goods-s_prism-line-fluxes_v1.1_catalog.fits"
))
jades_robust = jades_prism_spec[jades_prism_spec$z_Spec_flag == "A", ]
jades_prism_spec_match = coordmatch(
  coordref = jades_robust[, c("RA_NIRCam", "Dec_NIRCam")],
  coordcompare = phot_cat$ProFound[, c("RAmax", "Decmax")],
  rad = 0.3
)
jades_match = jades_robust[jades_prism_spec_match$bestmatch$refID, ]

zspec_zphot_plot(
  zphot = larson_default$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = sfhz_carnall$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = agn_blue_sfhz_no_medium$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  jades_match$z_Spec,
  zphot_lim = 0
)

idx_select = agn_blue_sfhz$eazy_stats$zbestPzIntegral[jades_prism_spec_match$bestmatch$compareID] >= 0.60
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID][idx_select],
  jades_match$z_Spec[idx_select],
  zphot_lim = 0
)

juodzbalis26 = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_agn_data_raw.csv"))
juodzbalis26_sed = data.frame(fread("~/Documents/Epochs/data/literature/AGN-Catalogues/juodzbalis26_sed_fits_raw.csv"))
juodzbalis26_combine = cbind(juodzbalis26, juodzbalis26_sed)
photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
match_juodzbalis26 = coordmatch(
  coordref = juodzbalis26[,c("RA", "Dec")],
  coordcompare = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth[, c("RA", "Dec")],
  rad = 0.3
)

h5createFile(h5file)
h5createGroup(
  file = h5file,
  "eazy"
)
h5createGroup(
  file = h5file,
  "eazy/eazy_stats"
)
h5createGroup(
  file = h5file,
  "eazy/eazy_pdfs"
)
h5createGroup(
  file = h5file,
  "eazy/eazy_norm_pdfs"
)

# h5delete(h5file, "eazy_stats")  
# h5delete(h5file, "eazy_pdfs")  
# h5delete(h5file, "eazy_norm_pdfs")  

h5write(
  obj = agn_blue_sfhz$eazy_zPDFs_normals$zgrid,
  file = h5file,
  name = "eazy/eazy_norm_pdfs/zgrid"
)
h5write(
  obj = agn_blue_sfhz$eazy_zPDFs$zgrid,
  file = h5file,
  name = "eazy/eazy_pdfs/zgrid"
)

h5write(
  obj = sfhz_carnall$eazy_stats,
  file = h5file,
  name = "eazy/eazy_stats/sfhz_carnall"
)
h5write(
  obj = as.matrix(sfhz_carnall$eazy_zPDFs[grepl("gal", names(sfhz_carnall$eazy_zPDFs))]),
  file = h5file,
  name = "eazy/eazy_pdfs/sfhz_carnall"
)
h5write(
  obj = as.matrix(sfhz_carnall$eazy_zPDFs_normals[grepl("gal", names(sfhz_carnall$eazy_zPDFs_normals))]),
  file = h5file,
  name = "eazy/eazy_norm_pdfs/sfhz_carnall"
)

h5write(
  obj = larson_default$eazy_stats,
  file = h5file,
  name = "eazy/eazy_stats/larson_fsps"
)
h5write(
  obj = as.matrix(larson_default$eazy_zPDFs[grepl("gal", names(larson_default$eazy_zPDFs))]),
  file = h5file,
  name = "eazy/eazy_pdfs/larson_fsps"
)
h5write(
  obj = as.matrix(larson_default$eazy_zPDFs_normals[grepl("gal", names(larson_default$eazy_zPDFs_normals))]),
  file = h5file,
  name = "eazy/eazy_norm_pdfs/larson_fsps"
)

h5write(
  obj = agn_blue_sfhz$eazy_stats,
  file = h5file,
  name = "eazy/eazy_stats/agn_blue_sfhz"
)
h5write(
  obj = as.matrix(agn_blue_sfhz$eazy_zPDFs[grepl("gal", names(agn_blue_sfhz$eazy_zPDFs))]),
  file = h5file,
  name = "eazy/eazy_pdfs/agn_blue_sfhz"
)
h5write(
  obj = as.matrix(agn_blue_sfhz$eazy_zPDFs_normals[grepl("gal", names(agn_blue_sfhz$eazy_zPDFs_normals))]),
  file = h5file,
  name = "eazy/eazy_norm_pdfs/agn_blue_sfhz"
)

h5delete(h5file, "specz")
h5createGroup(
  h5file, 
  group = "specz"
)
h5createGroup(
  h5file, 
  group = "specz/jades_match"
)
h5write(
  obj = jades_match,
  file = h5file,
  name = "specz/jades_match/specz"
)
h5write(
  obj = data.frame(agn_blue_sfhz$eazy_stats[jades_prism_spec_match$bestmatch$compareID, ], "ID" = jades_prism_spec_match$bestmatch$compareID),
  file = h5file,
  name = "specz/jades_match/agn_blue_sfhz_stats"
)

h5createGroup(
  h5file, 
  group = "specz/juodzbaliz_match"
)
h5write(
  obj = data.frame(juodzbalis26_combine[match_juodzbalis26$bestmatch$refID, ]),
  file = h5file,
  name = "specz/juodzbaliz_match/specz"
)
h5write(
  obj = data.frame(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth[match_juodzbalis26$bestmatch$compareID, ], "ID" = match_juodzbalis26$bestmatch$compareID),
  file = h5file,
  name = "specz/juodzbaliz_match/prospect_with_agn"
)


kk = 56
magplot(
  larson_default$eazy_zPDFs$zgrid,
  larson_default$eazy_zPDFs[[paste0("gal", jades_prism_spec_match$bestmatch$compareID[kk])]], 
  type = "l", 
  xlim = c(-1, 8)
)
lines(
  agn_blue_sfhz$eazy_zPDFs$zgrid,
  agn_blue_sfhz$eazy_zPDFs[[paste0("gal", jades_prism_spec_match$bestmatch$compareID[kk])]], 
  type = "l", 
  col = "red"
)
lines(
  sfhz_carnall$eazy_zPDFs$zgrid,
  sfhz_carnall$eazy_zPDFs[[paste0("gal", jades_prism_spec_match$bestmatch$compareID[kk])]], 
  type = "l", 
  col = "blue"
)
abline(
  v = jades_match$z_Spec[kk]
)

photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth_spec_astro_cat.csv"
))
photoZ_bc03_snorm_fritz_no_tau_birth = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_no_tau_birth_spec_astro_cat.csv"
))
photoZ_bc03_snorm_no_tau_birth = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_no_tau_birth_spec_astro_cat.csv"
))

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth$zQ50,
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_fritz_no_tau_birth$zQ50,
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_no_tau_birth$zQ50,
  jades_match$z_Spec,
  zphot_lim = 0
)

photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth_z1 = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth_z1_spec_astro_cat.csv"))
photoZ_bc03_snorm_fritz_no_tau_birth_z1 = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_no_tau_birth_z1_spec_astro_cat.csv"))

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth_z1$zQ50,
  jades_match$z_Spec,
  zphot_lim = 0
)
zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_fritz_no_tau_birth_z1$zQ50,
  jades_match$z_Spec,
  zphot_lim = 0
)

photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth_astro_cat.csv"))
photoZ_bc03_snorm_norm_eazy_no_tau_birth = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_norm_eazy_no_tau_birth_astro_cat.csv"))
photoZ_bc03_snorm_fritz_robust = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))

magplot(
  photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50,
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$DustMassTotalQ50), 
  xlim = c(1, 9),
  ylim = c(1,12)
)
magplot(
  photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50,
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50), 
  xlim = c(1, 9),
  ylim = c(35,50)
)

pdf("~/Documents/Epochs/plots/zphot_zspec_different_prospect_runs.pdf", width = 8, height = 5)
zspec_zphot_plot(
  zphot = eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  # zphot = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50[jades_prism_spec_match$bestmatch$compareID],
  # dja_spec_match$z, 
  jades_match$z_Spec,
  zphot_lim = 0
)
mtext(
  text = "EAZY zbest"
)

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_norm_eazy_no_tau_birth$zQ50,
  jades_match$z_Spec, 
  zphot_lim = 0
)
mtext(
  text = "ProSpect normal prior on eazy zbest"
)

# zspec_zphot_plot(
#   zphot = photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50,
#   jades_match$z_Spec, 
#   zphot_lim = 0
# )

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_no_tau_birth$zQ50,
  jades_match$z_Spec
)
mtext(
  text = "ProSpect full eazy prior"
)

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_no_eazy_no_tau_birth$zQ50,
  jades_match$z_Spec
)
mtext(
  text = "ProSpect no eazy prior"
)
dev.off()
# zspec_zphot_plot(
#   zphot = photoZ_bc03_snorm_fritz_no_eazy_no_tau_birth$zQ50,
#   jades_match$z_Spec
# )

zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_no_tau_birth$zQ50,
  jades_match$z_Spec
)
zspec_zphot_plot(
  zphot = photoZ_bc03_snorm_fritz_no_tau_birth$zQ50,
  jades_match$z_Spec
)

zspec_zphot_plot(
    zphot = eazy$z_ml[jades_prism_spec_match$bestmatch$compareID],
  jades_match$z_Spec
)



pdf("~/Desktop/compare_with_Juodzbalis26.pdf", width = 8, height = 6)
magplot(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  (photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50[match_juodzbalis26$bestmatch$compareID]), 
  xlab = "Spec z",
  ylab = "ProSpect Photoz", 
  xlim = c(3,6),
  ylim = c(3.6), 
  grid = FALSE,
  col = viridis::inferno(11), 
  pch = 16
)
magerr(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  unname(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$zQ50[match_juodzbalis26$bestmatch$compareID]),
  ylo = unname(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$z_lo[match_juodzbalis26$bestmatch$compareID]),
  yhi = unname(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$z_hi[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)

magplot(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlim = c(43, 46),
  ylim = c(35, 46),
  xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
  ylab = "ProSpect AGN Lbol [erg/s]", 
  grid = FALSE,
  col = viridis::rocket(11), 
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26$log_Lbol_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26$log_Lbol_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)

magplot(
  juodzbalis26_sed$log_Mstar_BEAGLE[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  xlim = c(6, 11.0),
  ylim = c(6, 11.0),
  xlab = "Juodzbalis+26 Stellar Mass [Msun]",
  ylab = "ProSpect Stellar Mass [Msun]",
  grid = FALSE,
  col = viridis::rocket(11), 
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26_sed$log_Mstar_BEAGLE[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26_sed$log_Mstar_BEAGLE_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26_sed$log_Mstar_BEAGLE_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)
dev.off()

magplot(
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlim = c(8, 10.0),
  ylim = c(8, 10.0),
  xlab = "Juodzbalis+26 Stellar Mass [Msun]",
  ylab = "ProSpect Stellar Mass [Msun]"
)
magerr(
  juodzbalis26_sed$log_Mstar_BEAGLE[match_idx$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_idx$bestmatch$compareID]),
  xlo = juodzbalis26_sed$log_Mstar_BEAGLE_lo[match_idx$bestmatch$refID],
  xhi = juodzbalis26_sed$log_Mstar_BEAGLE_hi[match_idx$bestmatch$refID],
  ylo = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_idx$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ16[match_idx$bestmatch$compareID]),
  yhi = log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ84[match_idx$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$StellarMassQ50[match_idx$bestmatch$compareID]),
)
abline(0,1, lty = 2)



# Make sure it is all consistent
# zbest_integrals = sapply(jades_prism_spec_match$bestmatch$compareID, function(x){compute_peak_integrals(x, agn_blue_sfhz$eazy_zPDFs, agn_blue_sfhz$eazy)})
# idx_select = zbest_integrals >= 0.6
# zspec_zphot_plot(
#   zphot = agn_blue_sfhz$eazy$z_ml[jades_prism_spec_match$bestmatch$compareID][idx_select],
#   jades_match$z_Spec[idx_select],
#   zphot_lim = 0
# )

magplot(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  (photoZ_bc03_snorm_fritz_robust$zQ50[match_juodzbalis26$bestmatch$compareID]), 
  xlab = "Spec z",
  ylab = "ProSpect Photoz", 
  xlim = c(3,6),
  ylim = c(3.6), 
  grid = FALSE,
  col = viridis::inferno(11), 
  pch = 16
)
magerr(
  juodzbalis26$z[match_juodzbalis26$bestmatch$refID],
  unname(photoZ_bc03_snorm_fritz_robust$zQ50[match_juodzbalis26$bestmatch$compareID]),
  ylo = unname(photoZ_bc03_snorm_fritz_robust$z_lo[match_juodzbalis26$bestmatch$compareID]),
  yhi = unname(photoZ_bc03_snorm_fritz_robust$z_hi[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)

magplot(
  juodzbalis26_sed$log_Mstar_BEAGLE[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  z = photoZ_bc03_snorm_fritz_robust$FitType == "StellarAGN",
  xlim = c(6, 12.0),
  ylim = c(6, 12.0),
  xlab = "Juodzbalis+26 Stellar Mass [Msun]",
  ylab = "ProSpect Stellar Mass [Msun]",
  grid = FALSE,
  # col = viridis::rocket(11), 
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26_sed$log_Mstar_BEAGLE[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26_sed$log_Mstar_BEAGLE_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26_sed$log_Mstar_BEAGLE_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_robust$StellarMassQ50[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1)

magplot(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_robust$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  z = photoZ_bc03_snorm_fritz_robust$FitType == "StellarAGN",
  xlim = c(43, 46),
  ylim = c(35, 46),
  xlab = "Juodzbalis+26 Ha broadline AGN Lbol [erg/s]",
  ylab = "ProSpect AGN Lbol [erg/s]", 
  grid = FALSE,
  # col = viridis::rocket(11), 
  pch = 16, 
  cex = 1.5
)
magerr(
  juodzbalis26$log_Lbol[match_juodzbalis26$bestmatch$refID],
  log10(photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
  xlo = juodzbalis26$log_Lbol_lo[match_juodzbalis26$bestmatch$refID],
  xhi = juodzbalis26$log_Lbol_hi[match_juodzbalis26$bestmatch$refID],
  ylo = log10(photoZ_bc03_snorm_fritz_robust$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_robust$AGNlumQ16[match_juodzbalis26$bestmatch$compareID]),
  yhi = log10(photoZ_bc03_snorm_fritz_robust$AGNlumQ84[match_juodzbalis26$bestmatch$compareID]) -
    log10(photoZ_bc03_snorm_fritz_robust$AGNlumQ50[match_juodzbalis26$bestmatch$compareID]),
)
abline(0,1, lty = 2)
