library(data.table)
library(Rfits)
library(rhdf5)
library(magicaxis)
library(celestial)
library(dplyr)

load_eazy = function(fname, fname_append){
  # h5ls(
  #   "~/Documents/JWST-MIRI-AGN/data/eazy/smiles/smiles_zphot.out.h5"
  # )
  eazy_zPDFs = data.frame(fread(
    paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF", fname_append, ".csv")
  ))
  eazy_ntemplates = h5read(
    paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_zphot", fname_append, ".out.h5"),
    'fit/coeffs_best'
  )
  eazy_ntemplates = colSums(
    eazy_ntemplates > 0
  )
  eazy = Rfits_read_table(
    paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_zphot", fname_append, ".out.zout.fits")
  )
  eazy$ntemplates = eazy_ntemplates
  
  if(!file.exists(paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv"))){
    eazy_stats = NULL
  }else{
    eazy_stats = data.frame(fread(
      paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv")
    ))
    eazy_stats$z_ml = eazy$z_ml
  }
  
  if(!file.exists(paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF_normal", fname_append, ".csv"))){
    eazy_zPDFs_normals = NULL
  }else{
    eazy_zPDFs_normals = data.frame(fread(
      paste0("~/Documents/JWST-MIRI-AGN/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF_normal", fname_append, ".csv")
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

agn_blue_sfhz = load_eazy("miri", "_agn_blue_sfhz_col")

filter_info = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/JWST-MIRI-AGN/data/catalogues/miri_profound_multiband_cat.fits", data.table = TRUE)

prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"))
# prospect_cat = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))

df = data.frame(
  "zspec" = agn_blue_sfhz$eazy$z_spec,
  "zeazy" = agn_blue_sfhz$eazy$z_ml,
  "zQ50" = prospect_cat$zQ50,
  "zQ16" = prospect_cat$zQ16,
  "zQ84" = prospect_cat$zQ84,
  "z_lo" = prospect_cat$z_lo,
  "z_hi" = prospect_cat$z_hi,
  "eazyChi2" = agn_blue_sfhz$eazy_stats$zbestChi2,
  "eazyHessianSD" = agn_blue_sfhz$eazy_stats$zbestHessianSD,
  "eazyPz80" = agn_blue_sfhz$eazy_stats$zbestPzIntegral
)

# df[
#   which(df$zspec < 4 & df$zeazy > 4 & df$zspec != -999),
# ]
# 
# prospect_cat[which(df$zspec < 4 & df$zeazy > 6 & df$zspec != -999), c("UV1500Q50", "StellarMassQ50", "FitType")]

h5file = "~/Documents/JWST-MIRI-AGN/data/analysis_data.h5"
h5createFile(h5file)
h5createGroup(
  file = h5file,
  "photoz"
)
# h5delete(h5file, "photoz/photoz")
h5write(
  obj = df,
  file = h5file,
  name = "photoz/photoz"
)
# h5delete(h5file, "photoz/eazyPz")
h5write(
  obj = as.matrix(agn_blue_sfhz$eazy_zPDFs[, !grepl("V1|zgrid", names(agn_blue_sfhz$eazy_zPDFs))]),
  file = h5file,
  name = "photoz/eazyPz"
)
# h5delete(h5file, "photoz/eazyPzNorm")
h5write(
  obj = as.matrix(agn_blue_sfhz$eazy_zPDFs_normals[, !grepl("V1|zgrid", names(agn_blue_sfhz$eazy_zPDFs_normals))]),
  file = h5file,
  name = "photoz/eazyPzNorm"
)
h5write(
  obj = agn_blue_sfhz$eazy_zPDFs$zgrid,
  file = h5file,
  name = "photoz/zgrid"
)

# idx = which(df$zeazy > 5 & df$zspec < 4 & df$zspec != -999)
# 
# df$zQ50[idx]
# df$zeazy[idx]
# 
# phot_cat$ProFound$flag[idx]
# 
# magplot(agn_blue_sfhz$eazy$dec, prospect_cat$Dec)
