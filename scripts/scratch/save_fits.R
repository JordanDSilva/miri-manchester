library(dplyr)
library(data.table)
library(Rfits)

for(i in c(1:109, 111:129, 131:144)){
  foo = readRDS(paste0("~/mnt/csf3/forAlbert/data/prospect_seds/photoZ_bc03_snorm_fritz_lbol_prior/galID_", i, ".rds"))
  
  save_fits_file = paste0("~/mnt/csf3/forAlbert/data/prospect_seds/photoZ_bc03_snorm_fritz_lbol_prior_fits/galID_", i, ".fits")
  
  df = data.frame(foo[c("obswave", "fluxQ50", "fluxQ16", "fluxQ84", "fluxQ05", "fluxQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 1,
    extname = "obs",
    create_file = TRUE,
    create_ext = TRUE
  )
  
  df = data.frame(foo[c("restWave", "lumQ50", "lumQ16", "lumQ84", "lumQ05", "lumQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 2,
    extname = "lum",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("StarsUnAttenWave", "StarsUnAttenQ50", "StarsUnAttenQ16", "StarsUnAttenQ84", "StarsUnAttenQ05", "StarsUnAttenQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 3,
    extname = "StarsUnAtten",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("DustWave", "DustEmitQ50", "DustEmitQ16", "DustEmitQ84", "DustEmitQ05", "DustEmitQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 4,
    extname = "Dust",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("AGNWave", "AGNQ50", "AGNQ16", "AGNQ84", "AGNQ05", "AGNQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 5,
    extname = "AGN",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("ageQ50", "ageQ16", "ageQ84", "ageQ05", "ageQ95","SFHQ50", "SFHQ16", "SFHQ84", "SFHQ05", "SFHQ95","ZHQ50", "ZHQ16", "ZHQ84", "ZHQ05", "ZHQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 6,
    extname = "evol",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
}

for(i in c(1:109, 111:129, 131:144)){
  foo = readRDS(paste0("~/mnt/csf3/forAlbert/data/prospect_seds/photoZ_bc03_snorm/galID_", i, ".rds"))
  
  save_fits_file = paste0("~/mnt/csf3/forAlbert/data/prospect_seds/photoZ_bc03_snorm_fits/galID_", i, ".fits")
  
  df = data.frame(foo[c("obswave", "fluxQ50", "fluxQ16", "fluxQ84", "fluxQ05", "fluxQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 1,
    extname = "obs",
    create_file = TRUE,
    create_ext = TRUE
  )
  
  df = data.frame(foo[c("restWave", "lumQ50", "lumQ16", "lumQ84", "lumQ05", "lumQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 2,
    extname = "lum",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("StarsUnAttenWave", "StarsUnAttenQ50", "StarsUnAttenQ16", "StarsUnAttenQ84", "StarsUnAttenQ05", "StarsUnAttenQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 3,
    extname = "StarsUnAtten",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("DustWave", "DustEmitQ50", "DustEmitQ16", "DustEmitQ84", "DustEmitQ05", "DustEmitQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 4,
    extname = "Dust",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  # df = data.frame(foo[c("AGNWave", "AGNQ50", "AGNQ16", "AGNQ84", "AGNQ05", "AGNQ95")])
  # Rfits_write_table(
  #   df,
  #   filename = save_fits_file,
  #   ext = 5,
  #   extname = "AGN",
  #   create_file = FALSE,
  #   create_ext = TRUE,
  #   overwrite_file = FALSE
  # )
  
  df = data.frame(foo[c("ageQ50", "ageQ16", "ageQ84", "ageQ05", "ageQ95","SFHQ50", "SFHQ16", "SFHQ84", "SFHQ05", "SFHQ95","ZHQ50", "ZHQ16", "ZHQ84", "ZHQ05", "ZHQ95")])
  Rfits_write_table(
    df,
    filename = save_fits_file,
    ext = 5,
    extname = "evol",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
}
# data_names = c(
#   "restWave", "lumQ50", "lumQ16", "lumQ84", "lumQ05", "lumQ95",
#   
#   "StarsUnAttenWave", "StarsUnAttenQ50", "StarsUnAttenQ16", "StarsUnAttenQ84", "StarsUnAttenQ05", "StarsUnAttenQ95",
#   
#   "StarsAttenWave", "StarsAttenQ50", "StarsAttenQ16", "StarsAttenQ84", "StarsAttenQ05", "StarsAttenQ95",
#   
#   "DustWave", "DustEmitQ50", "DustEmitQ16", "DustEmitQ84", "DustEmitQ05", "DustEmitQ95",
#   
#   "AGNWave", "AGNQ50", "AGNQ16", "AGNQ84", "AGNQ05", "AGNQ95",
#   
#   "obswave", "fluxQ50", "fluxQ16", "fluxQ84", "fluxQ05", "fluxQ95",
#   
#   "ageQ50", "ageQ16", "ageQ84", "ageQ05", "ageQ95",
#   
#   "SFHQ50", "SFHQ16", "SFHQ84", "SFHQ05", "SFHQ95",
#   
#   "ZHQ50", "ZHQ16", "ZHQ84", "ZHQ05", "ZHQ95"
# )


