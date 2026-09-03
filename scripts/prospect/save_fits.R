library(Rfits)
library(ProSpect)

save_fits = function(save_stub, fits_stub){
  
  foo = readRDS(paste0(save_stub))
  
  df = data.frame(foo[c("obswave", "fluxQ50", "fluxQ16", "fluxQ84", "fluxQ05", "fluxQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 1,
    extname = "obs",
    create_file = TRUE,
    create_ext = TRUE
  )
  
  df = data.frame(foo[c("restWave", "lumQ50", "lumQ16", "lumQ84", "lumQ05", "lumQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 2,
    extname = "lum",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("StarsUnAttenWave", "StarsUnAttenQ50", "StarsUnAttenQ16", "StarsUnAttenQ84", "StarsUnAttenQ05", "StarsUnAttenQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 3,
    extname = "StarsUnAtten",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("DustWave", "DustEmitQ50", "DustEmitQ16", "DustEmitQ84", "DustEmitQ05", "DustEmitQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 4,
    extname = "Dust",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("AGNWave", "AGNQ50", "AGNQ16", "AGNQ84", "AGNQ05", "AGNQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 5,
    extname = "AGN",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
  
  df = data.frame(foo[c("ageQ50", "ageQ16", "ageQ84", "ageQ05", "ageQ95","SFHQ50", "SFHQ16", "SFHQ84", "SFHQ05", "SFHQ95","ZHQ50", "ZHQ16", "ZHQ84", "ZHQ05", "ZHQ95")])
  Rfits_write_table(
    df,
    filename = fits_stub,
    ext = 6,
    extname = "evol",
    create_file = FALSE,
    create_ext = TRUE,
    overwrite_file = FALSE
  )
}

save_fits(
  save_stub = "~/Documents/Epochs/data/literature_lrds/Cliff_lrdmbb_sed.rds",
  fits_stub = "~/Documents/Epochs/data/literature_lrds/Cliff_lrdmbb_sed.fits"
)
save_fits(
  save_stub = "~/Documents/Epochs/data/literature_lrds/Cliff_fritz_sed.rds",
  fits_stub = "~/Documents/Epochs/data/literature_lrds/Cliff_fritz_sed.fits"
)
