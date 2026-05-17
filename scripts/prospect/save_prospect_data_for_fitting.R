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

## Make a .rda with all of the data I will need for fitting

main_dir = "~/Documents/Epochs/"
fname = "smiles"

smiles_miri_nirspec = Rfits_read_table(paste0(main_dir, "data/hlsp_smiles_jwst_nirspec_goodss_multi_v1.0_catalog.fits"))

todaysdate = Sys.Date()

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))

phot_cat = Rfits_read(
  paste0(main_dir, "/data/catalogues/", fname, "_profound_multiband_cat.fits")
)
## load in EAZY zPDFs and make approxfuns
## make sure to put zPDF prior back to prior func
eazy_zPDFs = data.frame(fread(
  paste0("~/Documents/Epochs/eazy/", fname, "/", fname, "_eazy_zPDF.csv")
))

pfunc=function(parm){
  
  # 100*erf(parm['mperiod']+2) - 100, ## Old Thorne+20 prior
  # -20*erf(parm['tau_screen']-2), ## Old Thorne+20 prior
  prior = sum(
    c(
      dnorm(parm['tau_birth'], mean=-0.2, sd=0.5, log = TRUE),
      100*erf(-5 * (parm['tau_screen'] - parm['tau_birth'])) - 100, ## Column density of dust in the ISM should be less than that of the birth clouds
      dnorm(parm['alpha_SF_birth'], mean=2, sd=1, log=TRUE),
      dnorm(parm['alpha_SF_screen'], mean=2, sd=1, log=TRUE),
      100*erf(-5 * (parm['mpeak'] - cosdistUniAgeAtz(10^parm['z'], ref = '737'))) - 100, ## disfavour mpeak above the age of the Uni at that epoch
      0
    ),
    na.rm = TRUE
  )
  return(prior)
}

pfunc_driver26_alpha=function(parm){
  prior = sum(
    c(
      dnorm(parm['tau_birth'],mean=-0.2, sd=0.5, log=TRUE),
      100*erf(-5 * (parm['tau_screen'] - parm['tau_birth'])) - 100, ## Column density of dust in the ISM should be less than that of the birth clouds
      100*erf(-5 * (parm['mpeak'] - cosdistUniAgeAtz(10^parm['z'], ref = '737'))) - 100, ## disfavour mpeak above the age of the Uni at that epoch
      dnorm(x = parm["alpha_SF_birth"], mean=pmax(2.2-10^parm['z']*0.28, 0.0625), sd = 1, log = TRUE),
      dnorm(x = parm["alpha_SF_screen"], mean=pmax(2.2-10^parm['z']*0.28, 0.0625), sd = 1, log = TRUE),
      0
    ),
    na.rm = TRUE
  )
  return(prior)
}

constraintfunc = function(parm, Data){
  uniage = cosdistUniAgeAtz(10^parm[grep('z', Data$parm.names)], ref = "737")
  if(parm[grep('mpeak', Data$parm.names)] > uniage){
    parm[grep('mpeak', Data$parm.names)] = uniage ## peak of star formation cannot be before the formation epoch of the galaxy (so older than the age of the Universe at the observed redshift)
  }
  
  if(parm[grep('tau_screen', Data$parm.names)] > parm[grep('tau_birth', Data$parm.names)]){ ## V-band optical depth in the ISM must be less than the birth clouds
    parm[grep('tau_screen', Data$parm.names)] = parm[grep('tau_birth', Data$parm.names)]
  }
  return(parm)
}

constraintfunc_driver26_alpha = function(parm, Data){
  ## Constrainits on alpha from Driver EBL paper ?
  uniage = cosdistUniAgeAtz(10^parm[grep('z', Data$parm.names)], ref = "737")
  if(parm[grep('mpeak', Data$parm.names)] > uniage){
    parm[grep('mpeak', Data$parm.names)] = uniage ## peak of star formation cannot be before the formation epoch of the galaxy (so older than the age of the Universe at the observed redshift)
  }
  
  if(parm[grep('tau_screen', Data$parm.names)] > parm[grep('tau_birth', Data$parm.names)]){ ## V-band optical depth in the ISM must be less than the birth clouds
    parm[grep('tau_screen', Data$parm.names)] = parm[grep('tau_birth', Data$parm.names)]
  }
  
  parm[grep('alpha_SF_screen', Data$parm.names)] = pmax(0.0625, 2.2 - 0.28*10^parm[grep('z', Data$parm.names)])
  parm[grep('alpha_SF_birth', Data$parm.names)] = pmax(0.0625, 2.2 - 0.28*10^parm[grep('z', Data$parm.names)])
  
  return(parm)
}

save(
  todaysdate,
  frame_info, filters_info, filter_transforms, 
  smiles_miri_nirspec,
  phot_cat, 
  pfunc, pfunc_driver26_alpha, 
  constraintfunc, constraintfunc_driver26_alpha, 
  eazy_zPDFs,
  file = "~/Documents/Epochs/data/prospect_data_for_fitting.rda"
)

save(
  todaysdate,
  frame_info, filters_info, filter_transforms, 
  smiles_miri_nirspec, 
  phot_cat,
  pfunc, pfunc_driver26_alpha, 
  constraintfunc, constraintfunc_driver26_alpha, 
  file = "~/Documents/Epochs/source_inject/prospect_data_for_fitting.rda"
)
