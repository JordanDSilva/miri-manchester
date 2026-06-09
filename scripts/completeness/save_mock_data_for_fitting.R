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

## Make a .rda with all of the data I will need for fitting mock photometry

main_dir = "~/Documents/Epochs/"
fname = "smiles"

## load mock photometry 
mock_photometry = data.frame(
  fread(paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv"))
)
mock_photometry_zPDFs = data.frame(
  fread(paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_mock_photometry_eazy_zPDF.csv"))
)
mock_eazy_zPDFs_normals = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_zPDF_normal.csv")
))
mock_eazy_stats = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_stats.csv")
))

todaysdate = Sys.Date()

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))

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
pfunc_no_tau_birth=function(parm){
  prior = sum(
    c(
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
  mock_photometry, mock_photometry_zPDFs,
  mock_eazy_zPDFs_normals, mock_eazy_stats,
  pfunc, pfunc_driver26_alpha, pfunc_no_tau_birth,
  constraintfunc, constraintfunc_driver26_alpha, 
  file = "~/Documents/Epochs/source_inject/prospect_mock_data_for_fitting.rda"
)
file.copy(
  "~/Documents/Epochs/source_inject/prospect_mock_data_for_fitting.rda",
  "~/Documents/Epochs/data/prospect/prospect_mock_data_for_fitting.rda", 
  overwrite = TRUE
)
file.copy(
  "~/Documents/Epochs/source_inject/prospect_mock_data_for_fitting.rda",
  "~/Documents/Epochs/data/prospect_mock_data_for_fitting.rda", 
  overwrite = TRUE
)