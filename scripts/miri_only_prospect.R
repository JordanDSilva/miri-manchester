library(ProSpect)
library(nloptr)
library(Highlander)
library(ParmOff)
library(celestial)
library(data.table)
library(foreach)
library(doParallel)
library(matrixStats)
library(stringr)
library(scales)
library(RColorBrewer)
library(magicaxis)

## fit miri only sources with my own photometry that I got myself
set.seed(666)

combined_cat = data.frame(fread(
  paste0("~/Documents/Epochs/data/catalogues/smiles_profound_MIRI_matched_EPOCHSv2_MIRI_only_eyeballed_photometry_cat.csv")
))
nircam_south = fread("~/Documents/Epochs/data/catalogues/miri_only_eyeballed_nircam_south.csv")
nircam_east = fread("~/Documents/Epochs/data/catalogues/miri_only_eyeballed_nircam_east.csv")
nircam_combine = data.frame(rbind(
  nircam_east[c(1,2,3,4,5,6,7,8,9), ],
  nircam_south[c(10), ], 
  nircam_east[c(11), ], 
  fill = TRUE
))

filter_transforms_all = readRDS("~/Documents/Epochs/data/throughputs/filt_funcs.rds")
filters_info_all = fread("~/Documents/Epochs/data/throughputs/wavelengths.txt", col.names = c("wavelengths", "filters"))

filters_info = filters_info_all[
  grepl("NIRCAM|MIRI", filters_info_all$filters), 
]
filter_transforms = filter_transforms_all[grepl("NIRCAM|MIRI", filters_info_all$filters)]

df_flux = cbind(
  nircam_combine[, paste0(str_remove_all(filters_info$filters[grepl("NIRCAM", filters_info$filters)], "JWST_NIRCAM_"), "_flux_aper32")],
  combined_cat[, paste0(filters_info$filters[grepl("MIRI", filters_info$filters)], "_flux_aper32")]
)
df_flux_err = cbind(
  nircam_combine[, paste0(str_remove_all(filters_info$filters[grepl("NIRCAM", filters_info$filters)], "JWST_NIRCAM_"), "_flux_aper32_err_5pc")],
  combined_cat[, paste0(filters_info$filters[grepl("MIRI", filters_info$filters)], "_flux_aper32_err_5pc")]
)

do_fit = function(i){
  kk = i
  
  flux_data = data.frame(
    "filter" = filters_info$filters,
    "cenwave" = filters_info$wavelengths,
    "flux" = unname(unlist(df_flux[kk, ])),
    "fluxerr" = unname(unlist(df_flux_err[kk, ]))
  )
  flux_data$flux[flux_data$flux == 0] = NA
  
  magplot(
    NA,
    log = "x",
    ylim = c(min(flux_data$flux, na.rm = TRUE)/1.5, max(flux_data$flux, na.rm = TRUE)*1.5),
    xlim = c(0.4, 30),
    xlab = "Wavelength [Micron]",
    ylab = "Flux [Jy]"
  )
  points(
    filters_info$wavelengths/1e4,
    flux_data$flux,
    col = "red"
  )
  magerr(
    filters_info$wavelengths/1e4,
    flux_data$flux,
    ylo = (flux_data$fluxerr),
    col = "red"
  )
  
  pfunc=function(parm){
    prior = sum(
      c(
        dnorm(parm['tau_birth'],mean=-0.2,sd=0.5,log=TRUE),
        100*erf(-5 * (parm['tau_screen'] - parm['tau_birth'])) - 100, ## Column density of dust in the ISM should be less than that of the birth clouds
        # -20*erf(parm['tau_screen']-2), ## Old Thorne+20 prior
        dnorm(parm['alpha_SF_birth'],mean=2,sd=1,log=TRUE),
        dnorm(parm['alpha_SF_screen'],mean=2,sd=1,log=TRUE),
        # 100*erf(parm['mperiod']+2) - 100, ## Old Thorne+20 prior
        100*erf(-5 * (parm['mpeak'] - cosdistUniAgeAtz(10^parm['z'], ref = '737'))) - 100, ## disfavour mpeak above the age of the Uni at that epoch
        0
      ),
      na.rm = TRUE
    )
    return(prior)
  }
  
  Data = list(
    flux=flux_data,
    mode = "photom",
    arglist=list(
      massfunc=massfunc_snorm_trunc,
      Z=Zfunc_massmap_lin,
      
      emission = T,
      photoz = T,
      IGMfunc = "Inoue14",
      
      mtrunc = 2,
      
      AGNrm = 60,
      AGNbe = -0.5,
      AGNal = 4.0,
      AGNct = 100,
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7
    ),
    speclib=BC03hr,
    Dale=Dale_NormTot,
    filtout=filter_transforms,
    
    SFH=SFHfunc,
    AGN = Fritz,
    Dale_M2L_func=function(alpha_SF){Dale_M2L_variableDTH_func(alpha_SF = alpha_SF, qPAH_VSG = 0.05)},
    parm.names=c('z', 'mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen','Zfinal','AGNlum', 'AGNan', 'AGNta'), 
    mon.names=c("LP","dustmass.birth", "dustmass.screen", "dustmass.total", "dustlum.birth", "dustlum.screen", "dustlum.total", "SFRburst", paste("flux.",filters_info$filters,sep='')),
    logged=c(
      T,
      T, F, T, F,
      T,T,F,F,
      T,
      T, F, T
    ),
    intervals=list(
      lo=c(
        -2,
        -3, -2, log10(0.3), -0.5,
        -2.5, -5, 0, 0, 
        -4, 
        35, 0.001, -1
      ), 
      hi=c(
        log10(40), 
        4, 13.46728, 2, 1, 
        1.5, 1, 4, 4, 
        -1.3, 
        49, 89.990, 1
      )
    ), 
    fit = 'optim',
    N=dim(flux_data)[1],
    prior = pfunc,
    constraints = NULL,
    verbose=FALSE
  )
  
  constraintfunc = function(parm){
    uniage = cosdistUniAgeAtz(10^grep('z', Data$parm.names), ref = "737")
    if(parm[grep('mpeak', Data$parm.names)] > uniage){
      parm[grep('mpeak', Data$parm.names)] = uniage ## peak of star formation cannot be before the formation epoch of the galaxy (so older than the age of the Universe at the observed redshift)
    }
    
    if(parm[grep('tau_screen', Data$parm.names)] > parm[grep('tau_birth', Data$parm.names)]){ ## V-band optical depth in the ISM must be less than the birth clouds
      parm[grep('tau_screen', Data$parm.names)] = parm[grep('tau_birth', Data$parm.names)]
    }
    return(parm)
  }
  Data$constraints = constraintfunc
  
  startpoint = (Data$intervals$lo+Data$intervals$hi)/2
  names(startpoint) = Data$parm.names
  startpoint["AGNlum"] = 44
  startpoint["z"] = log10(15) ## start point in redshift that most galaxies are probably close to z=1 ? 
  
  Data$fit = 'optim'
  # opt = optim(
  #   par = startpoint,
  #   fn = ProSpectSEDlike,
  #   Data = Data,
  #   method = "Nelder-Mead"
  # )
  opt=nloptr(
    x0 = startpoint,
    eval_f = ProSpectSEDlike,
    Data = Data,
    opts = list(
      "algorithm" ="NLOPT_LN_COBYLA",
      "xtol_rel" = 1.0e-8
    )
  )
  
  # print(opt$par)
  Data$fit = 'check'
  bestfit_opt = ProSpectSEDlike(opt$solution, Data=Data)
  plot(bestfit_opt, ylim = c(1e-15, 1e-1))
  # plot(bestfit_opt$SEDout)
  
  Data$fit = "LD"
  # start = Sys.time()
  highout = Highlander(
    bestfit_opt$parm, 
    Data, 
    ProSpectSEDlike,
    Niters = c(1000, 1000), 
    NfinalMCMC = 1000,
    optim_iters = 2, 
    seed = 666,
    lower=Data$intervals$lo, 
    upper=Data$intervals$hi
  )
  # message(Sys.time() - start)
  
  Data$fit = 'check'
  bestfit = ProSpectSEDlike(highout$par, Data=Data)
  plot(bestfit, ylim = c(1e-10, 1e-1))
  plot(bestfit$SEDout, ylim = c(1e1, 1e10))
  # 10^bestfit$parm['z']
  ret_ = list(
    "bestfit" = bestfit,
    "highout" = highout
  )
}

start = Sys.time()
registerDoParallel(cores = 4)
foobar = foreach(i = 1:4, .combine = c) %dopar% {
  do_fit(i)
}
stopImplicitCluster()
print(Sys.time() - start)

start = Sys.time()
plan(multicore, workers = 4)
foobar = future_map_dbl(
  1:4,
  do_fit, 
  .options = furrr_options(seed = TRUE)
)
plan(sequential)
print(Sys.time() - start)

# 1.4826 * median(
#   abs(
#     ( (foobar - nirspec_z) - median(foobar - nirspec_z) ) / (1 + nirspec_z)
#   )
# )

# plot(bestfit, ylim = c(1e-10, 1e-1))
# plot(bestfit$SEDout)
# 
parm = bestfit$parm
# parm[bestfit$Data$logged] = 10^parm[bestfit$Data$logged]
# 
post = highout$LD_last$Posterior1
# # post[, bestfit$Data$logged] = 10^post[, bestfit$Data$logged]
# 
png("~/Desktop/3.png", width = 20, height = 20, units = "in", res = 240)
par(mar = c(3.5, 3.5, 3.5, 3.5), cex = 0.5)
magtri(post, cex = 0.5)
dev.off()
# 
# plot(
#   density(10^post[,'z'])
# )
# 
sed_samples = as.matrix(foreach(i = 1:dim(post)[1], .combine = rbind) %do% {
  sample_ = ProSpectSEDlike(post[i,], Data=Data)
  sample_$SEDout$FinalFlux$flux
})

pdf("~/Documents/Epochs/plots/test_sed_fit.pdf", width = 7, height = 7)
magplot(
  NA,
  log = "xy",
  ylim = c(1e-7, max(flux_data$flux, na.rm = TRUE)*1.5),
  xlim = c(0.4, 30),
  xlab = "Wavelength [Micron]",
  ylab = "Flux [Jy]"
)
# lines(
#   bestfit$SEDout$FinalFlux$wave/1e4,
#   bestfit$SEDout$FinalFlux$flux,
#   lwd = 3
# )
magerr(
  bestfit$SEDout$FinalFlux$wave/1e4,
  colQuantiles(sed_samples, probs = 0.5),
  ylo = rowDiffs(colQuantiles(sed_samples, probs = c(0.5,0.16))),
  yhi = rowDiffs(colQuantiles(sed_samples, probs = c(0.84,0.5))),
  col = alpha('grey', 0.5),
  poly = TRUE
)
points(
  filters_info$wavelengths/1e4,
  flux_data$flux,
  col = "red"
)
magerr(
  filters_info$wavelengths/1e4,
  flux_data$flux,
  ylo = (flux_data$fluxerr),
  col = "red"
)
legend(
  x = "bottomright",
  legend = c(paste0("zphot = ",round(10^bestfit$parm["z"],3)), paste0("zspec = ",round(nirspec_z[1],3)))
)

plot(
  bestfit$SEDout, 
  ylim = c(1e3, 1e10)
)
dev.off()




## need these for posterior calculations 
# mass_func_args_idx = bestfit$Data$parm.names %in% names(formals(bestfit$Data$arglist$massfunc))
# massfunc_args = parm[mass_func_args_idx]
# 
# do.call('integrate', c(list(
#   f = bestfit$Data$arglist$massfunc, lower = 0, upper = 1e7
# ), massfunc_args))$value / 1e7
# 
# SMstar = ParmOff(.func = SMstarfunc, #the function we want to run
#                           .args = c(parm, bestfit$Data$arglist), #the superset of potential matching parameters
#                           massfunc = bestfit$Data$arglist$massfunc,
#                           speclib = BC03hr,
#                           Z = Zfunc_massmap_lin
# )

.UV_tophat = function(wave){
  throughput = rep(0L, length(wave))
  throughput[wave >= (1500 - 50) & wave <= (1500 + 50)] = 1L
  return(throughput)
}
.UVLum1500 = function(out){
  out$Data$fit = "check"
  
  UV1500 = photom_lum(
    wave = out$SEDout$FinalLum$wave,
    lum = out$SEDout$FinalLum$lum,
    z = 0.0, 
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref,
    filters = list(.UV_tophat)
  )
  
  return(UV1500)
}
.UVLumAGN1500 = function(out){
  out$Data$fit = "check"
  
  UV1500 = photom_lum(
    wave = out$SEDout$AGN$wave,
    lum = out$SEDout$AGN$lum,
    z = 0.0, 
    LumDist_Mpc = out$Data$arglist$LumDist_Mpc,
    ref = out$Data$arglist$ref,
    filters = list(.UV_tophat)
  )
  
  return(UV1500)
}
.stellar_mass = function(out){
  
  out$Data$fit = "check"
  
  parm = out$parm
  parm[out$Data$logged] = 10^parm[out$Data$logged]
  
  SMstar = ParmOff(
    .func = SMstarfunc, #the function we want to run
    .args = c(parm, out$Data$arglist), #the superset of potential matching parameters
    massfunc = out$Data$arglist$massfunc,
    speclib = BC03hr,
    Z = Zfunc_massmap_lin
  )
  return(as.numeric(SMstar['TotSMstar']))
}
.SFR10 = function(out){
  
  out$Data$fit = "check"
  
  parm = out$parm
  parm[out$Data$logged] = 10^parm[out$Data$logged]
  
  mass_func_args_idx = out$Data$parm.names %in% names(formals(out$Data$arglist$massfunc))
  massfunc_args = c(parm[mass_func_args_idx], "magemax" = out$Data$arglist$magemax)
  
  sfr10 = do.call('integrate', c(list(
    f = out$Data$arglist$massfunc, lower = 0, upper = 1e7
  ), massfunc_args))$value / 1e7
  
  return(sfr10)
}

nparm = length(bestfit$Data$parm.names)
post = cbind(highout$LD_last$Posterior1, highout$LD_last$Monitor)
post = data.table(post)
post = post[order(-LP),]

## calculating a chisq cut based on the number of parameters to isolate a 1sigma range
chisq_cut = max(post$LP) - qchisq(0.68,nparm)/2
toppost = post[post$LP > chisq_cut,]
toppost[,LP:=NULL]

registerDoParallel(cores = 8)
SED_samples = foreach(i = 1:dim(toppost)[1], .combine = rbind) %dopar% {
  if(i %% 100 == 0){
    message(i)
  }
  
  out = ProSpectSEDlike(highout$LD_last$Posterior1[i,], Data)
  
  z = 10^as.numeric(highout$LD_last$Posterior1[i,'z'])
  stellar_mass = .stellar_mass(out)
  sfr10 = .SFR10(out)
  UV1500 = .UVLum1500(out)
  c(
    "z" = z, 
    "UV1500" = UV1500,
    "StellarMass" = stellar_mass,
    "SFR10" = sfr10
  )
}
stopImplicitCluster()

SED_samples = data.frame(SED_samples)

foobar = SED_samples
foobar$StellarMass = log10(foobar$StellarMass)
foobar$SFR10 = log10(foobar$SFR10)
foobar$StellarMass[is.infinite(foobar$StellarMass)] = 0
foobar$SFR10[is.infinite(foobar$SFR10)] = 0

magtri(
  as.matrix(foobar[foobar$StellarMass >= 6 & foobar$z <= 10, ])
)
