#scratch 

do_calc = function(i){
  kk = i
  save_stub = paste0(
    main_dir, "/prospect/photoZ_bc03_snorm_fritz/galID_",kk, ".rds" 
  )
  
  out = readRDS(save_stub)
  
  Data = out$bestfit$Data
  Data$speclib = BC03lr
  Data$AGN = Fritz
  Data$Dale = Dale_NormTot
  
  nparm = length(out$bestfit$parm)
  post = cbind(out$highout$LD_last$Posterior1, out$highout$LD_last$Monitor)
  post = data.table(post)
  post = post[order(-LP),]
  
  ## calculating a chisq cut based on the number of parameters to isolate a 1sigma range
  chisq_cut = max(post$LP) - qchisq(0.68,nparm)/2
  toppost = post[post$LP > chisq_cut,]
  toppost[,LP:=NULL]
  
  flux_samp = matrix(0, nrow = dim(toppost)[1], ncol = length(out$bestfit$SEDout$FinalFlux$wave))
  flux_starsUnAtten_samp = matrix(0, nrow = dim(toppost)[1], ncol = length(out$bestfit$SEDout$FinalFlux$wave))
  flux_starsAtten_samp = matrix(0, nrow = dim(toppost)[1], ncol = length(out$bestfit$SEDout$FinalFlux$wave))
  flux_starsAGN_samp = matrix(0, nrow = dim(toppost)[1], ncol = length(out$bestfit$SEDout$FinalFlux$wave))
  flux_starsDust_samp = matrix(0, nrow = dim(toppost)[1], ncol = length(out$bestfit$SEDout$FinalFlux$wave))
  
  astro = foreach(jj = 1:dim(toppost)[1], .combine = rbind) %do% {
    
    samp = ProSpectSEDlike(parm = unlist(toppost[jj, 1:nparm]), Data = Data)
    
    obs_samp = .fitSEDObs(samp)
    flux_samp[jj,] = obs_samp$flux
    flux_starsUnAtten_samp[jj,] = obs_samp$starsUnAtten
    flux_starsAtten_samp[jj,] = obs_samp$starsAtten
    flux_starsAGN_samp[jj,] = obs_samp$AGN
    flux_starsDust_samp[jj,] = obs_samp$Dust
    
    df = c(
      "StellarMass" = .stellar_mass(out = samp),
      "SFR10" = .SFR10(out = samp),
      "SFR100" = samp$SEDout$Stars$SFRburst,
      "UV1500" = .UVLum1500(out = samp),
      "UV1500AGN" = .UVLumAGN1500(out = samp),
      "DustLumBirth" = as.numeric(samp$SEDout$dustlum['birth']),
      "DustLumScreen" = as.numeric(samp$SEDout$dustlum['screen']),
      "DustLumTotal" = as.numeric(samp$SEDout$dustlum['total']),
      "LP" = samp$LP
    )
    extra_parm = samp$parm 
    extra_parm[samp$Data$logged] = 10^extra_parm[samp$Data$logged]
    c(extra_parm, df)
  }
  
  ## astrophysics
  out$astroPost = astro
  
  astro_quantiles = colQuantiles(
    as.matrix(astro), probs = c(0.5, 0.16, 0.84)
  )
  params = rownames(astro_quantiles)
  
  summary_dt = data.table(
    Parameter = params,
    Q50       = astro_quantiles[, "50%"],
    Q16       = astro_quantiles[, "16%"],
    Q84       = astro_quantiles[, "84%"],
    Lo        = astro_quantiles[, "50%"] - astro_quantiles[, "16%"],  # lower uncertainty
    Hi        = astro_quantiles[, "84%"] - astro_quantiles[, "50%"]   # upper uncertainty
  )
  out$astro = summary_dt
  
  ## astronomy
  fitted_sed = lapply(
    list(
      "flux" = flux_samp,
      "starsUnAtten" = flux_starsUnAtten_samp,
      "starsAtten" = flux_starsAtten_samp,
      "AGN" = flux_starsAGN_samp,
      "Dust" = flux_starsDust_samp
    ),
    function(x){
      ret = colQuantiles(as.matrix(x), probs = c(0.50, 0.16, 0.84))
      
    }
  )
  out$fittedSED = fitted_sed
  
  saveRDS(
    out,
    save_stub
  )
}


Data$fit = 'optim'
opt = nloptr(
  x0 = startpoint,
  eval_f = ProSpectSEDlike,
  lb = Data$intervals$lo,
  ub = Data$intervals$hi,
  Data = Data,
  opts = list(
    "algorithm"   = "NLOPT_LN_BOBYQA",   # faster than COBYLA with bounds
    "xtol_rel"    = 1.0e-6,              # loosen slightly — 1e-8 is overkill for SED fitting
    "ftol_rel"    = 1.0e-6,              # stop when objective stops improving too
    "ftol_abs"    = 1.0e-6,             # absolute objective tolerance (log-likelihood units)
    "maxeval"     = 10000               # cap iterations so it doesn't run forever
  )
)