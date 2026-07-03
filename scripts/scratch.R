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



# completeness crap
## Same method as 2 from Leethochawalit et al 2022 MNRAS
mag_bins = seq(10, 40, 0.5)
input_mag_hist = maghist(
  x = samples_df$mag, 
  breaks = mag_bins, 
  plot = F, 
  verbose = F
)

input_match_mag_hist = maghist(
  # x = samples_df$mag[match_idx$bestmatch$refID],
  x = pro_fake$segstats$mag[match_idx$bestmatch$compareID],
  breaks = mag_bins, 
  plot = F, 
  verbose = F
)

input_mag_contaminants_hist = maghist(
  # x = samples_df$mag[match_idx$bestmatch$refID],
  x = pro_fake$segstats$mag[!(1:dim(pro_fake$segstats)[1] %in% match_idx$bestmatch$compareID)],
  breaks = mag_bins, 
  plot = F, 
  verbose = F
)
compl_mag = (input_match_mag_hist$counts)/(input_mag_hist$counts)
compl_pois = sqrt( compl_mag^2 * 
                     ( (sqrt(input_match_mag_hist$counts)/input_match_mag_hist$counts)^2 + 
                         (sqrt(input_mag_hist$counts)/input_mag_hist$counts)^2 ) )
compl_pois[compl_pois == 0] = 0.01 

compl_contaminants_mag = (input_mag_contaminants_hist$counts)/(input_mag_hist$counts)
compl_contaminants_pois = sqrt( compl_contaminants_mag^2 * 
                                  ( (sqrt(input_mag_contaminants_hist$counts)/input_mag_contaminants_hist$counts)^2 + 
                                      (sqrt(input_mag_hist$counts)/input_mag_hist$counts)^2 ) )
compl_pois[compl_pois == 0] = 0.01 

x = input_mag_hist$mids
# x = c(input_mag_hist$mids, 28.2)
# compl_mag = c(compl_mag, 0)
# compl_pois = c(compl_pois, 1e-99)

## logisitic function
logistic = function(x, parm){
  L = parm[1]
  k = parm[2]
  x0 = parm[3]
  
  yy = L / (1 + exp(-1*k * (x-x0)))
  return(yy)
}
LL = function(parm, Data){
  
  like = dnorm(
    x = Data$y,
    mean = logistic(Data$x, parm),
    sd = Data$yerr,
    log = TRUE
  )
  return(
    sum(like, na.rm = TRUE)
  )
}
opt = optim(
  par = c(0.8, 1, 26),
  fn = LL, 
  Data = list(
    x = c(x[!is.nan(compl_mag)]),
    y = c(compl_mag[!is.nan(compl_mag)]),
    yerr = c(compl_pois[!is.nan(compl_mag)])
  ),
  control = list(fnscale = -1)
)

pdf(
  file = paste0(main_dir, "/", fname, "/", fname, "_detection_completeness.pdf"),
  width = 12, height = 5
)
par(mfrow = c(1,2), oma = rep(1.0, 4), mar = c(2.5, 2.5, 1.5, 1.5))
magplot(
  x,
  compl_mag,
  xlim = c(15, 35),
  ylim = c(0, 1),
  xlab = "Detect magnitude",
  ylab = "Detection completeness"
)
magerr(
  x,
  compl_mag,
  ylo = compl_pois
)
abline(
  v = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9, 
  col = "cornflowerblue", 
  lwd = 2
)
lines(
  x,
  logistic(x, opt$par),
  col = "blue",
  lwd = 2
)
legend(
  x = "topright",
  legend = c("Logisitic fit"),
  col = "blue",
  lty = 1
)

magplot(
  x,
  compl_contaminants_mag,
  xlim = c(15, 35),
  # ylim = c(0, 1),
  xlab = "Detect magnitude",
  ylab = "Detection contamination"
)
magerr(
  x,
  compl_contaminants_mag,
  ylo = compl_contaminants_pois
)
abline(
  v = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9, 
  col = "cornflowerblue", 
  lwd = 2
)

dev.off()

compl_df = data.frame(
  "x" = x,
  "compl_mag" = compl_mag,
  "compl_mag_err" = compl_pois, 
  "compl_contaminants" = compl_contaminants_mag,
  "compl_contaminants_err" = compl_contaminants_pois,
  "compl_fit" = logistic(x, opt$par)
)

fwrite(
  compl_df, 
  paste0(main_dir, "/", fname, "/", fname, "_det_mag_completeness_stats.csv")
)

## From completeness scripts
calc_selection_3d = function(
    zbins, zmids, ztrue, 
    xtrue, xbins, xmids,
    ytrue, ybins, ymids,
    selection_idx,
    ztest = NULL,
    xtest = NULL,
    ytest = NULL
){
  Ndim = c(length(zmids), length(xmids), length(ymids))
  completeness_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  purity_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  
  if(is.null(ztest)){
    ztest = rnorm(length(ztrue), ztrue, 0.06)
  }
  if(is.null(xtest)){
    xtest = rnorm(length(xtrue), xtrue, 0.3)
    # xtest = xtrue
  }
  if(is.null(ytest)){
    ytest = rnorm(length(ytrue), ytrue, 0.3)
    # ytest = ytrue
  }
  
  for(i in seq_along(zmids)){
    for(j in seq_along(xmids)){
      for(k in seq_along(ymids)){
        
        id_array_true = which(
          (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) &
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1])
        )
        
        id_array_sel = which(
          (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) &
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1]) & 
            (selection_idx)
        )
        
        id_array_det_sel = which(
          (ztest >= zbins[i]) & (ztest < zbins[i+1]) &
            (xtest >= xbins[j]) & (xtest < xbins[j+1]) &
            (ytest >= ybins[k]) & (ytest < ybins[k+1]) &
            (selection_idx)
        )
        
        id_in_bin = which(
          abs(ztest[id_array_det_sel] - ztrue[id_array_det_sel])/(1 + ztrue[id_array_det_sel]) <= 0.15 & 
            abs(xtest[id_array_det_sel] - xtrue[id_array_det_sel]) <= (xbins[j+1] - xbins[j]) & 
            abs(ytest[id_array_det_sel] - ytrue[id_array_det_sel]) <= (ybins[k+1] - ybins[k])
        )
        
        id_array_det_sel_exact = intersect(
          id_array_sel, id_array_det_sel
        )
        
        if(length(id_array_true) == 0){
          completeness_array[i,j,k] = NA
          purity_array[i,j,k] = NA
        }else{
          # completeness_array[i,j,k] = pmin(length(id_array_det_sel)/length(id_array_true),1)
          completeness_array[i,j,k] = length(id_in_bin)/length(id_array_true)
          purity_array[i,j,k] = length(id_in_bin)/length(id_array_det_sel)
        }
      }
    }
    
    cont_mat = 1 - as.matrix(purity_array[i,,])
    # cont_mat[cont_mat > 0.5] = 2
    compl_mat = as.matrix(completeness_array[i,,])
    
    par(mfrow = c(1,2))
    image(
      x = xmids,
      y = ymids,
      z = compl_mat,
      xlab = "Mstar",
      ylab = "SFR",
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      xlim = range(xbins),
      ylim = range(ybins), 
    )
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      legend.lab = "Completeness"
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3)), 
        paste0("Min = ", round(min(compl_mat, na.rm = TRUE), 3)),
        paste0("Max = ", round(max(compl_mat, na.rm = TRUE), 3))
      )
    )
    image(
      x = xmids,
      y = ymids,
      z = cont_mat,
      xlab = "Mstar",
      ylab = "SFR",
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      xlim = range(xbins),
      ylim = range(ybins), 
    )
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      legend.lab = "Contamination"
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3)), 
        paste0("Min = ", round(min(cont_mat, na.rm = TRUE), 3)),
        paste0("Max = ", round(max(cont_mat, na.rm = TRUE), 3))
      )
    )
  }
  return( completeness_array )
}
