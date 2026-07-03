library(ProFound)
library(ProSpect)
library(ProSpectUtils)
library(ProFit)
library(Rfits)
library(celestial)
library(magicaxis)
library(data.table)
library(pracma)
library(MASS)
library(foreach)
library(doParallel)
library(matrixStats)
library(nloptr)
library(Highlander)
library(furrr)
library(purrr)
library(stringr)
library(minpack.lm)
library(numDeriv)

# main_dir = "~/Documents/Epochs/source_inject/"
main_dir = "~/mnt/morgan/raid/"
fname = "smiles"

profound_cat = Rfits_read(paste0("~/Documents/Epochs/data/catalogues/", fname, "_profound_multiband_cat.fits"))
# det_image = Rfits_read(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"))
# pro_detect = readRDS(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound.rds"))

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/smiles_survey/smiles_survey_mock_cat.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/smiles_survey/smiles_survey_mock_eazy_photometry.csv")
))

Nsample = 2000

Re_z = function(z){
  # From Ormerod+26 in kpc then convert to arsec
  rnorm(length(z), 4.50, 1.32)*(1+z)^rnorm(length(z), -0.71, 0.19) / (cosdistAngScale(z = z, ref = "737"))
}

source_inject = function(main_dir, Nsample, astro_cat, fname = 'smiles'){
  
  frame_dir = paste0(
    main_dir, "data/JP/ProFound/Data/", fname, "/", fname, "/"
  )
  detect_dir = paste0(
    main_dir, "data/JP/ProFound/Detects/", fname, "/", fname, "/"
  )
  
  ## make the fake sources
  det_F560W_flux_true = as.numeric(phot_cat$F560W_true_flux)
  det_F770W_flux_true = as.numeric(phot_cat$F770W_true_flux)
  
  det_mag_true = -2.5*log10(colMeans(
    rbind(det_F560W_flux_true, det_F770W_flux_true)
  ))+8.9
  
  eligible_idx <- which(det_mag_true > 10 & det_mag_true < 31)
  
  sampled_idx <- sample(
    x    = eligible_idx,
    size = Nsample
  )
  
  mag_df = -2.5*log10(phot_cat[sampled_idx, grepl("true", names(phot_cat))])+8.9
  
  astro_cat_detect_sample <- astro_cat[sampled_idx, ]
  
  detect_image = Rfits_read(
    paste0(detect_dir, fname, "_", fname, "_", fname, "_profound_stack.fits")
  )
  
  non_na_coords = which(!is.na(detect_image$image[,]$imDat), arr.ind = T) ## populate fake galaxies in observed image
  sample_idx = sample(1:dim(non_na_coords)[1], size = Nsample)
  xcen = non_na_coords[sample_idx, 1]
  ycen = non_na_coords[sample_idx, 2]
  samples_df$xcen = xcen
  samples_df$ycen = ycen
  
  synthetic_mags = det_mag_true[sampled_idx]
  synthetic_R50 = Re_z(astro_cat_detect_sample$z)
  samples_df = list()
  samples_df$mag <- synthetic_mags
  samples_df$R50 <- synthetic_R50
  samples_df$nser = rep(1, Nsample)
  samples_df$ang = runif(Nsample, 0.01, 89.99)
  samples_df$axrat = runif(Nsample, 0.2, 1.0)
  
  source_list=list(
    sersic = list(
      xcen = xcen,
      ycen = ycen,
      mag = samples_df$mag,
      re = samples_df$R50/pixscale(det_image$image),
      nser = samples_df$nser,
      ang = samples_df$ang,
      axrat = samples_df$axrat,
      box = 0
    ),
    sky = list(bg = 0)
  )
  
  ## Use profit to inject sources into frame
  start = Sys.time()
  model_sources = profitMakeModel(
    source_list,
    dim = dim(det_image$image), 
    magzero = det_image$image$keyvalues$MAGZERO, 
    rough = TRUE
  )
  print(Sys.time() - start)
  
  
  fits_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources.fits")
  mag_samples_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_mag_samples.rds")
  pro_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_profound.rds")
  phot_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_photometry.rds")
  
  # obs_mag_density <- density(profound_cat$ProFound$mag, bw = "nrd0")
  obs_mag_density <- density(c(profound_cat$ProFound$mag, runif(2000, 28, 40)), bw = "nrd0")
  
  # Normalise to get a probability
  p_obs <- approxfun(obs_mag_density$x, obs_mag_density$y / max(obs_mag_density$y))
  
  kde2 <- kde2d(
    profound_cat$ProFound$mag, 
    log10(profound_cat$ProFound$R50), 
    n = 200,
    lims = c(range(profound_cat$ProFound$mag), range(log10(profound_cat$ProFound$R50)))
  )
  
  det_F560W_flux_true = as.numeric(phot_cat$F560W_true_flux)
  det_F770W_flux_true = as.numeric(phot_cat$F770W_true_flux)
  
  det_F560W_flux_test = as.numeric(phot_cat$F560W_flux)
  det_F770W_flux_test = as.numeric(phot_cat$F770W_flux)
  
  det_F560W_flux_err = as.numeric(phot_cat$F560W_flux_err)
  det_F770W_flux_err = as.numeric(phot_cat$F770W_flux_err)
  
  det_mag_true = -2.5*log10(colMeans(
    rbind(det_F560W_flux_true, det_F770W_flux_true)
  ))+8.9
  det_mag_test = -2.5*log10(colMeans(
    rbind(pmax(0,det_F560W_flux_test), pmax(0,det_F770W_flux_test)), 
  ))+8.9
  
  eligible_idx <- which(det_mag_true > 10 & det_mag_true < 31)
  
  sampled_idx <- sample(
    x    = eligible_idx,
    size = Nsample
  )
  
  astro_cat_detect_sample <- astro_cat[sampled_idx, ]
  
  synthetic_mags = det_mag_true[sampled_idx]
  
  synthetic_R50 <- mapply(sample_R50_given_mag,
                          mag_val = synthetic_mags,
                          MoreArgs = list(kde2 = kde2, n = 1))
  
  samples_df = list()
  samples_df$mag <- synthetic_mags
  samples_df$R50 <- synthetic_R50
  samples_df$nser = rep(1, Nsample)
  samples_df$ang = runif(Nsample, 0.01, 89.99)
  samples_df$axrat = runif(Nsample, 0.2, 1.0)
  
  non_na_coords = which(!is.na(det_image$image[,]$imDat), arr.ind = T) ## populate fake galaxies in observed image
  sample_idx = sample(1:dim(non_na_coords)[1], size = Nsample)
  xcen = non_na_coords[sample_idx, 1]
  ycen = non_na_coords[sample_idx, 2]
  samples_df$xcen = xcen
  samples_df$ycen = ycen
  
  source_list=list(
    sersic = list(
      xcen = xcen,
      ycen = ycen,
      mag = samples_df$mag,
      re = samples_df$R50/pixscale(det_image$image),
      nser = samples_df$nser,
      ang = samples_df$ang,
      axrat = samples_df$axrat,
      box = 0
    ),
    sky = list(bg = 0)
  )
  
  ## Use profit to inject sources into frame
  start = Sys.time()
  model_sources = profitMakeModel(
    source_list,
    dim = dim(det_image$image), 
    magzero = det_image$image$keyvalues$MAGZERO, 
    rough = TRUE
  )
  print(Sys.time() - start)
  
  fake_sources = Rfits_create_image(
    data = model_sources$z + det_image$image[,]$imDat, 
    keyvalues = det_image$image[,]$keyvalues
  )
  plot(fake_sources)
  
  mask = pro_detect$mask == 1L | pro_detect$objects_redo == 1L
  pro_fake = profoundProFound(
    image = fake_sources,
    mask = mask,
    
    # skyRMS = pro_detect$skyRMS,
    magzero = fake_sources$keyvalues$MAGZERO,
    
    sky = 0,
    skycut = 1.5, 
    pixcut = 7.0,
    
    tolerance = 4.0,
    reltol = 0,
    cliptol = 30,
    
    iters = 4,
    
    sigma = 1.0,
    box = 100,
    grid = 100,
    boxiters = 0,
    
    roughpedestal = FALSE,
    pixelcov = FALSE,
    boundstats = TRUE,
    groupstats = TRUE,
    nearstats = TRUE, 
    redosky = TRUE,
    rem_mask = TRUE,
    # redoskysize = 11,
    
    fluxtype = "Jansky",
    pixscale = pixscale(fake_sources, unit = "asec"),
    
    verbose = FALSE,
  )
  
  detect_segstats = pro_fake$segstats
  rm(pro_fake)
  gc()
  
  samples_coords = cbind(Rwcs::Rwcs_p2s(cbind(samples_df$xcen, samples_df$ycen), keyvalues = fake_sources$keyvalues))
  match_idx = celestial::coordmatch(
    coordref = samples_coords,
    coordcompare = detect_segstats[,c("RAmax", "Decmax")], 
    rad = 0.3, 
    radunit = "asec"
  )
  
  hh_true = maghist(
    samples_df$mag, plot = FALSE
  )
  hh_true_detect = maghist(
    samples_df$mag[match_idx$bestmatch$refID], plot = FALSE, breaks = hh_true$breaks
  )
  hh_detect = maghist(
    detect_segstats$mag[match_idx$bestmatch$compareID], plot = FALSE, breaks = hh_true$breaks
  )
  
  hh_true = maghist(
    log10(astro_cat_detect_sample$StellarMass), plot = FALSE
  )
  hh_true_detect = maghist(
    log10(astro_cat_detect_sample$StellarMass)[match_idx$bestmatch$refID] + (detect_segstats$mag[match_idx$bestmatch$compareID] - samples_df$mag[match_idx$bestmatch$refID]), plot = FALSE, breaks = hh_true$breaks
  )
  hh_detect = maghist(
    detect_segstats$mag[match_idx$bestmatch$compareID], plot = FALSE, breaks = hh_true$breaks
  )
  
  magplot(
    hh_true$mids,
    hh_true_detect$counts/hh_true$counts, 
    ylim = c(0,1)
  )
  
  out = list(
    samples_df = samples_df,
    detect_segstats = detect_segstats,
    matched_segstats = detect_segstats[match_idx$bestmatch$compareID, ],
    match_idx = match_idx
  )
  return(
    out
  )
}

out = source_inject(det_image, profound_cat, pro_detect, Nsample, fname = 'smiles')

magplot(
  NA, 
  xlim = c(10, 40),
  ylim = c(0,1),
  type = "l",
  xlab = "Magnitude (recovered)",
  ylab = "Completeness",
)

mag_in = out$samples_df$mag
mag_rec = out$detect_segstats$mag
mag_match = out$matched_segstats$mag
good_recovery = out$detect_segstats$segID %in% out$detect_segstats$segID[out$match_idx$bestmatch$compareID]

mag_bins <- seq(10, 40, 0.5)
hh1 = maghist(
  mag_match, 
  breaks = mag_bins,
  plot = FALSE
)
hh2 = maghist(
  mag_in, 
  breaks = mag_bins,
  plot = FALSE
)

bin_df <- data.frame(
  mag_mid = hh1$mids,
  comp    = hh1$counts / hh2$counts,
  err     = (hh1$counts / hh2$counts) * sqrt(1/hh1$counts + 1/hh2$counts)
)

bin_df = bin_df[!is.na(bin_df$comp) & !is.infinite(bin_df$comp) & !is.na(bin_df$err), ]

comp_nls <- nlsLM(
  comp ~ C_max / (1 + exp(-b1 * (mag_mid - m50))),
  data    = bin_df,
  start   = list(C_max = 0.7, b1 = -1, m50 = 26),
  lower   = c(C_max = 0, b1 = -5.5, m50 = 10),
  upper   = c(C_max = 1, b1 = 0,   m50 = 40),
  weights = 1 / bin_df$err^2
)

coef(comp_nls)

mag_seq  <- seq(10, 40, 0.1)
comp_pred <- predict(comp_nls, newdata = data.frame(mag_mid = mag_seq))

plot(bin_df$mag_mid, bin_df$comp, pch = 16,
     xlab = "Magnitude (recovered)", ylab = "Completeness",
     ylim = c(0, 1))
arrows(bin_df$mag_mid, bin_df$comp - bin_df$err,
       bin_df$mag_mid, bin_df$comp + bin_df$err,
       length = 0.02, angle = 90, code = 3)
lines(mag_seq, predict(comp_nls, newdata = data.frame(mag_mid = mag_seq)),
      col = "steelblue", lwd = 2)

# get variance-covariance matrix of parameters
vcov_mat <- vcov(comp_nls)

# predicted values
comp_pred <- predict(comp_nls, newdata = data.frame(mag_mid = mag_seq))

pred_func <- function(params) {
  params["C_max"] / (1 + exp(params["b1"] * (mag_seq - params["m50"])))
}

jac <- jacobian(pred_func, coef(comp_nls))

# variance of prediction at each mag point
pred_var <- rowSums((jac %*% vcov_mat) * jac)
pred_sd  <- sqrt(pred_var)

# 16th and 84th percentile bands (1 sigma)
pred_lo <- comp_pred - pred_sd
pred_hi <- comp_pred + pred_sd

# plot
plot(bin_df$mag_mid, bin_df$comp, pch = 16,
     xlab = "Magnitude (recovered)", ylab = "Completeness",
     xlim = c(10, 40),
     ylim = c(0, 1))
arrows(bin_df$mag_mid, bin_df$comp - bin_df$err,
       bin_df$mag_mid, bin_df$comp + bin_df$err,
       length = 0.02, angle = 90, code = 3)

# shaded confidence band
polygon(c(mag_seq, rev(mag_seq)),
        c(pred_lo, rev(pred_hi)),
        col = adjustcolor("steelblue", alpha.f = 0.2),
        border = NA)

# main curve
lines(mag_seq, comp_pred, col = "steelblue", lwd = 2)

# 16/84 lines
lines(mag_seq, pred_lo, col = "steelblue", lwd = 1, lty = 2)
lines(mag_seq, pred_hi, col = "steelblue", lwd = 1, lty = 2)

# m50 and C_max reference lines
abline(v = coef(comp_nls)["m50"],          lty = 2, col = "red")
abline(h = coef(comp_nls)["C_max"] / 2,    lty = 2, col = "grey50")
abline(h = coef(comp_nls)["C_max"],         lty = 2, col = "grey50")

saveRDS(
  list(
    "samples" = out, 
    "nls" = comp_nls,
    "compl_curves" = data.frame(
      "mag" = mag_seq,
      "predQ50" = comp_pred,
      "predQ16" = comp_pred - pred_sd,
      "predQ84" = comp_pred + pred_sd
    ),
    "compl_fn" = approxfun(
      x = mag_seq,
      y = comp_pred,
      rule = 2
    )
  ),
  paste0(main_dir, fname, "/", fname, "_mock_1detection_samples.rds")
)



# compl_repeats = lapply(1:10, function(x){source_inject(det_image, profound_cat, pro_detect, Nsample, fname = 'smiles')})
# compl_repeats[[1]]

# magbins = seq(10, 40, 1.0)
# cols = hcl.colors(10, "Set 2")
# magplot(NA, xlim = c(10, 40), ylim = c(0,2.0))
# for(i in 1:10){
#   hh_true = maghist(
#     compl_repeats[[i]]$samples_df$mag,
#     breaks = magbins,
#     plot = FALSE
#   )
#   hh_detect = maghist(
#     compl_repeats[[i]]$matched_segstats$mag,
#     breaks = magbins,
#     plot = FALSE
#   )
#   points(
#     hh_true$mids,
#     hh_detect$counts/hh_true$counts, 
#     col = cols[i]
#   )
# }
# mag_true = unlist(lapply(1:10, function(i){compl_repeats[[i]]$samples_df$mag}))
# mag_detect = unlist(lapply(1:10, function(i){compl_repeats[[i]]$matched_segstats$mag}))
# mag_detect_true = unlist(lapply(1:10, function(i){compl_repeats[[i]]$samples_df$mag[compl_repeats[[i]]$match_idx$bestmatch$refID]}))
# 
# hh_true = maghist(
#   mag_true,
#   breaks = magbins,
#   plot = FALSE
# )
# hh_detect = maghist(
#   mag_detect,
#   breaks = magbins,
#   plot = FALSE
# )
# hh_detect_true = maghist(
#   mag_detect_true,
#   breaks = magbins,
#   plot = FALSE
# )
# points(
#   hh_true$mids,
#   hh_detect$counts/hh_true$counts, 
#   pch = 16
# )
# points(
#   hh_true$mids,
#   (hh_detect$counts-hh_detect_true$counts)/hh_true$counts, 
#   pch = 15
# )
# magerr(
#   hh_true$mids,
#   hh_detect$counts/hh_true$counts,
#   sqrt( hh_detect$counts/hh_true$counts^2 * ((sqrt(hh_detect$counts)/hh_detect$counts)^2 + (sqrt(hh_true$counts)/hh_true$counts)^2) )
# )