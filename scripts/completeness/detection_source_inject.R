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

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles"

profound_cat = Rfits_read(paste0("~/Documents/Epochs/data/catalogues/", fname, "_profound_multiband_cat.fits"))
det_image = Rfits_read(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"))
pro_detect = readRDS(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound.rds"))

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/smiles_survey/smiles_survey_mock_cat.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/smiles_survey/smiles_survey_mock_eazy_photometry.csv")
))

dir.create(
  paste0(main_dir, fname, "/")
)

load(
  paste0(
    main_dir, "/prospect_mock_data_for_fitting.rda"
  )
)

Nsample = 10000

sample_R50_given_mag <- function(mag_val, kde2, n = 1) {
  
  # Find nearest mag grid column
  
  ix <- which.min(abs(kde2$x - mag_val))
  
  # Extract the conditional density slice P(R50 | mag = mag_val)
  cond_density <- kde2$z[ix, ]          # column = fixed mag, row = R50
  cond_density <- pmax(cond_density, 0) # remove tiny negatives
  cond_density <- cond_density / sum(cond_density)  # normalise to sum=1
  
  # Sample R50 grid points with those probabilities
  idx <- sample(length(kde2$y), size = n, replace = TRUE, prob = cond_density)
  
  # Add jitter within half a grid cell
  dy  <- diff(kde2$y)[1] / 2
  10^(kde2$y[idx] + runif(n, -dy, dy))
}

source_inject = function(det_image, profound_cat, pro_detect, Nsample, astro_cat, fname = 'smiles'){
  
  fits_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources.fits")
  mag_samples_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_mag_samples.rds")
  pro_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_profound.rds") 
  phot_stub = paste0(main_dir, "/", fname, "/", fname, "_fake_sources_photometry.rds")
  
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
  phot_cat_detect_sample <- phot_cat[sampled_idx, ]
  
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
  
  fwrite(
    phot_cat_detect_sample[match_idx$bestmatch$refID, ], 
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_uniform_detect_mock_eazy_photometry_recovered.csv")
  )
  fwrite(
    astro_cat_detect_sample[match_idx$bestmatch$refID, ], 
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_uniform_detect_mock_cat_recovered.csv")
  )
  
  fwrite(
    phot_cat_detect_sample, 
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_uniform_detect_mock_eazy_photometry.csv")
  )
  fwrite(
    astro_cat_detect_sample, 
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_uniform_detect_mock_cat.csv")
  )

  out = list(
    samples_df = samples_df,
    detect_segstats = detect_segstats,
    matched_segstats = detect_segstats[match_idx$bestmatch$compareID, ],
    match_idx = match_idx,
    model_sources = model_sources
  )
  
  saveRDS(
    out,
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_uniform_detect.rds")
  )
  Rfits_write_image(
    fake_sources, 
    paste0(main_dir, "/", fname, "_uniform_detect/", fname, "_injection_image.fits")
  )
  return(
    out
  )
}

out = source_inject(det_image, profound_cat, pro_detect, Nsample, astro_cat = astro_cat, fname = 'smiles')
