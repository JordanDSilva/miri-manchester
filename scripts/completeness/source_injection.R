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

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles"

profound_cat = Rfits_read(paste0("~/Documents/Epochs/data/catalogues/", fname, "_profound_multiband_cat.fits"))
det_image = Rfits_read(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound_stack.fits"))
pro_detect = readRDS(paste0("~/Documents/Epochs/data/images/", fname, "/", fname, "_", fname, "_", fname, "_profound.rds"))

dir.create(
  paste0(main_dir, fname, "/")
)

load(
  paste0(
    main_dir, "/prospect_data_for_fitting.rda"
  )
)

Nsample = 10000

source_inject = function(det_image, profound_cat, pro_detect, Nsample, fname = 'smiles'){
  
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

  logged = c(
    TRUE,
    TRUE, FALSE, TRUE, FALSE,
    TRUE, TRUE, FALSE, FALSE,
    TRUE,
    TRUE, FALSE, TRUE
  )
  intervals = list(
    lo = c(
      -2,
      -3, -2, log10(0.3), -0.5,
      -2.5, -5, 0, 0, 
      -4, 
      35, 0.001, -1
    ), 
    hi = c(
      log10(40), 
      4, 13.46728, 2, 1, 
      1.5, 1, 4, 4, 
      -1.3, 
      49, 89.990, 1
    )
  )
  

  ## do rejection sampling to get mags 
  if(!file.exists(phot_stub)){
    # parm_samples <- runif(
    #   n   = Nsample * length(intervals$lo),
    #   min = rep(intervals$lo, each = Nsample),
    #   max = rep(intervals$hi, each = Nsample)
    # ) |> matrix(nrow = Nsample, ncol = length(intervals$lo))
    
    Naccept = 0
    count = 1
    prospect_samples = list()
    while(Naccept < Nsample){
      if(count %% 1000 == 0){
        print(Naccept)
      }
      count = count + 1
      
      parm <- runif(
        n   = length(intervals$lo),
        min = rep(intervals$lo),
        max = rep(intervals$hi)
      )
      
      parm[7] = runif(n = 1, min = intervals$lo[7], max = parm[6]) ## tau_screen < tau_birth
      parm[logged] = 10^parm[logged]
      
      magemax = celestial::cosdistUniAgeAtz(z = parm[1], ref = '737')
      agemax = magemax * 1e9
      Zagemax = magemax
      
      mpeak_parm = runif(1, min = -2, max = magemax)
      parm[3] = mpeak_parm
      
      phot = ProSpectSED(
        z = parm[1],
        mSFR = parm[2],
        mpeak = parm[3],
        mperiod = parm[4],
        mskew = parm[5],
        tau_birth = parm[6],
        tau_screen = parm[7],
        alpha_SF_birth = parm[8],
        alpha_SF_screen = parm[9],
        Zfinal = parm[10],
        AGNlum = parm[11],
        AGNan = parm[12],
        AGNta = parm[13], 
        
        agemax = agemax,
        magemax = magemax,
        Zagemax = Zagemax,
        
        AGN = Fritz,
        speclib = BC03hr,
        Dale = Dale_NormTot,
        Dale_M2L_func = FALSE,
        emission = TRUE,
        filtout = filter_transforms,
        filters = filters_info$filter_names,
        returnall = TRUE,
        
        SFH = SFHfunc,
        SMstar = TRUE,
        massfunc = massfunc_snorm_trunc,
        Z = Zfunc_massmap_lin,
        IGMabsorb = 'Inoue14',
        
        mtrunc = 2,
        
        AGNrm = 60,
        AGNbe = -0.5,
        AGNal = 4.0,
        AGNct = 100,
        
        ref='737',
        H0 = 70,
        OmegaM = 0.3,
        OmegaL = 0.7, 
      )
      
      mag_synthetic = 
        -2.5*log10(
          mean(phot$Photom[grepl("F560W|F770W", filters_info$filters)])
        ) + 8.9
      
      p_accept <- p_obs(mag_synthetic)
      p_accept[is.na(p_accept)] <- 0  # outside observed range = reject

      u        <- runif(length(mag_synthetic))
      keep     <- u < p_accept ## so this allows you to get the relative PDF mass in a given magnitude bin/x
      
      if(keep){
        Naccept = Naccept + 1
        UV1500 = ProSpect::photom_lum(
          wave = phot$FinalLum$wave, 
          lum = phot$FinalLum$lum, 
          z = 0, 
          LumDist_Mpc = celestial::cosdistLumDist(z = phot$z, ref = '737'), 
          ref = '737', 
          filters = list(UV_tophat)
        )
        UV1500AGN = ProSpect::photom_lum(
          wave = phot$AGN$wave, 
          lum = phot$AGN$lum, 
          z = 0, 
          LumDist_Mpc = celestial::cosdistLumDist(z = phot$z, ref = '737'), 
          ref = '737', 
          filters = list(UV_tophat)
        )
        massfunc_args = c(
          mSFR = parm[2],
          mpeak = parm[3],
          mperiod = parm[4],
          mskew = parm[5],
          magemax = magemax,
          mtrunc = 2
        )
        SFR10 = do.call("integrate", c(list(f = massfunc_snorm_trunc, lower = 0, upper = 1e+07), massfunc_args))$value/1e+07
        
        df = list(
          'astro' = c(
            "z" = parm[1],
            "StellarMass" = unname(phot$Stars$SMstar['TotSMstar']),
            "SFR10" = SFR10,
            "SFR100" = phot$Stars$SFRburst,
            "AGNlum" = parm[11],
            "UV1500" = UV1500,
            "UV1500AGN" = UV1500AGN,
            "detmag" = mag_synthetic
          ),
          parm = parm,
          phot = phot$Photom
        )
        prospect_samples = c(
          prospect_samples, list(df)
        )
      }
    }
    saveRDS(prospect_samples, phot_stub)
  }else{
    message("Reading ProSpect samples")
    prospect_samples = readRDS(phot_stub)
  }
  
  synthetic_mags = c(sapply(
    prospect_samples,
    function(x)as.numeric(x$astro['detmag'])
  ))
  
  synthetic_R50 <- mapply(sample_R50_given_mag,
                          mag_val = synthetic_mags,
                          MoreArgs = list(kde2 = kde2, n = 1))
  
  if(!file.exists(fits_stub)){
    
    ## Make 2D histogram samples of mag and R50
    # kde <- kde2d(profound_cat$ProFound$mag, log10(profound_cat$ProFound$R50), n = 200)
    # cdf_matrix <- cumsum(kde2$z) |> matrix(nrow = nrow(kde2$z))
    # cdf_matrix <- cdf_matrix / max(cdf_matrix)  # normalise to [0,1]
    # u <- runif(Nsample)
    # cdf_vec <- as.vector(cdf_matrix)
    # samples <- lapply(u, function(ui) {
    #   idx <- which.min(abs(cdf_vec - ui))  # nearest cell
    #   ix  <- ((idx - 1) %% nrow(kde2$z)) + 1
    #   iy  <- ((idx - 1) %/% nrow(kde2$z)) + 1
    #   c(mag = kde2$x[ix], R50 = kde2$y[iy])
    # })
    # samples_df <- do.call(rbind, samples) |> as.data.frame()
    # dx <- diff(kde2$x)[1] / 2
    # dy <- diff(kde2$y)[1] / 2
    samples_df = list()
    samples_df$mag <- synthetic_mags
    samples_df$R50 <- synthetic_R50
    
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
        nser = rep(1, Nsample),
        ang = runif(Nsample, 0.01, 89.99),
        axrat = runif(Nsample, 0.2, 1.0),
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
    
    Rfits_write_image(
      data = fake_sources,
      filename = fits_stub
    )
    saveRDS(
      samples_df,
      mag_samples_stub
    )
  }else{
    message("Reading fake sources")
    fake_sources = Rfits_read_image(fits_stub)
    samples_df = readRDS(mag_samples_stub)
  }

  if(!file.exists(pro_stub)){
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
    
    saveRDS(
      pro_fake,
      file = pro_stub
    )
  }else{
    message("Reading fake profound")
    pro_fake = readRDS(pro_stub)
  }
  
  samples_coords = cbind(Rwcs::Rwcs_p2s(cbind(samples_df$xcen, samples_df$ycen), keyvalues = fake_sources$keyvalues))
  astro_cat$RA = samples_coords[,1]
  astro_cat$DEC = samples_coords[,2]
  
  match_idx = celestial::coordmatch(
    coordref = astro_cat[,c("RA", "DEC")],
    coordcompare = pro_fake$segstats[,c("RAmax", "Decmax")], 
    rad = 0.3, 
    radunit = "asec"
  )
  test_sep = match_idx$sep
  test_sep[is.na(test_sep)] = Inf
  idx = test_sep > 0.3
  match_idx$sep[idx] = NA
  trim_seps = rowMins(match_idx$sep, na.rm = TRUE)
  sum(!is.infinite(trim_seps))
  
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
  
  frame_sigma = unlist(unname(frame_info[, paste0(filters_info$filters, "_5sigma_depth_32asec")]))/5.0
  Nflux = length(filters_info$filters)
  
  df_flux_raw = do.call(
    rbind, 
    lapply(prospect_samples, function(x)x$phot)
  )
  df_flux_err = data.frame(matrix(
    frame_sigma, 
    ncol = Nflux,
    nrow = dim(df_flux_raw)[1], 
    byrow = TRUE
  ))
  
  df_flux_shuffle = matrix(rnorm(n = length(df_flux_raw), mean=as.matrix(df_flux_raw), sd=as.matrix(df_flux_err)), nrow(df_flux_raw))
  
  df_flux = data.frame(df_flux_shuffle)
  
  names(df_flux) = paste0(filters_info$filters, "_flux")
  names(df_flux_err) = paste0(filters_info$filters, "_flux_err")
  
  # A_scale = (2 * pro_fake$segstats$flux[match_idx$bestmatch$compareID] - df_flux$F560W_flux[match_idx$bestmatch$compareID]) / df_flux$F770W_flux[match_idx$bestmatch$compareID]
  # df_flux$F770W_flux[match_idx$bestmatch$compareID] = df_flux$F770W_flux[match_idx$bestmatch$compareID] * A_scale
  
  recovered_det_flux = rowMeans(cbind(df_flux$F560W_flux, df_flux$F770W_flux))
  
  flag_unrealistically_negative_fluxes = (df_flux + 4*df_flux_err) <= 0
  
  df_flux[flag_unrealistically_negative_fluxes] = -999

  ztrue = sapply(prospect_samples, function(x)as.numeric(x$astro['z']))
  
  astro_cat = data.frame(do.call(rbind, lapply(prospect_samples, function(x){x$astro[1:length(x$astro)]})))
  astro_cat$detected = rep(0, Nsample)
  astro_cat$detected[match_idx$bestmatch$refID] = 1L
  astro_cat$recoveredmag = rep(NA, Nsample)
  astro_cat$recoveredmag[match_idx$bestmatch$refID] = pro_fake$segstats$mag[match_idx$bestmatch$compareID]
  
  
  pro_fake_all = pro_fake$segstats
  pro_fake_recovered = pro_fake$segstats[match_idx$bestmatch$compareID, ]
  pro_fake_contaminants = pro_fake$segstats[!(1:dim(pro_fake$segstats)[1] %in% match_idx$bestmatch$compareID), ]
  
  ## worse case scenario - only fit the bands with the most coverage...
  bands_to_keep = paste0(c(
    "F435W", "F606W", "F775W", "F814W", "F850LP",
    "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W", 
    "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W"
  ), "_flux")
  bands_to_set_NA = names(df_flux) %notin% bands_to_keep 
  df_flux[,bands_to_set_NA] = -999
  
  flux_df = cbind(
    ID = 1:Nsample, df_flux, df_flux_err, zspec = ztrue
  )
  
  fwrite(
    flux_df, 
    paste0(main_dir, "/", fname, "/", fname, "_mock_eazy_photometry.csv")
  )
  fwrite(
    astro_cat, 
    paste0(main_dir, "/", fname, "/", fname, "_mock_cat.csv")
  )
  fwrite(
    pro_fake_all, 
    paste0(main_dir, "/", fname, "/", fname, "_mock_all_segstats.csv")
  )
  fwrite(
    pro_fake_recovered, 
    paste0(main_dir, "/", fname, "/", fname, "_mock_recovered_segstats.csv")
  )
  fwrite(
    pro_fake_contaminants, 
    paste0(main_dir, "/", fname, "/", fname, "_mock_contaminants_segstats.csv")
  )
  
  eazy_translate = read.table(paste0("~/Documents/Epochs/eazy/eazy_", fname, "_zphot.translate"))
  eazy_translate$V1 = str_remove_all(str_remove_all(eazy_translate$V1, "_app_0.32"), "_scaled")
  fwrite(
    eazy_translate, 
    paste0(main_dir, "/eazy_zphot_mock.translate"),
    sep = " ", 
    col.names = FALSE,
    row.names = FALSE
  )
  
  return(
    NULL
  )
}

compl = source_inject(det_image, profound_cat, pro_detect, Nsample, fname = 'smiles')

