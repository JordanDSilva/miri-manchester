library(ProSpect)
library(matrixStats)
library(foreach)
library(doParallel)
library(data.table)

ncores = 1 
max_slurm_jobs = 84 # should be the maximum of your job array
slurm_paralel_job_array = FALSE
do_lrd = FALSE

args <- commandArgs(trailingOnly = TRUE)
if (length(args[1]) == 0 | is.na(args[1])) {
  idx <- 1
} else if (grepl(":", args[1])) {
  bounds <- as.integer(strsplit(args[1], ":")[[1]])
  idx <- bounds[1]:bounds[2]
} else if (grepl(",", args[1])) {
  idx <- as.integer(strsplit(args[1], ",")[[1]])
}else {
  idx <- as.integer(args[1])
}

if (length(args[2]) == 0 | is.na(args[2])) {
  main_dir = "~/Documents/Epochs/"
} else {
  main_dir <- args[2]
}
if (length(args[3]) == 0 | is.na(args[3])) {
  fittype = "photoZ_bc03_snorm_fritz_norm_eazy"
} else {
  fittype <- args[3]
}

prospect_dir = paste0(
  main_dir, "/data/prospect/", fittype, "/" 
)
save_dir = paste0(
  main_dir, "/data/prospect_seds/", fittype, "/" 
)
dir.create(
  save_dir, recursive = TRUE
)

load(
  paste0(main_dir, "/data/prospect_data_for_fitting.rda")
)

index_parallel_job_array = function(INDEX, ncores, N, start_Nindex, end_index){
  # will find optimal lots of ncores to fit N galaxies starting at start_Nindex:... inside of a 1 to max_index's slurm JOB array
  # Example N = 200000, ncores = 7, max_index = 800
  # At the minimum I will need floor(20000/(800*7)) = 35 lots of 7 indices 
  # JOB = 1 indices will be 1:7, then 5601:5608, then 11201:11207 ... 35 lots of these < 200000
  # 800*7*35 = 196000 so there will be 4000 remainders 
  # floor(4000/7) = 571 so tack on 7 lots of remainders for JOB = 1-571, and the last 3 stragglers on to JOB=572
  # I could fit 1:100000 galaxies N=100000, then fit 100001:200000 but need to set start_Nindex=100001
  
  maxjobs = end_index * ncores
  how_many_groups = floor(N/maxjobs)
  N_remaining = N - maxjobs*how_many_groups
  
  a = INDEX*ncores - ncores 
  b = INDEX*ncores - 1
  
  if(INDEX == 1){ 
    a = 1 
    b = ncores
    c = 0
  }else{
    c = 1
  }
  job_vec_init = (a:b) + c
  
  job_vec = c()
  for(k in 1:how_many_groups){
    job_vec = c(
      job_vec, job_vec_init + (k-1)*maxjobs
    )
  }
  
  if(N_remaining > 0 ){
    remaining_indices = (maxjobs*how_many_groups + 1):N
    job_vec = as.numeric(na.omit( c(job_vec, remaining_indices[job_vec_init]) )) + (start_Nindex-1)
  }else{
    job_vec = as.numeric(na.omit( c(job_vec) )) + (start_Nindex-1)
  }
  
  return( job_vec ) 
}

message(
  paste0("Total number of objects = ", dim(phot_cat$PhotTot)[1])
)

if(do_lrd){
  lrd_ids = unlist(
    fread(
      paste0(main_dir, "/data/LRD_colour_selection_IDs.txt")
    )
  )
}

if(slurm_paralel_job_array){
  if(do_lrd){
    job_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = length(lrd_ids), 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
    fit_idx = lrd_ids[job_idx]
  }else{
    fit_idx = index_parallel_job_array(
      INDEX = idx, 
      ncores = ncores, 
      N = dim(df_flux)[1], 
      start_Nindex = 1,
      end_index = max_slurm_jobs
    )
  }
}else{
  if(do_lrd){
    fit_idx = lrd_ids[idx]
  }else{
    fit_idx = idx
  }
}

plot_spectrum_ribbon <- function(spectrum,
                                 key,
                                 xvar   = "restWave",
                                 col    = "black",
                                 alpha  = 0.3,
                                 add    = FALSE,
                                 ylim   = NULL,
                                 xlim   = NULL,
                                 xlab   = "Wavelength",
                                 ylab   = "Luminosity",
                                 ...) {
  
  # build the Q16/Q50/Q84/Q95 column names from the key
  q16_name <- paste0(key, "Q16")
  q50_name <- paste0(key, "Q50")
  q84_name <- paste0(key, "Q84")
  q95_name <- paste0(key, "Q95")   # used only for auto ylim
  
  x   <- spectrum[[xvar]]
  q16 <- spectrum[[q16_name]]
  q50 <- spectrum[[q50_name]]
  q84 <- spectrum[[q84_name]]
  
  if (is.null(x) || is.null(q16) || is.null(q50) || is.null(q84)) {
    stop(sprintf("Missing one or more required columns: %s, %s, %s, %s",
                 xvar, q16_name, q50_name, q84_name))
  }
  
  # keep only rows that are finite (drops NA/NaN/Inf) AND strictly positive
  # (needed both to filter bad data and because log-scale can't plot 0/negative anyway)
  keep <- is.finite(x)   & x   > 0 &
    is.finite(q16) & q16 > 0 & 
    is.finite(q50) & q50 > 0 & 
    is.finite(q84) & q84 > 0 
  
  x   <- x[keep]
  q16 <- q16[keep]
  q50 <- q50[keep]
  q84 <- q84[keep]
  
  if (length(x) == 0) {
    stop(sprintf("No valid (finite, positive) data remaining for key '%s' after filtering", key))
  }
  
  # default axis ranges
  if (is.null(xlim)) xlim <- range(x, na.rm = TRUE)
  if (is.null(ylim)) {
    q95 <- spectrum[[q95_name]]
    upper_ref <- if (!is.null(q95)) q95 else q84
    upper_ref <- upper_ref[is.finite(upper_ref) & upper_ref > 0]
    ylim <- pmax(1e-1, range(upper_ref, na.rm = TRUE))
  }
  
  # convert named colour to RGBA with transparency for the fill
  poly_col <- adjustcolor(col, alpha.f = alpha)
  
  if (!add) {
    plot(NA, xlim = xlim, ylim = ylim, log = "xy",
         xlab = xlab, ylab = ylab, ...)
  }
  
  # shaded Q16-Q84 ribbon as a polygon
  polygon(
    c(x, rev(x)),
    c(q16, rev(q84)),
    col    = poly_col,
    border = NA
  )
  
  # median line on top
  lines(x, q50, col = adjustcolor(col, 0.5), lwd = 2)
  
  invisible(list(x = x, q16 = q16, q50 = q50, q84 = q84))
}

plot_sfh_zh <- function(spectrum,
                        key,
                        xvar,
                        col     = "black",
                        alpha   = 0.3,
                        ylab    = NULL,
                        ylab2   = NULL,
                        xlab    = "Age (Gyr)",
                        add     = FALSE,
                        lwd     = 2,
                        ...) {
  
  # build column names from the key
  q16_name <- paste0(key, "Q16")
  q50_name <- paste0(key, "Q50")
  q84_name <- paste0(key, "Q84")
  
  x   <- spectrum[[xvar]]
  q16 <- spectrum[[q16_name]]
  q50 <- spectrum[[q50_name]]
  q84 <- spectrum[[q84_name]]
  
  if (is.null(x) || is.null(q16) || is.null(q50) || is.null(q84)) {
    stop(sprintf("Missing required column(s): %s, %s, %s, %s",
                 xvar, q16_name, q50_name, q84_name))
  }
  
  # drop NA/Inf rows across all series, keep everything aligned
  keep <- is.finite(x) & is.finite(q16) & is.finite(q50) & is.finite(q84)
  x   <- x[keep]
  q16 <- q16[keep]
  q50 <- q50[keep]
  q84 <- q84[keep]
  
  if (length(x) == 0) {
    stop(sprintf("No valid data remaining for key '%s' after filtering", key))
  }
  
  poly_col <- adjustcolor(col, alpha.f = alpha)
  
  if (!add) {
    # first call: draw primary series + ribbon on the standard left axis
    ylim <- range(c(q16, q84), na.rm = TRUE)
    
    plot(
      x, q50,
      type = 'l',
      lwd  = lwd,
      col  = col,
      xlab = xlab,
      ylab = ylab,
      ylim = ylim,
      ...
    )
    
    polygon(
      c(x, rev(x)),
      c(q16, rev(q84)),
      col    = poly_col,
      border = NA
    )
    lines(x, q50, col = col, lwd = lwd)   # redraw median on top of ribbon
    
  } else {
    # second call: overlay this series + ribbon on a new right-hand axis
    y_max <- max(q84, na.rm = TRUE)
    y_min <- min(0, min(q16, na.rm = TRUE))
    
    par(usr = c(
      par()$usr[1:2],
      y_min - y_max * 0.04,
      y_max * 1.04
    ))
    
    polygon(
      c(x, rev(x)),
      c(q16, rev(q84)),
      col    = poly_col,
      border = NA
    )
    lines(x, q50, col = col, lwd = lwd)
    axis(4, col = col, col.axis = col)
    
    if (!is.null(ylab2)) {
      mtext(ylab2, side = 4, line = 3, col = col)
    }
  }
  
  invisible(list(x = x, q16 = q16, q50 = q50, q84 = q84))
}

registerDoParallel(cores = ncores)
for(i in fit_idx){
  message(i)
  
  file_stub = paste0(
    prospect_dir, "galID_", i, ".rds"
  )
  
  save_stub = paste0(
    save_dir, "galID_", i, ".rds"
  )
  
  if(file.exists(file_stub) & !file.exists(save_stub)){

    out = readRDS(file_stub)
    
    post = out$post
    
    Data = out$bestfit$Data
    Data$speclib = BC03hr
    
    if(!is.null(out$bestfit$SEDout$AGN)){
      if(do_lrd){
        if(grepl("MBB", fittype)){
          Data$AGN = LRDBBSKIRTOR
        }else{
          stop("Unknown fittype")
        }
      }else{
        Data$AGN = Fritz
      }
    }else{
      Data$AGN = NULL
    }
    Data$Dale = Dale_NormTot
    
    lum_samples = foreach(j = 1:dim(post)[1]) %dopar% {
      
      post_sample = out$highout$LD_last$Posterior1[j, names(out$bestfit$parm)]
      out_sample = ProSpectSEDlike(
        post_sample,
        Data
      )
      
      if(is.null(Data$AGN)){
        AGN = rep(0,10)
      }else{
        AGN = out_sample$SEDout$AGN$lum
      }
      
      if(!is.null(Data$AGN)){
        photom_AGN_source = photom_lum(
          wave = out_sample$SEDout$AGN$wave,
          lum = out_sample$SEDout$AGN$lum,
          filters = filter_transforms, 
          outtype = "Jy", 
          z = 10^out_sample$parm[1], 
          ref = "737"
        )
        AGN_Dale_frac = sum( out_sample$SEDout$AGN$lum[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4] * c(0, diff(out_sample$SEDout$AGN$wave[out_sample$SEDout$AGN$wave >= 5e4 & out_sample$SEDout$AGN$wave <= 20e4])) ) / 
          sum( out_sample$SEDout$FinalLum$lum[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4] * c(0, diff(out_sample$SEDout$FinalLum$wave[out_sample$SEDout$FinalLum$wave >= 5e4 & out_sample$SEDout$FinalLum$wave <= 20e4])) )
      }else{
        photom_AGN_source = rep(0, length(filter_transforms))
        AGN_Dale_frac = 0
      }
      
      photom_all_source = photom_lum(
        wave = out_sample$SEDout$FinalLum$wave,
        lum = out_sample$SEDout$FinalLum$lum,
        filters = filter_transforms, 
        outtype = "Jy", 
        z = 10^out_sample$parm[1], 
        ref = "737"
      )
      
      photom_all = photom_all_source * Inoue14_IGM(wave = filters_info$cenwave, z = 10^out_sample$parm[1], Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
      photom_AGN = photom_AGN_source *  Inoue14_IGM(wave = filters_info$cenwave, z = 10^out_sample$parm[1], Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
      
      names(photom_AGN) = paste0(filters_info$filters, "_fit_AGN_flux")
      names(photom_all) = paste0(filters_info$filters, "_fit_flux")
      names(photom_AGN_source) = paste0(filters_info$filters, "_source_fit_AGN_flux")
      names(photom_all_source) = paste0(filters_info$filters, "_source_fit_flux")
      
      list(
        "rest" = out_sample$SEDout$FinalLum$lum,
        "StarsUnAtten" = out_sample$SEDout$StarsUnAtten$lum,
        "StarsAtten" = out_sample$SEDout$StarsAtten$lum,
        "DustEmit" = out_sample$SEDout$DustEmit$lum,
        "AGN" = AGN,
        "obs" = out_sample$SEDout$FinalFlux$flux, 
        "agevec" =  out_sample$SEDout$Stars$agevec/1e9,
        "sfhvec" = out_sample$SEDout$Stars$SFR,
        "Zvec" = out_sample$SEDout$Stars$Zvec,
        "photom" = c(
          photom_all,
          photom_AGN,
          photom_all_source,
          photom_AGN_source,
          "AGN_Dale_frac" = AGN_Dale_frac
        )
      )
    }
    
    rest_samples         = do.call(rbind, lapply(lum_samples, function(x) x$rest))
    StarsUnAtten_samples = do.call(rbind, lapply(lum_samples, function(x) x$StarsUnAtten))
    StarsAtten_samples   = do.call(rbind, lapply(lum_samples, function(x) x$StarsAtten))
    DustEmit_samples     = do.call(rbind, lapply(lum_samples, function(x) x$DustEmit))
    AGN_samples          = do.call(rbind, lapply(lum_samples, function(x) x$AGN))
    obs_samples          = do.call(rbind, lapply(lum_samples, function(x) x$obs))
    phot_samples         = do.call(rbind, lapply(lum_samples, function(x) x$photom))
    
    med_phot = colMedians(as.matrix(phot_samples))
    Q16_phot = colQuantiles(as.matrix(phot_samples), probs = 0.16)
    Q84_phot = colQuantiles(as.matrix(phot_samples), probs = 0.84)
    Q05_phot = colQuantiles(as.matrix(phot_samples), probs = 0.05)
    Q95_phot = colQuantiles(as.matrix(phot_samples), probs = 0.95)
    
    names(med_phot) = paste0(names(med_phot), "Q50")
    names(Q16_phot) = paste0(names(Q16_phot), "Q16")
    names(Q84_phot) = paste0(names(Q84_phot), "Q84")
    names(Q05_phot) = paste0(names(Q05_phot), "Q05")
    names(Q95_phot) = paste0(names(Q95_phot), "Q95")
    
    row_phot <- c(ID = as.integer(i), med_phot, Q16_phot, Q84_phot, Q05_phot, Q95_phot)
    
    max_age_len <- max(sapply(lum_samples, function(x) length(x$agevec)))
    agevec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
      v <- x$agevec
      length(v) <- max_age_len   # pads with NA automatically if v is shorter
      v
    }))
    sfhvec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
      v <- x$sfhvec
      length(v) <- max_age_len   # pads with NA automatically if v is shorter
      v
    }))
    Zvec_samples <- do.call(rbind, lapply(lum_samples, function(x) {
      v <- x$Zvec
      length(v) <- max_age_len   # pads with NA automatically if v is shorter
      v
    }))
    
    spectrum = list()
    
    # fitted photom
    spectrum$phot = row_phot
    spectrum$phot_post = phot_samples

    # rest frame total
    spectrum$restWave     = out$bestfit$SEDout$FinalLum$wave
    spectrum$lumQ50       = colMedians(as.matrix(rest_samples), na.rm = TRUE)
    spectrum$lumQ16       = colQuantiles(as.matrix(rest_samples), probs = 0.16, na.rm = TRUE)
    spectrum$lumQ84       = colQuantiles(as.matrix(rest_samples), probs = 0.84, na.rm = TRUE)
    spectrum$lumQ05       = colQuantiles(as.matrix(rest_samples), probs = 0.05, na.rm = TRUE)
    spectrum$lumQ95       = colQuantiles(as.matrix(rest_samples), probs = 0.95, na.rm = TRUE)

    # unattenuated stars
    spectrum$StarsUnAttenWave = out$bestfit$SEDout$StarsUnAtten$wave
    spectrum$StarsUnAttenQ50 = colMedians(as.matrix(StarsUnAtten_samples), na.rm = TRUE)
    spectrum$StarsUnAttenQ16 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.16, na.rm = TRUE)
    spectrum$StarsUnAttenQ84 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.84, na.rm = TRUE)
    spectrum$StarsUnAttenQ05 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.05, na.rm = TRUE)
    spectrum$StarsUnAttenQ95 = colQuantiles(as.matrix(StarsUnAtten_samples), probs = 0.95, na.rm = TRUE)

    # attenuated stars
    spectrum$StarsAttenWave = out$bestfit$SEDout$StarsAtten$wave
    spectrum$StarsAttenQ50   = colMedians(as.matrix(StarsAtten_samples), na.rm = TRUE)
    spectrum$StarsAttenQ16   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.16, na.rm = TRUE)
    spectrum$StarsAttenQ84   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.84, na.rm = TRUE)
    spectrum$StarsAttenQ05   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.05, na.rm = TRUE)
    spectrum$StarsAttenQ95   = colQuantiles(as.matrix(StarsAtten_samples), probs = 0.95, na.rm = TRUE)

    # dust emission
    spectrum$DustWave = out$bestfit$SEDout$DustEmit$wave
    spectrum$DustEmitQ50     = colMedians(as.matrix(DustEmit_samples), na.rm = TRUE)
    spectrum$DustEmitQ16     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.16, na.rm = TRUE)
    spectrum$DustEmitQ84     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.84, na.rm = TRUE)
    spectrum$DustEmitQ05     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.05, na.rm = TRUE)
    spectrum$DustEmitQ95     = colQuantiles(as.matrix(DustEmit_samples), probs = 0.95, na.rm = TRUE)

    # AGN
    if(!is.null(out$bestfit$SEDout$AGN)){
      spectrum$AGNWave = out$bestfit$SEDout$AGN$wave
      spectrum$AGNQ50          = colMedians(as.matrix(AGN_samples), na.rm = TRUE)
      spectrum$AGNQ16          = colQuantiles(as.matrix(AGN_samples), probs = 0.16, na.rm = TRUE)
      spectrum$AGNQ84          = colQuantiles(as.matrix(AGN_samples), probs = 0.84, na.rm = TRUE)
      spectrum$AGNQ05          = colQuantiles(as.matrix(AGN_samples), probs = 0.05, na.rm = TRUE)
      spectrum$AGNQ95          = colQuantiles(as.matrix(AGN_samples), probs = 0.95, na.rm = TRUE)
    }else{
      spectrum$AGNWave = 1:10
      spectrum$AGNQ50 = rep(1e-99, 10)
      spectrum$AGNQ16 = rep(1e-99, 10)
      spectrum$AGNQ84 = rep(1e-99, 10)
      spectrum$AGNQ05 = rep(1e-99, 10)
      spectrum$AGNQ95 = rep(1e-99, 10)
    }

    # observed frame flux
    spectrum$obswave         = out$bestfit$SEDout$FinalFlux$wave
    spectrum$fluxQ50         = colMedians(as.matrix(obs_samples), na.rm = TRUE)
    spectrum$fluxQ16         = colQuantiles(as.matrix(obs_samples), probs = 0.16, na.rm = TRUE)
    spectrum$fluxQ84         = colQuantiles(as.matrix(obs_samples), probs = 0.84, na.rm = TRUE)
    spectrum$fluxQ05         = colQuantiles(as.matrix(obs_samples), probs = 0.05, na.rm = TRUE)
    spectrum$fluxQ95         = colQuantiles(as.matrix(obs_samples), probs = 0.95, na.rm = TRUE)

    # observed frame flux
    spectrum$obswave         = out$bestfit$SEDout$FinalFlux$wave
    spectrum$fluxQ50         = colMedians(as.matrix(obs_samples), na.rm = TRUE)
    spectrum$fluxQ16         = colQuantiles(as.matrix(obs_samples), probs = 0.16, na.rm = TRUE)
    spectrum$fluxQ84         = colQuantiles(as.matrix(obs_samples), probs = 0.84, na.rm = TRUE)
    spectrum$fluxQ05         = colQuantiles(as.matrix(obs_samples), probs = 0.05, na.rm = TRUE)
    spectrum$fluxQ95         = colQuantiles(as.matrix(obs_samples), probs = 0.95, na.rm = TRUE)

    # observed frame flux
    spectrum$ageQ50         = colMedians(as.matrix(agevec_samples), na.rm = TRUE)
    spectrum$ageQ16         = colQuantiles(as.matrix(agevec_samples), probs = 0.16, na.rm = TRUE)
    spectrum$ageQ84         = colQuantiles(as.matrix(agevec_samples), probs = 0.84, na.rm = TRUE)
    spectrum$ageQ05         = colQuantiles(as.matrix(agevec_samples), probs = 0.05, na.rm = TRUE)
    spectrum$ageQ95         = colQuantiles(as.matrix(agevec_samples), probs = 0.95, na.rm = TRUE)

    spectrum$SFHQ50         = colMedians(as.matrix(sfhvec_samples), na.rm = TRUE)
    spectrum$SFHQ16         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.16, na.rm = TRUE)
    spectrum$SFHQ84         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.84, na.rm = TRUE)
    spectrum$SFHQ05         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.05, na.rm = TRUE)
    spectrum$SFHQ95         = colQuantiles(as.matrix(sfhvec_samples), probs = 0.95, na.rm = TRUE)

    spectrum$ZHQ50         = colMedians(as.matrix(Zvec_samples), na.rm = TRUE)
    spectrum$ZHQ16         = colQuantiles(as.matrix(Zvec_samples), probs = 0.16, na.rm = TRUE)
    spectrum$ZHQ84         = colQuantiles(as.matrix(Zvec_samples), probs = 0.84, na.rm = TRUE)
    spectrum$ZHQ05         = colQuantiles(as.matrix(Zvec_samples), probs = 0.05, na.rm = TRUE)
    spectrum$ZHQ95         = colQuantiles(as.matrix(Zvec_samples), probs = 0.95, na.rm = TRUE)
    
    pdf(
      file = paste0(dirname(save_stub), "/plot_", i, ".pdf"), 
      width = 10, 
      height = 12,
      pointsize = 16
    )
    par(
      mfrow = c(3,1), 
      mar = c(3.5, 3.5, 1.0, 3.5),   # bottom, left, top, right — trimmed top/right, kept left for y-labels
      oma = c(2.5, 2.5, 1.0, 2.5),   # shared bottom margin for a common x-axis label across all 3 panels
      cex = 1.1,                      # dialed back from 2.0 -- try this first
      cex.axis = 1.1,
      cex.lab = 1.1,
      lwd = 2,
      mgp = c(2.2, 0.7, 0)            # pulls axis title/labels/ticks closer to the plot
    )
    obs_xlo = min(Data$flux$cenwave) - 100
    obs_xhi = max(Data$flux$cenwave) + 1000
    obs_ylo = min(pmax(1e-10, Data$flux$flux - 3*Data$flux$fluxerr, na.rm = TRUE)) * 0.5
    obs_yup = max(Data$flux$flux + 3*Data$flux$fluxerr, na.rm = TRUE) * 1.5
    
    plot_spectrum_ribbon(
      spectrum = spectrum, key = "flux", xvar = "obswave", xlim = c(obs_xlo, obs_xhi), ylim = c(obs_ylo, obs_yup), xlab = "Wavelength [micron]", ylab = "Flux [Jy]"
    )
    for(index in 1:length(Data$flux$cenwave)){
      x = Data$flux$cenwave[index]
      
      q50 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ50")])
      q16 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ16")])
      q84 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_fluxQ84")])
      
      polygon(
        c(x-1e3, x-1e3, x+1e3, x+1e3),
        c(q16, q84, q84, q16),
        col    = adjustcolor(col = "forestgreen", 0.4),
        border = NA
      )
      points(
        x, q50, col = "forestgreen", pch = 16
      )
      
      q50 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ50")])
      q16 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ16")])
      q84 = as.numeric(row_phot[paste0(Data$flux$filter[index], "_fit_AGN_fluxQ84")])
      
      polygon(
        c(x-1e3, x-1e3, x+1e3, x+1e3),
        c(q16, q84, q84, q16),
        col    = adjustcolor(col = "magenta", 0.4),
        border = NA
      )
      points(
        x, q50, col = "magenta", pch = 16
      )
    }
    points(
      Data$flux$cenwave,
      Data$flux$flux,
      col = "black",
      pch = 1
    )
    magerr(
      Data$flux$cenwave,
      Data$flux$flux,
      ylo = Data$flux$fluxerr, 
      col = "black"
    )
    legend(
      x = "bottomright",
      c(
        paste0("z = ", round(out$astro$zQ50, 3)), 
        paste0("RA = ", round(out$galInfo$RAmax, 5)), 
        paste0("Dec = ", round(out$galInfo$Decmax, 5)), 
        "Data",
        "Fit",
        "Fit photom",
        "Fit AGN photom"
      ),
      pch = c(NA, NA, NA, 1, NA, 16, 16),
      lty = c(NA, NA, NA, NA, 1, NA, NA),
      lwd = c(NA, NA, NA, NA, 5, NA, NA),
      col = c(NA, NA, NA, "black", "grey", "forestgreen", "magenta"),
      cex = 0.5
    )
    
    plot_spectrum_ribbon(spectrum, key = "StarsUnAtten", xvar = "StarsUnAttenWave", col = "blue", xlim = c(10, 1e7))
    plot_spectrum_ribbon(spectrum, key = "StarsAtten", xvar = "StarsUnAttenWave", col = "cornflowerblue", add = TRUE)
    plot_spectrum_ribbon(spectrum, key = "AGN", xvar = "AGNWave", col = "purple", add = TRUE)
    plot_spectrum_ribbon(spectrum, key = "DustEmit", xvar = "DustWave", col = "brown", add = TRUE)
    plot_spectrum_ribbon(spectrum, key = "lum", xvar = "restWave", col = "black", add = TRUE, alpha = 0.2)
    
    plot_sfh_zh(spectrum, key = "SFH", xvar = "ageQ50", col = "blue", ylab = "SFR [Msun yr^-1]", xlab = "Lookback time since observation [Gyr]", log = "x")
    plot_sfh_zh(spectrum, key = "ZH", xvar = "ageQ50", col = "red", ylab2 = "Z [mass fraction]", add = TRUE, log = "x")
    dev.off()
    
    dir.create(
      dirname(save_stub), 
      recursive = TRUE
    )
    saveRDS(
      spectrum, 
      save_stub
    )
  }
}
stopImplicitCluster()
