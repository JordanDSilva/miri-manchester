library(data.table)
library(magicaxis)
library(ProSpect)
library(celestial)
library(Rfits)
library(fields)

main_dir = "~/Documents/Epochs/"
fname = "smiles"

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))
phot_cat = Rfits_read(
  paste0(main_dir, "/data/catalogues/", fname, "_profound_multiband_cat.fits")
)
prospect_cat = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_combined_catalogue.csv"
))

## magnitude completeness
detection_compl = readRDS(
  "~/Documents/Epochs/source_inject/smiles/smiles_mock_1detection_samples.rds"
)
selection_compl = readRDS(
  "~/Documents/Epochs/source_inject/smiles_uniform_detect/smiles_uniform_detect_selection_completeness.rds"
)

lbt_bins = c(seq(8.0, 13.0, 0.5), 13.1)
zbins = approx(
  x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
  y = seq(0,40,1),
  xout = lbt_bins
)$y
# zbins = seq(1.0, 12.0, 0.5)
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

MUV_bins = seq(-30, -10, 0.5)
MUV_mids = MUV_bins[1:(length(MUV_bins)-1)] + diff(MUV_bins)/2.0

MAB_bins = seq(10, 40, 1.0)
MAB_mids = MAB_bins[1:(length(MAB_bins)-1)] + diff(MAB_bins)/2.0

kk = 8

zidx = which(prospect_cat$zQ50 >= zbins[kk] & prospect_cat$zQ50 < zbins[kk + 1] & phot_cat$ProFound$flag == 0)

Vmax_mag = function(parm, Data){
  
  z=Data$z
  wave=Data$wave
  lum=Data$lum
  ApMag_lim=Data$ApMag_lim
  filter=Data$filter
  
  photom_mag_values = photom_lum(
    wave = wave, 
    lum = lum, 
    z = z+parm, 
    outtype = 'magAB', 
    filters = filter,
    ref = "737"
  ) 
  
  photom_flux_values = 10^( -0.4 * (photom_mag_values - 8.9) ) * Inoue14_IGM(wave = c(56076.384, 76128.766), z = z+parm, Inoue14_LAFcoef = Inoue14_LAFcoef, Inoue14_DLAcoef = Inoue14_DLAcoef)
  
  mean_flux = mean(photom_flux_values, na.rm = TRUE)
  mean_mag = -2.5*log10(mean_flux) + 8.9
  
  return(mean_mag)
}

zmax = foreach(i = zidx, .combine = c) %do% {
  message(i)
  
  if(prospect_cat$FitType[i] == "StellarAGN"){
    sed_stub = paste0(
      "~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
    file_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
  }else{
    sed_stub = paste0(
      "~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
    file_stub = paste0(
      "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth/galID_", i, ".rds"
    )
  }
  
  if(file.exists(sed_stub)){
    fit_sed = readRDS(
      sed_stub
    )
    
    fit = readRDS(
      file_stub
    )
    
    det_limit = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9
    Data = list(
      wave = fit_sed$restwave,
      lum = fit_sed$lumQ50,
      z = prospect_cat$zQ50[i], 
      ApMag_lim = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9,
      filter = list(filter_transforms$JWST_MIRI_F560W, filter_transforms$JWST_MIRI_F770W)
    )
    
    parm_zgrid = seq(
      0.0, zbins[kk+1]-prospect_cat$zQ50[i], 0.001
    )
    
    redshift_mags = sapply(
      parm_zgrid, 
      function(z){
        Vmax_mag(z, Data)
      }
    )
    
    zmax = parm_zgrid[redshift_mags >= det_limit][1] + prospect_cat$zQ50[i]
    if(is.na(zmax)){
      zmax = zbins[kk+1]
    }
    return(zmax)
  }else{
    return(0)
  }
}

zmax_zbin = pmin(zbins[kk+1], zmax)
Vmax_zbin = frame_info$detect_effarea_amin2 * 60^-2 * (180/pi)^-2 / (4*pi) * 1e9 * (cosdistCoVol(z = zmax_zbin, ref = "737") - cosdistCoVol(z = zbins[kk], ref = "737"))

detection_compl_zbin = detection_compl$compl_fn(phot_cat$ProFound$mag[zidx])
selection_compl_zbin = selection_compl$fns$compl_fn(prospect_cat$zQ50[zidx], prospect_cat$UV1500Q50[zidx], selection_compl$compl_interp)

Ndim = c(length(MUV_mids))
temp = array(data = rep(0, prod(Ndim)), dim = Ndim)
for(j in seq_along(MAB_mids)){
    idx = prospect_cat$UV1500Q50[zidx] >= MUV_bins[j] & prospect_cat$UV1500Q50[zidx] < MUV_bins[j+1]
  # idx = phot_cat$ProFound$mag[zidx] >= MAB_bins[j] & phot_cat$ProFound$mag[zidx] < MAB_bins[j+1]
  # temp[j] = sum(1 / (Vmax_zbin[idx] * selection_compl_zbin[idx]), na.rm = TRUE)
    temp[j] = sum(1 / (Vmax_zbin[idx] * selection_compl$compl[kk,j]), na.rm = TRUE)
}

magplot(
  MUV_mids,
  temp, 
  log = 'y', 
  ylim = c(1e-5, 0.1),
  xlim = c(-22, -12)
)

magplot(
  MUV_mids, 
  selection_compl$compl[3,]
)


compute_mass_completeness <- function(log_mstar, mag, redshift, mag_limit, z_bins) {
  
  # step 1: select faintest 20%
  faint_idx <- mag > quantile(mag, 0.80)
  
  log_mstar_f <- log_mstar[faint_idx]
  mag_f       <- mag[faint_idx]
  z_f         <- redshift[faint_idx]
  
  # step 2: scale mass to limiting magnitude
  log_M_lim <- log_mstar_f - 0.4 * (mag_limit - mag_f)
  
  # step 3: 90th percentile in redshift bins
  z_bin <- cut(z_f, breaks = z_bins)
  M_lim_z <- tapply(log_M_lim, z_bin, quantile, probs = 0.80, na.rm = TRUE)
  
  return(M_lim_z)
}

mass_completeness = compute_mass_completeness(
  log_mstar = prospect_cat$UV1500Q50,
  mag = phot_cat$ProFound$mag,
  redshift = prospect_cat$zQ50,
  mag_limit = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9,
  z_bins = zbins
)
