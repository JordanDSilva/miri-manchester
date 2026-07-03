library(data.table)
library(magicaxis)
library(ProSpect)
library(celestial)

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
compl_repeats = readRDS(
  "~/Documents/Epochs/source_inject/smiles/smiles_mock_10detection_samples.rds"
)

logistic_completeness <- function(mag, parm) {
  phi_max   <- parm[1]  # maximum completeness (e.g. 1)
  mag_break <- parm[2]  # magnitude at 50% of maximum
  k         <- parm[3]  # steepness of the drop (larger = sharper)
  
  phi <- phi_max / (1 + exp(k * (mag - mag_break)))
  
  return(phi)
}

mag_true = unlist(lapply(1:10, function(i){compl_repeats[[i]]$samples_df$mag}))
mag_detect = unlist(lapply(1:10, function(i){compl_repeats[[i]]$matched_segstats$mag}))
mag_detect_true = unlist(lapply(1:10, function(i){compl_repeats[[i]]$samples_df$mag[compl_repeats[[i]]$match_idx$bestmatch$refID]}))

magbins = seq(10, 40, 1.0)

hh_true = maghist(
  mag_true,
  breaks = magbins,
  plot = FALSE
)
hh_detect = maghist(
  mag_detect,
  breaks = magbins,
  plot = FALSE
)
hh_detect_true = maghist(
  mag_detect_true,
  breaks = magbins,
  plot = FALSE
)

magnitude_completeness_fn_opt = function(parm, Data){
  LL = sum(abs(
    (Data$yy - logistic_completeness(Data$xx, parm))
  ), na.rm = TRUE)
  return(LL)
}
magnitude_completeness_opt = optim(
  par = c(0.8, 25, 1),
  fn = magnitude_completeness_fn_opt, 
  Data = list(
    "yy" = hh_detect$counts/hh_true$counts,
    "xx" = magbins
  )
)

magplot(
  hh_true$mids,
  hh_detect$counts/hh_true$counts, 
  pch = 16
)
points(
  hh_true$mids,
  (hh_detect$counts-hh_detect_true$counts)/hh_true$counts, 
  pch = 15
)
magerr(
  hh_true$mids,
  hh_detect$counts/hh_true$counts,
  ylo = sqrt( hh_detect$counts/hh_true$counts^2 * ((sqrt(hh_detect$counts)/hh_detect$counts)^2 + (sqrt(hh_true$counts)/hh_true$counts)^2) )
)
lines(
  magbins,
  logistic_completeness(
    magbins,
    parm = magnitude_completeness_opt$par
  )
)

foo = readRDS(
  "~/Desktop/mstar_sfr_compl_test.rds"
)

find_nearest_cell <- function(mstar_obs, sfr_obs, mstar_grid, sfr_grid) {
  m_clamped <- pmax(pmin(mstar_obs, max(mstar_grid)), min(mstar_grid))
  s_clamped <- pmax(pmin(sfr_obs,   max(sfr_grid)),   min(sfr_grid))
  
  i_mstar <- which.min(abs(mstar_grid - m_clamped))
  i_sfr   <- which.min(abs(sfr_grid   - s_clamped))
  return(c(i_mstar, i_sfr))
}

lbt_bins = c(
  c(seq(8.0, 13.0, 0.5), 13.1)
  # seq(cosdistTravelTime(z = 0, ref = "737"), cosdistTravelTime(z = 2, ref = "737"), 1.0),
  # seq(cosdistTravelTime(z = 2.4, ref = "737"), cosdistTravelTime(z = 40, ref = "737"), 0.5)
)
zbins = approx(
  x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
  y = seq(0,40,1),
  xout = lbt_bins
)$y
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2
print(zmids)

MUV_bins = seq(-30, -10, 0.5)
MUV_mids = MUV_bins[1:(length(MUV_bins)-1)] + diff(MUV_bins)/2.0

SFR_bins = seq(-5, 5, 0.5)
SFR_mids = SFR_bins[1:(length(SFR_bins)-1)] + diff(SFR_bins)/2.0

Mstar_bins = seq(6, 15, 0.5)
Mstar_mids = Mstar_bins[1:(length(Mstar_bins)-1)] + diff(Mstar_bins)/2.0

kk = 2
z_vol = frame_info$detect_effarea_amin2 * 60^-2 * (180/pi)^-2 / (4*pi) * 1e9 * (cosdistCoVol(z = zbins[kk+1], ref = "737") - cosdistCoVol(z = zbins[kk], ref = "737"))

zidx = which(prospect_cat$zQ50 >= zbins[kk] & prospect_cat$zQ50 < zbins[kk+1] & prospect_cat$selected == 1L & phot_cat$ProFound$flag == 0L)
zmax_zbin = pmin(zmax[zidx], zbins[kk+1])
# zmax_zbin = zmax[zidx]
Vmax_zbin = frame_info$detect_effarea_amin2 * 60^-2 * (180/pi)^-2 / (4*pi) * 1e9 * (cosdistCoVol(z = zmax_zbin, ref = "737") - cosdistCoVol(z = zbins[kk], ref = "737"))

# # get completeness for each object
# compl_values <- mapply(function(m, s) {
#   
#   idx <- find_nearest_cell(m, s, Mstar_mids, SFR_mids)
#   compl = foo$compl[kk, idx[1], idx[2]]
#   if(is.na(compl)){
#     compl = NA
#   }
#   
#   return(pmin(1,compl))
# }, log10(prospect_cat$StellarMassQ50[zidx]), log10(prospect_cat$SFR10Q50[zidx]))
# get completeness for each object
compl_values <- mapply(function(m, s) {
  
  idx <- find_nearest_cell(m, s, Mstar_mids, MUV_mids)
  compl = foo$compl[kk, idx[1], idx[2]]
  if(is.na(compl)){
    compl = NA
  }
  
  return(pmin(1,compl))
}, log10(prospect_cat$StellarMassQ50[zidx]), prospect_cat$UV1500Q50[zidx])

magnitude_completeness_values = logistic_completeness(
  phot_cat$ProFound$mag[zidx],
  parm = magnitude_completeness_opt$par
)


Ndim = c(length(Mstar_mids), length(MUV_mids))
temp = array(data = rep(0, prod(Ndim)), dim = Ndim)
for(i in seq_along(Mstar_mids)){
  for(j in seq_along(MUV_mids)){
    # idx = log10(prospect_cat$StellarMassQ50[zidx]) >= Mstar_bins[i] & log10(prospect_cat$StellarMassQ50[zidx]) < Mstar_bins[i+1] & log10(prospect_cat$SFR10Q50[zidx]) >= SFR_bins[j] & log10(prospect_cat$SFR10Q50[zidx]) < SFR_bins[j+1]
    idx = log10(prospect_cat$StellarMassQ50[zidx]) >= Mstar_bins[i] & log10(prospect_cat$StellarMassQ50[zidx]) < Mstar_bins[i+1] & prospect_cat$UV1500Q50[zidx] >= MUV_bins[j] & prospect_cat$UV1500Q50[zidx] < MUV_bins[j+1]
    compl_val = foo$compl[kk, i, j]
    cont_val = foo$cont[kk, i, j]
    temp[i,j] = sum(1 / (Vmax_zbin[idx] * foo$compl[kk, i, j] * magnitude_completeness_values[idx]), na.rm = TRUE)
    
    # if(!is.na(cont_val) & !is.na(compl_val)){
    #   if(cont_val >= 0.5 | compl_val <= 0.5){
    #     temp[i,j] = NA
    #   }else{
    #     temp[i,j] = sum(1 / (Vmax_zbin[idx] * foo$compl[kk, i, j] * magnitude_completeness_values[idx]), na.rm = TRUE)
    #   }
    # }
  }
}

double_schechter <- function(mstar, parm) {
  log_mstar_star <- parm[1]  # log10 of characteristic mass M*
  alpha1         <- parm[2]  # faint-end slope 1
  alpha2         <- parm[3]  # faint-end slope 2
  phi1_star      <- parm[4]  # normalisation 1
  phi2_star      <- parm[5]  # normalisation 2
  
  x <- 10^(mstar - log_mstar_star)
  
  phi <- log(10) * exp(-x) * x * (phi1_star * x^alpha1 + phi2_star * x^alpha2)
  
  return(phi)
}

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
  log_mstar = log10(prospect_cat$StellarMassQ50),
  mag = phot_cat$ProFound$mag,
  redshift = prospect_cat$zQ50,
  mag_limit = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9,
  z_bins = zbins
)

sfr_completeness = compute_mass_completeness(
  log_mstar = log10(prospect_cat$SFR10Q50),
  mag = phot_cat$ProFound$mag,
  redshift = prospect_cat$zQ50,
  mag_limit = -2.5*log10(frame_info$detect_5sigma_depth_32asec)+8.9,
  z_bins = zbins
)

magplot(
  log10(prospect_cat$StellarMassQ50[zidx]),
  log10(prospect_cat$SFR10Q50[zidx]), 
  xlim = c(8, 11),
  ylim = c(-3, 3)
)
ttemp = temp 
ttemp[is.na(ttemp)] = 0
contour(
  ttemp, 
  x = Mstar_mids,
  y = SFR_mids,
  add = TRUE, 
  levels = 10^seq(-6, -2, 0.5)
  # nlevels = 100
)
abline(v = mass_completeness[kk])
abline(h = sfr_completeness[kk])
legend(
  x = "topleft",
  legend = c(
    paste0("z = ", round(zmids[kk], 3))
  )
)

smf = rowSums(temp, na.rm = TRUE)/(diff(Mstar_bins))
magplot(
  Mstar_mids, smf, 
  log = 'y', 
  xlim = c(8.5, 11.5), 
  ylim = c(1e-5, 1e-1), 
  xlab = "Stellar Mass",
  ylab = "Phi(Mstar)"
)
abline(
  v = mass_completeness[kk], 
  lty = 2
)
legend(
  x = "topright", 
  legend = c(
    paste0("z = ", round(zmids[kk], 3))
  )
)
sfrf = colSums(temp, na.rm = TRUE)/(diff(SFR_bins))
magplot(
  SFR_mids, sfrf, 
  log = 'y', 
  xlim = c(-5.5, 2.5), 
  ylim = c(1e-5, 1e-1), 
  xlab = "Stellar Mass",
  ylab = "Phi(Mstar)"
)
abline(
  v = sfr_completeness[kk], 
  lty = 2
)


uvlf = colSums(temp, na.rm = TRUE)/(diff(MUV_bins))
magplot(
  MUV_mids, uvlf, 
  log = 'y', 
  xlim = c(-21.5, -10.5), 
  ylim = c(1e-5, 1e-1), 
  xlab = "Stellar Mass",
  ylab = "Phi(Mstar)"
)
abline(
  v = sfr_completeness[kk], 
  lty = 2
)

Ndim = c(length(Mstar_mids), length(MUV_mids))
temp = rep(0, length(MUV_mids))
for(j in seq_along(MUV_mids)){
  idx = prospect_cat$UV1500Q50[zidx] >= MUV_bins[j] & prospect_cat$UV1500Q50[zidx] < MUV_bins[j+1]
  compl_val = foo$completeness[kk, j]
  temp[j] = sum(1 / (Vmax_zbin[idx] * magnitude_completeness_values * compl_val), na.rm = TRUE)
}
magplot(
  MUV_mids, temp, pch = 16, log = 'y', ylim = c(1e-8, 1e-1), xlim = c(-12, -25)
)
