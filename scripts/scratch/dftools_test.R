library(celestial)
library(ProSpect)
library(Rfits)
library(data.table)
library(magicaxis)
library(foreach)
library(Highlander)

main_dir = "~/Documents/Epochs/"
fname = "smiles"

filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))
filter_transforms = readRDS(paste0(main_dir, "data/throughputs/filt_funcs.rds"))

frame_info = data.frame(fread(paste0(main_dir, "data/depths/", fname, "_profound_frame_stats.csv")))

combined_cat = data.frame(fread(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_robust_catalogue.csv"
))
eazy_stats = data.frame(fread(
  "~/Documents/Epochs/data/eazy/smiles_agn_blue_sfhz/smiles_eazy_stats_agn_blue_sfhz.csv"
))

Vmax_data = data.frame(fread(paste0(main_dir, "/data/depths/", "zmax_5sigma_data.csv")))

# lbt_bins = c(seq(8.0, 13.0, 0.5), 13.2)
# zbins = approx(
#   x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
#   y = seq(0,40,1),
#   xout = lbt_bins
# )$y
# # zbins = seq(0.0, 12.0, 0.5)
# zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2
# 
# Mstar_bins = seq(6.0, 12.0, 0.5)

zbins = c(6.5, 7.50)
kk = 1
zidx = which(combined_cat$zQ50 >= zbins[kk] & combined_cat$zQ50 < zbins[kk + 1])
Vmax = (cosdistCoVol(z = Vmax_data$zmax[zidx], ref = "737")
 - cosdistCoVol(z = zbins[kk], ref = "737")) * 1e9 * (frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2)/(4*pi)

message("Compute Veff from integral of dz/[dVc dOmega]")
OmegaM = 0.30
OmegaL = 0.70
H0 = 70
Ez = function(z){
  sqrt(
    OmegaM * (1+z)^3 + OmegaL
  )
}
DH = 3000 * (H0/100)^-1
dVc = function(z){
  DH * (1+z)^2 * (cosdistAngDist(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = 0))^2 / Ez(z)
} ## differntial comoving volume element per dz and dOmega[steradian]

interp_func_load = readRDS("~/Desktop/zooweemama.rds")
interp_func = function(z_in, M_in){
  compl = interp_func_load$fn(
    z_in,
    M_in,
    zmids = interp_func_load$zmids,
    Mstar_mids = interp_func_load$Mstar_mids,
    mat = interp_func_load$compl_mat
  )
  
  cont = interp_func_load$fn(
    z_in,
    M_in,
    zmids = interp_func_load$zmids,
    Mstar_mids = interp_func_load$Mstar_mids,
    mat = interp_func_load$cont_mat
  )
  
  return((1-cont)/compl)
}


single_schechter <- function(x, parm) {
  M0    <- parm[1]
  phi1  <- 10^parm[2]
  alpha <- parm[3]
  
  log10e <- log10(exp(1))
  
  phi <- log10e * exp(-10^(x - M0)) * 10^(x - M0) * phi1 * 10^(alpha * (x - M0))
  
  return(phi)
}

double_schechter <- function(x, parm) {
  M0    <- parm[1]
  phi1  <- 10^parm[2]
  phi2  <- 10^parm[3]
  alpha <- parm[4]
  beta  <- parm[5]
  
  log10e <- log10(exp(1))
  
  phi <- log10e * exp(-10^(x - M0)) * 10^(x - M0) * (phi1 * 10^(alpha * (x - M0)) + phi2 * 10^(beta * (x - M0)))
  
  return(phi)
}
schechter_parms_z <- function(z) {
  
  # C0 + z*C1 + z^2*C2 for each parameter
  M0    <- 10.963 + 0.211*z - 0.094*z^2
  alpha <- -1.627 + 0.428*z - 0.154*z^2
  beta  <- -0.539 + 0.053*z - 0.017*z^2
  phi1  <- -3.212 + 0.059*z - 0.105*z^2  # log10
  phi2  <- -2.592 - 0.515*z + 0.045*z^2  # log10
  
  return(c(M0, phi1, phi2, alpha, beta))
}
double_schechter_z <- function(x, z) {
  parm <- schechter_parms_z(z)
  double_schechter(x, parm)
}

Veff <- function(x) {
  sapply(x, function(xi) {
    
    fuck_this_shit = function(z){
      ifelse(is.na(interp_func(z,xi)), 0.0, interp_func(z,xi)) * dVc(z) * frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2
    }
    # fuck_this_shit = function(z){
    #   1.0 * dVc(z) * frame_info$F444W_effarea_amin2 * 60^-2 * (pi/180)^2
    # }
    pracma::integral(
      fun = fuck_this_shit, xmin = zbins[kk], xmax = zbins[kk+1]
    )
  })
}

magplot(
  Mstar_bins,
  Veff(Mstar_bins)
)
magplot(
  Mstar_bins,
  interp_func(zmids[1],Mstar_bins)
)
Mlim = approx(x = interp_func(zmids[kk], seq(6, 12, 0.01)), y = seq(6, 12, 0.01), xout = 0.8)$y
# Mlim = 1.0
Mlim = interp_func_load$fn(
  zmids[kk], seq(6, 12, 0.01),
  zmids = interp_func_load$zmids,
  Mstar_mids = interp_func_load$Mstar_mids,
  mat = interp_func_load$cont_mat
)
obs.selection = function(x){
  ifelse(
    x < seq(6, 12, 0.01)[Mlim <= 0.6][1] | x > 11.5,
    0,
    1
  )
}

mstar_zbin = log10(combined_cat$StellarMassQ50[zidx])
mstar_err_zbin = 0.5*(log10(combined_cat$StellarMassQ84[zidx]) - log10(combined_cat$StellarMassQ16[zidx]))
# [mstar_zbin > Mlim]
df = dftools::dffit(
  x = mstar_zbin, 
  x.err = mstar_err_zbin,
  selection = Vmax, 
  # obs.selection = obs.selection,
  gdf = single_schechter,
  dx = 0.01, 
  # xmin = 6,
  # xmax = 12,
  n.bootstrap = 100,
  n.iterations = 1000,
  p.initial = c(10, -4, -0.5), 
  method = "Nelder-Mead", 
  prior = function(parm){
    logprior = sum(
      c(
        # dnorm(parm[1], mean = 10.5, sd = 0.3, log = TRUE),
        0
      ),
      na.rm = TRUE
    )
    return(logprior)
  }
)
dftools::dfplot(df, xlim = c(5, 12), ylim = c(1e-7, 1e-1))
df$fit$p.best

# points(
#   c(8.125, 8.875, 9.625, 10.375, 11.125),
#   c(49.65, 7.37, 1.32, 0.64, 0.04)*1e-4*0.75*0.5,
#   pch = 8, col = "black"
# )
# lines(
#   seq(6, 12, 0.1),
#   single_schechter(
#     seq(6, 12, 0.1),
#     parm = c(11.89, -6.33, -1.97)
#   ), 
#   col = "red"
# )

points(
  c(8.125, 8.875, 9.625),
  c(16.06, 1.64, 0.17)*1e-4*0.75*0.5,
  pch = 8, col = "black"
)
lines(
  seq(6, 12, 0.1),
  single_schechter(
    seq(6, 12, 0.1),
    parm = c(10.82, -6.31, -2.12)
  ), 
  col = "red"
)

