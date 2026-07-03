library(Rfits)
library(magicaxis)
library(ProSpect)
library(celestial)
library(scales)
library(dftools)
library(foreach)
library(hyper.fit)
library(matrixStats)

epochs = Rfits_read_table("~/Documents/Epochs/data/EPOCHS_MASTER_public_final.fits")

epochs_As = sum(
  c(
    57.32, 3.90, 12.3, 4.00, 66.40, 6.08, 4.31, 9.76, 1.28, 4.03, 22.98
  )
)
epochs_area = epochs_As * (pi/180)^2 * 60^-2

UV2SFR = function(MUV){
  magAB2CGS(MUV) * 1.15e-28 * (3.086e+19)^2 * 4*pi
}
SFR2UV = approxfun(
  UV2SFR(seq(-13, -40, -0.01)),
  seq(-13, -40, -0.01)
)
magplot(
  epochs$M_UV_pipes_zgauss,
  epochs$M_UV,
  xlim = c(-23, -12), 
  ylim = c(-23, -12)
  # z = epochs$zbest, 
  # xlim = c(1e-3, 1e3), 
  # ylim = c(1e-3, 1e3)
)
abline(0, 1)

png("~/Documents/Epochs/plots/rough-plots/LUV_SFR_compare.png", width = 8, height = 6, res = 240, units = "in")
par(mar = c(3.5, 3.5, 1.0, 1.0))
magplot(
  UV2SFR(epochs$M_UV_pipes_zfix),
  epochs$sfr_10Myr_pipes_zfix, 
  log = "xy", 
  pch = 16, 
  z = epochs$sfr_10Myr_pipes_zfix/epochs$sfr_100Myr_pipes_zfix, 
  xlim = c(0.1, 250), 
  ylim = c(0.1, 250), 
  xlab = "LUV * K_UV (1.15*10^-28)",
  ylab = "SFR_10Myr",
  title = "SFR10/SFR100",
  position = "centreright"
)
abline(0, 1)
dev.off()

schechter_function_UV = function(x, p){
  
  phi1 = 10^p[1]
  MUV_star = p[2]
  alpha = p[3]
  
  mm = 10^(0.4 * (MUV_star - x))
  
  ret = (0.4 * log(10)) * phi1 * mm^(alpha+1) * exp(-1 * mm)
  return(ret)
}
schechter_function = function(x, p){
  
  phi1 = 10^p[1]
  Mstar = p[2]
  alpha = p[3]
  
  mm = 10^((x - Mstar))
  
  ret = (log(10)) * phi1 * mm^(alpha+1) * exp(-1 * mm)
  return(ret)
}

adams24_z8UVLF = data.frame(
  "MUV" = c(-21.35, -20.85, -20.35, -19.85, -19.35, -18.60),
  "phi" = c(2.277, 9.974, 13.12, 28.64, 54.04, 69.35)*1e-5
)
adams24_z9UVLF = data.frame(
  "MUV" = c(-22.05, -21.55, -21.05, -20.55, -20.05, -19.55, -18.80),
  "phi" = c(0.628, 0.628, 1.257, 6.427, 10.76, 18.22, 42.45)*1e-5
)
adams24_z8UVLF_fit = schechter_function(
  x = seq(-13, -30, -0.1), 
  p = c(-3.19, -20, -1.74)
)

muv_bw = 0.5
muv_bins = seq(-25, -16, muv_bw)
muv_mids = muv_bins[-length(muv_bins)] + diff(muv_bins) / 2

zmin = 7.5
zmax = 8.5
veff = (cosdistCoVol(z = zmax, ref = "737") - cosdistCoVol(z = zmin, ref = "737")) * 1e9 * (epochs_area)/(4*pi)
zidx = epochs$zbest >= zmin & epochs$zbest < zmax & epochs$FIELDNAME == "JADES-Deep-GS"

hh_muv = maghist(
  x = epochs$M_UV[zidx], 
  breaks = muv_bins, 
  plot = FALSE
)
hh_muv_zfix = maghist(
  x = epochs$M_UV_pipes_zfix[zidx], 
  breaks = muv_bins, 
  plot = FALSE
)
hh_muv_zgauss = maghist(
  x = epochs$M_UV_pipes_zgauss[zidx], 
  breaks = muv_bins, 
  plot = FALSE
)

magplot(
  muv_mids,
  hh_muv$counts/(veff*muv_bw), 
  log = "y", 
  pch = 0,
  xlim = c(-22.0, -17.0), 
  ylim = c(1e-6, 5e-3), 
  xlab = "MUV [mag]",
  ylab = "Phi [Mpc^-3 dex^-1]"
)
points(
  muv_mids,
  hh_muv_zfix$counts/(veff*muv_bw), 
  pch = 1
)
points(
  muv_mids,
  hh_muv_zgauss$counts/(veff*muv_bw), 
  pch = 2
)
# points(
#   adams24_z8UVLF$MUV, 
#   adams24_z8UVLF$phi, 
#   pch = 2, col = "cornflowerblue"
# )
# points(
#   adams24_z9UVLF$MUV, 
#   adams24_z9UVLF$phi, 
#   pch = 2, col = "blue"
# )
# lines(
#   (seq(-13, -30, -0.1)), 
#   (adams24_z8UVLF_fit)
# )
legend(
  x = "topleft",
  pch = c(0,1,2),
  legend = c("MUV", "MUV zfix", "MUV zgauss")
)


hh_muv_from_sfr = maghist(
  x = SFR2UV(epochs$sfr_10Myr_pipes_zgauss[zidx]), 
  breaks = muv_bins, 
  plot = FALSE
)

# magplot(
#   muv_mids,
#   hh_muv$counts/(veff*muv_bw), 
#   log = "y", 
#   ylim = c(1e-6, 1e-2), 
#   xlim = c(-22.0, -17.0)
# )
# points(
#   muv_mids,
#   hh_muv_from_sfr$counts/(veff*muv_bw), 
#   pch = 16
# )

zmin = 6.5
zmax = 7.5
veff = (cosdistCoVol(z = zmax, ref = "737") - cosdistCoVol(z = zmin, ref = "737")) * 1e9 * (epochs_area)/(4*pi)
zidx = epochs$zbest >= zmin & epochs$zbest < zmax

x_sfr = seq(-3, 3, 0.01)
sfr_bw = 0.25
sfr_bins = seq(-3, 3, sfr_bw)
sfr_mids = sfr_bins[-length(sfr_bins)] + diff(sfr_bins) / 2

hh_sfr10 = maghist(
  x = log10(epochs$sfr_10Myr_pipes_zfix[zidx]), 
  breaks = sfr_bins, 
  plot = FALSE
)
hh_sfr100 = maghist(
  x = log10(epochs$sfr_100Myr_pipes_zfix[zidx]), 
  breaks = sfr_bins, 
  plot = FALSE
)
hh_sfr_from_muv = maghist(
  x = log10(UV2SFR(epochs$M_UV[zidx])), 
  breaks = sfr_bins, 
  plot = FALSE
)

png("~/Documents/Epochs/plots/rough-plots/SFRF_compare.png", width = 8, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0))
magplot(
  sfr_mids,
  hh_sfr10$counts/(veff*sfr_bw), 
  log = "y", 
  pch = 16,
  cex = 1,
  ylim = c(1e-6, 5e-3), 
  xlim = c(-1.0, 2.5),
  xlab = "SFR [Msol yr^-1]",
  ylab = "Phi [Mpc^-3 dex^-1]",
  main = paste0(zmin, "< z <", zmax)
)
points(
  sfr_mids,
  hh_sfr_from_muv$counts/(veff*sfr_bw), 
  pch = 1, 
  cex = 1.5
)
points(
  sfr_mids,
  hh_sfr100$counts/(veff*sfr_bw), 
  pch = 17, 
  cex = 1.0
)
sfrf_LL = function(p, Data){
  
  idx = Data$x > Data$xcut & Data$yerr > 0
  yy = schechter_function(
    x = Data$x[idx], 
    p = p
  )
  LL = dnorm(
    x = Data$y[idx],
    mean = yy, 
    sd = (Data$yerr[idx]+1e-99)^2, 
    log = TRUE
  )
  return(sum(LL, na.rm = TRUE))
}

opt_sfr = optim(
  par = c(-3.0, 1.0, -1.0),
  fn = sfrf_LL, 
  Data = list(
    "x" = sfr_mids,
    "y" = hh_sfr10$counts/(veff*sfr_bw),
    "yerr"= sqrt(hh_sfr$counts)/(veff*sfr_bw), 
    "xcut" = 0.4
  ), 
  control = list(fnscale = -1)
)
opt_sfr_from_MUV = optim(
  par = c(-3.0, 1.0, -1.0),
  fn = sfrf_LL, 
  Data = list(
    "x" = sfr_mids,
    "y" = hh_sfr_from_muv$counts/(veff*sfr_bw),
    "yerr"= sqrt(hh_sfr_from_muv$counts)/(veff*sfr_bw), 
    "xcut" = 0.4
  ), 
  control = list(fnscale = -1)
)
# lines(
#   x_sfr, 
#   schechter_function(
#     x_sfr, 
#     p = opt_sfr$par
#   )
# )
# lines(
#   x_sfr, 
#   schechter_function(
#     x_sfr, 
#     p = opt_sfr_from_MUV$par
#   ), 
#   lty = 2
# )
# lines(
#   log10(UV2SFR(seq(-13, -30, -0.1))),
#   (adams24_z8UVLF_fit), 
#   col = "blue"
# )
# points(
#   log10(UV2SFR(adams24_z8UVLF$MUV)),
#   adams24_z8UVLF$phi,
#   pch = 2, col = "cornflowerblue"
# )
legend(
  x = "topright",
  pch = c(1, 16, 17),
  legend = c("LUV->SFR", "SFR10", "SFR100")
)
dev.off() 

## Lower limits
zbins = c(6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 13.5, 15.5)
zmids = zbins[-length(zbins)] + diff(zbins) / 2

sfr_data = foreach(i = 1:(length(zbins)-1), .combine = rbind)%do%{
  zmin = zbins[i]
  zmax = zbins[i+1]
  veff = (cosdistCoVol(z = zmax, ref = "737") - cosdistCoVol(z = zmin, ref = "737")) * 1e9 * (epochs_area)/(4*pi)
  
  zidx = epochs$zbest >= zmin & epochs$zbest < zmax
  LUVdens = sum(magAB2CGS(epochs$M_UV_pipes_zgauss[zidx]) * (3.086e+19)^2 * 4*pi)/veff
  SFRUVdens = LUVdens * 1.15e-28
  SFR10dens = sum(epochs$sfr_10Myr_pipes_zgauss[zidx])/veff
  SFR100dens = sum(epochs$sfr_100Myr_pipes_zgauss[zidx])/veff
  data.frame(
    LUVdens, SFRUVdens, SFR10dens, SFR100dens, median(epochs$zbest[zidx])
  )
}
adams24_rhoUV = -0.36*sfr_data[,5] + 28.55
correction_factors = 10^adams24_rhoUV/sfr_data[,1]

png("~/Documents/Epochs/plots/rough-plots/csfh_lower_limit_SFR_compare.png", width = 8, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0))
magplot(
  sfr_data[,5], 
  10^adams24_rhoUV * 1.15e-28, 
  log = "y", 
  type = "l",
  lty = 2,
  ylim = c(5e-6, 5e-2), 
  xlim = c(6.5, 15.5),
  xlab = "Redshift",
  ylab = "CSFH [Msun yr^-1 Mpc^-3]"
)
points(
  sfr_data[,5],
  sfr_data[,2], 
  cex = 1.5,
  pch = 1
)
points(
  sfr_data[,5],
  sfr_data[,3], 
  pch = 16
)
points(
  sfr_data[,5],
  sfr_data[,4], 
  pch = 17
)
legend(
  x = "topright",
  pch = c(1, 16, 17, NA),
  lty = c(NA, NA, NA, 2),
  legend = c("LUV->SFR", "SFR10", "SFR100", "Adams+24 fit")
)
dev.off()
# points(
#   c(6.450, 8.231, 10.126, 12.583),
#   10^c(-1.707, -1.809, -2.451, -2.596), 
#   pch = 1, col = "blue"
# )

maghist(
  log10(epochs$sfr_10Myr_pipes_zgauss) - log10(magAB2CGS(epochs$M_UV[zidx]) * (3.086e+19)^2 * 4*pi), 
  breaks = 30, log = "y"
)
abline(
  v = log10(1.15e-28)
)

lookback_time = cosdistTravelTime(
  z = epochs$redshift_pipes_zgauss,
  ref = "737"
)


# magplot(0.5*((log10(epochs$sfr_10Myr_pipes_u1_zfix[zidx][declining_idx] + epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx])) + (log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx] - epochs$sfr_10Myr_pipes_l1_zfix[zidx][declining_idx]))),
# 
# 0.5*(epochs$sfr_10Myr_pipes_u1_zfix[zidx][declining_idx] + epochs$sfr_10Myr_pipes_l1_zfix[zidx][declining_idx])/ (log(10)*epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx])
# )
# abline(0,1)


## Mstar-SFR
zmin = 7.5
zmax = 8.5
veff = (cosdistCoVol(z = zmax, ref = "737") - cosdistCoVol(z = zmin, ref = "737")) * 1e9 * (epochs_area)/(4*pi)
message("Time between redshifts = ", round(1e3*(cosdistTravelTime(z = zmax, ref = "737") - cosdistTravelTime(z = zmin, ref = "737")),4), " Myr")
zidx = epochs$zbest >= zmin & epochs$zbest < zmax

Mstar_lim = 8

declining_idx = epochs$sfr_10Myr_pipes_zfix[zidx]/(epochs$sfr_100Myr_pipes_zfix[zidx]) < 1 & epochs$sfr_10Myr_pipes_zfix[zidx] > 1e-2 & epochs$sfr_100Myr_pipes_zfix[zidx] > 1e-2 & epochs$stellar_mass_pipes_zfix[zidx] > Mstar_lim

intermediate_idx = epochs$sfr_10Myr_pipes_zfix[zidx]/(epochs$sfr_100Myr_pipes_zfix[zidx]) >= 1 &
  epochs$sfr_10Myr_pipes_zfix[zidx]/(epochs$sfr_100Myr_pipes_zfix[zidx]) < 2 & epochs$sfr_10Myr_pipes_zfix[zidx] > 1e-2 & epochs$sfr_100Myr_pipes_zfix[zidx] > 1e-2 & epochs$stellar_mass_pipes_zfix[zidx] > Mstar_lim

rising_idx = epochs$sfr_10Myr_pipes_zfix[zidx]/(epochs$sfr_100Myr_pipes_zfix[zidx]) >= 2 & epochs$sfr_10Myr_pipes_zfix[zidx] > 1e-2 & epochs$sfr_100Myr_pipes_zfix[zidx] > 1e-2 & epochs$stellar_mass_pipes_zfix[zidx] > Mstar_lim

below_mstar_idx = epochs$stellar_mass_pipes_zfix[zidx] < Mstar_lim
# maghist(
#   epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx], log = "xy", 
# )
# maghist(
#   epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx], log = "xy", 
# )

sfh_idx = rep(1, length(zidx))
sfh_idx[intermediate_idx] = 2L
sfh_idx[declining_idx] = 3L
sfh_idx[below_mstar_idx] = 4L

# magplot(
#   epochs$stellar_mass_pipes_zfix[zidx], 
#   epochs$sfr_10Myr_pipes_zfix[zidx], 
#   z = rising_idx,
#   log = "y", 
#   ylim = c(0.005, 500), 
#   xlim = c(6.5, 11.0), 
#   zcol = c("red", "blue"), 
#   n = 0
# )
# magplot(
#   epochs$stellar_mass_pipes_zfix[zidx], 
#   epochs$sfr_10Myr_pipes_zfix[zidx], 
#   z = intermediate_idx,
#   log = "y", 
#   ylim = c(0.005, 500), 
#   xlim = c(6.5, 11.0), 
#   zcol = c("red", "blue"), 
#   n = 0
# )
# magplot(
#   epochs$stellar_mass_pipes_zfix[zidx], 
#   epochs$sfr_10Myr_pipes_zfix[zidx], 
#   z = declining_idx,
#   log = "y", 
#   ylim = c(0.005, 500), 
#   xlim = c(6.5, 11.0), 
#   zcol = c("red", "blue"), 
#   n = 0
# )
# magplot(
#   epochs$stellar_mass_pipes_zfix[zidx], 
#   epochs$sfr_10Myr_pipes_zfix[zidx], 
#   z = epochs$sfr_10Myr_pipes_zfix[zidx]/epochs$sfr_100Myr_pipes_zfix[zidx],
#   log = "y", 
#   ylim = c(0.005, 500), 
#   xlim = c(6.5, 11.0)
# )

png("~/Documents/Epochs/plots/rough-plots/SFR-UV-Mstar.png", width = 12, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0), mfrow = c(1,2))
magplot(
  epochs$stellar_mass_pipes_zfix[zidx], 
  UV2SFR(epochs$M_UV_pipes_zfix[zidx]), 
  z = sfh_idx,
  log = "y", 
  zcol = alpha(c("cornflowerblue", "forestgreen", "red", "darkgrey"), 0.5),
  dobar = FALSE, 
  pch = 16,
  xlim = c(6.5, 11.0), 
  ylim = c(0.1, 200),
  xlab = "log10(Mstar[Msol])",
  ylab = "LUV -> SFR [Msol yr^-1]", 
  main = paste(zmin, "< z < ", zmax)
)
legend(
  x = "topleft",
  col = c("cornflowerblue", "forestgreen", "red", "darkgrey"), 
  pch = 16,
  legend = c("SFR10/SFR100>2", "1<SFR10/SFR100<2", "SFR10/SFR100<1", "logMstar<8")
)

# png("~/Documents/Epochs/plots/rough-plots/SFR-Mstar_fit.png", width = 8, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0))
magplot(
  epochs$stellar_mass_pipes_zfix[zidx], 
  epochs$sfr_10Myr_pipes_zfix[zidx], 
  z = sfh_idx,
  log = "y", 
  zcol = alpha(c("cornflowerblue", "forestgreen", "red", "darkgrey"), 0.5),
  dobar = FALSE, 
  pch = 16,
  xlim = c(6.5, 11.0), 
  ylim = c(0.1, 200),
  xlab = "log10(Mstar[Msol])",
  ylab = "SFR10 [Msol yr^-1]", 
  main = paste(zmin, "< z < ", zmax)
)
magerr(
  epochs$stellar_mass_pipes_zfix[zidx][rising_idx],
  epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx],
  xlo = epochs$stellar_mass_pipes_l1_zfix[zidx][rising_idx],
  xhi = epochs$stellar_mass_pipes_u1_zfix[zidx][rising_idx],
  ylo = epochs$sfr_10Myr_pipes_l1_zfix[zidx][rising_idx],
  yhi = epochs$sfr_10Myr_pipes_u1_zfix[zidx][rising_idx], 
  col =alpha("cornflowerblue", 0.4)
)
magerr(
  epochs$stellar_mass_pipes_zfix[zidx][intermediate_idx],
  epochs$sfr_10Myr_pipes_zfix[zidx][intermediate_idx],
  xlo = epochs$stellar_mass_pipes_l1_zfix[zidx][intermediate_idx],
  xhi = epochs$stellar_mass_pipes_u1_zfix[zidx][intermediate_idx],
  ylo = epochs$sfr_10Myr_pipes_l1_zfix[zidx][intermediate_idx],
  yhi = epochs$sfr_10Myr_pipes_u1_zfix[zidx][intermediate_idx], 
  col =alpha("forestgreen", 0.4)
)
magerr(
  epochs$stellar_mass_pipes_zfix[zidx][declining_idx],
  epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx],
  xlo = epochs$stellar_mass_pipes_l1_zfix[zidx][declining_idx],
  xhi = epochs$stellar_mass_pipes_u1_zfix[zidx][declining_idx],
  ylo = epochs$sfr_10Myr_pipes_l1_zfix[zidx][declining_idx],
  yhi = epochs$sfr_10Myr_pipes_u1_zfix[zidx][declining_idx], 
  col =alpha("red", 0.4)
)
legend(
  x = "topleft",
  col = c("cornflowerblue", "forestgreen", "red", "darkgrey"), 
  pch = 16,
  legend = c("SFR10/SFR100>2", "1<SFR10/SFR100<2", "SFR10/SFR100<1", "logMstar<8")
)
abline(
  v = Mstar_lim, 
  lty = 2
)

hf_sfh_declining = hyper.fit(
  cbind(
    epochs$stellar_mass_pipes_zfix[zidx][declining_idx], 
    log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx])
  ), 
  vars = cbind(
    0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx][declining_idx] + epochs$stellar_mass_pipes_l1_zfix[zidx][declining_idx]),

    0.5*((log10(epochs$sfr_10Myr_pipes_u1_zfix[zidx][declining_idx] + epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx])) + (log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx] - epochs$sfr_10Myr_pipes_l1_zfix[zidx][declining_idx]+1e-99)))
  ), 
  algo.func = "LD", 
  algo.method = "CHARM",
  Specs=list(alpha.star = 0.44)
)
hf_sfh_rising = hyper.fit(
  cbind(
    epochs$stellar_mass_pipes_zfix[zidx][rising_idx], 
    log10(epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx])
  ),
  vars = cbind(
    0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx][rising_idx] + epochs$stellar_mass_pipes_l1_zfix[zidx][rising_idx]),

    0.5*((log10(epochs$sfr_10Myr_pipes_u1_zfix[zidx][rising_idx] + epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx])) + (log10(epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx]) - log10(epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx] - epochs$sfr_10Myr_pipes_l1_zfix[zidx][rising_idx]+1e-323)))
  ),
  algo.func = "LD",
  algo.method = "CHARM",
  Specs=list(alpha.star = 0.44)
)
# plot(hf_sfh_declining)
# plot(hf_sfh_rising)

sm_bins = seq(4.5, 12.5, 0.01)
for(i in 1:1000){
  lines(
    sm_bins, 10^(hf_sfh_declining$fit$Posterior1[i,1]*(sm_bins) + hf_sfh_declining$fit$Posterior1[i,2]), 
    col = alpha("darkred", 0.01)
  )
  
  lines(
    sm_bins, 10^(hf_sfh_rising$fit$Posterior1[i,1]*(sm_bins) + hf_sfh_rising$fit$Posterior1[i,2]), 
    col = alpha("blue", 0.01)
  )
}
# lines(
  #   sm_bins, 10^(hf_sfh_declining$parm[1]*(sm_bins) + hf_sfh_declining$parm[2]), col = alpha("darkred", 0.05)
  # )
# lines(
#   sm_bins, 10^(hf_sfh_rising$parm[1]*(sm_bins) + hf_sfh_rising$parm[2]), col = alpha("blue", 0.05)
# )
dev.off()



harvey24_smf = read.table(
  "~/Documents/Epochs/data/smf/FiducialBagpipesGSMF_SchechterParams.ecsv",
  header = TRUE,
  sep = "",          # any whitespace
  fill = TRUE,       # handles ragged rows like "9-2.26..."
  stringsAsFactors = FALSE, 
  skip = 24
)
z7_smf = schechter_function(sm_bins, p = c(harvey24_smf$logphistar_50[1], harvey24_smf$logMstar50[1], harvey24_smf$alpha50[1]))
magplot(
  sm_bins, 
  z7_smf, 
  log = "y"
)

z7_sfrd_declining = trapz(
  sm_bins[sm_bins>8], z7_smf[sm_bins>8] * 10^(hf_sfh_declining$parm[1]*(sm_bins[sm_bins>8]) + hf_sfh_declining$parm[2])
)
z7_sfrd_rising = trapz(
  sm_bins[sm_bins>8], z7_smf[sm_bins>8] * 10^(hf_sfh_rising$parm[1]*(sm_bins[sm_bins>8]) + hf_sfh_rising$parm[2])
)
# z7_sfrd_declining = trapz(
#   sm_bins[sm_bins>=6 & sm_bins < 8], z7_smf[sm_bins>=6 & sm_bins < 8] * 10^(hf_sfh_declining$parm[1]*(sm_bins[sm_bins>=6 & sm_bins < 8]) + hf_sfh_declining$parm[2])
# )
# z7_sfrd_rising = trapz(
#   sm_bins[sm_bins>=8 & sm_bins < 12], z7_smf[sm_bins>=8 & sm_bins < 12] * 10^(hf_sfh_rising$parm[1]*(sm_bins[sm_bins>=8 & sm_bins < 12]) + hf_sfh_rising$parm[2])
# )


magplot(
  epochs$zbest,
  epochs$sfr_10Myr_pipes_zgauss, 
  log = "y", 
  ylim = c(1e-5, 1e4)
)
abline(0,1)


# x_sfr = seq(-3, 3, 0.01)
sm_bw = 0.25
sm_bins = seq(4.5, 12.5, sm_bw)
sm_mids = sm_bins[-length(sm_bins)] + diff(sm_bins) / 2

calc_Mstar_IQRSFR = function(sfr, sfr_err, mstar, mstar_err, sm_bins, Nsamples){
  ## Interquartile range?
  
  samples_ = foreach(i = 1:Nsamples, .combine = rbind) %do% {
    if(i %% 100 == 0){
      message(i)
    }
    x = rnorm(n = length(mstar), mean = mstar, sd = (mstar_err+1e-99))
    y = rnorm(n = length(sfr), mean = sfr, sd = (sfr_err+1e-99))
  
    calc_ = foreach(m = 1:(length(sm_bins)-1), .combine = c) %do% {
      midx = x >= sm_bins[m] & x < sm_bins[m+1] & y > 0
      
      quantile(log10(y[midx]), probs = 0.75, na.rm = TRUE) - quantile(log10(y[midx]), probs = 0.25, na.rm = TRUE)
      
      # sd(log10(y[midx]), na.rm = TRUE)
    }
    
  }
  
  df = data.frame(
    "Q50" = colQuantiles(samples_, probs = 0.5),
    "Q16" = colQuantiles(samples_, probs = 0.16),
    "Q84" = colQuantiles(samples_, probs = 0.84)
  )
  return( df )
}

iqr_mstar_rising = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zfix[zidx][rising_idx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zfix[zidx][rising_idx] + epochs$sfr_10Myr_pipes_l1_zfix[zidx][rising_idx]), 
  epochs$stellar_mass_pipes_zfix[zidx][rising_idx], 
  0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx][rising_idx] + epochs$stellar_mass_pipes_l1_zfix[zidx][rising_idx]),
  sm_bins, 
  100
)
iqr_mstar_declining = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zfix[zidx][declining_idx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zfix[zidx][declining_idx] + epochs$sfr_10Myr_pipes_l1_zfix[zidx][declining_idx]), 
  epochs$stellar_mass_pipes_zfix[zidx][declining_idx], 
  0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx][declining_idx] + epochs$stellar_mass_pipes_l1_zfix[zidx][declining_idx]),
  sm_bins, 
  100
)
iqr_mstar = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zfix[zidx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zfix[zidx] + epochs$sfr_10Myr_pipes_l1_zfix[zidx]), 
  epochs$stellar_mass_pipes_zfix[zidx], 
  0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx] + epochs$stellar_mass_pipes_l1_zfix[zidx]),
  sm_bins, 
  100
)


png("~/Documents/Epochs/plots/rough-plots/Mstar-sigmaSFR.png", width = 8, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0), mfrow = c(1,1))
magplot(
  sm_mids, iqr_mstar_rising$Q50, col = "blue", 
  xlim = c(7.5, 10.5), ylim = c(0, 1.0), 
  xlab = "logMstar", ylab = "IQR SFR10",
  main = paste(zmin, "< z < ", zmax)
)
magerr(
  sm_mids, iqr_mstar_rising$Q50, 
  ylo = iqr_mstar_rising$Q50 - iqr_mstar_rising$Q16,
  yhi = iqr_mstar_rising$Q84 - iqr_mstar_rising$Q50, 
  col = alpha("blue", 0.6)
)
points(
  sm_mids, iqr_mstar_declining$Q50, col = "red"
)
magerr(
  sm_mids, iqr_mstar_declining$Q50, 
  ylo = iqr_mstar_declining$Q50 - iqr_mstar_declining$Q16,
  yhi = iqr_mstar_declining$Q84 - iqr_mstar_declining$Q50, 
  col = alpha("red", 0.6)
)
lines(
  sm_mids, iqr_mstar$Q50
)
lines(
  sm_mids, iqr_mstar$Q16, lty = 2
)
lines(
  sm_mids, iqr_mstar$Q84, lty = 2
)
legend(
  x = "topleft",
  col = c("cornflowerblue", "red", "black"), 
  pch = 16,
  legend = c("SFR10/SFR100>2", "SFR10/SFR100<1", "Total")
)
dev.off()


