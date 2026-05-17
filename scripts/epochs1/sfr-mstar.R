library(Rfits)
library(magicaxis)
library(ProSpect)
library(celestial)
library(scales)
library(dftools)
library(foreach)
library(hyper.fit)
library(matrixStats)

## scatter in SFSs

epochs = Rfits_read_table("~/Documents/Epochs/data/EPOCHS_MASTER_public_final.fits")

## Mstar-SFR
zmin = 6.5
zmax = 7.5

zidx = epochs$redshift_pipes_zgauss >= zmin & epochs$redshift_pipes_zgauss < zmax

Mstar_lim = 1

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
    x = rnorm(n = length(mstar), mean = mstar, sd = (mstar_err))
    y = rnorm(n = length(sfr), mean = log10(sfr), sd = (sfr_err/(log(10)*sfr)))
    
    calc_ = foreach(m = 1:(length(sm_bins)-1), .combine = c) %do% {
      midx = x >= sm_bins[m] & x < sm_bins[m+1]
      
      quantile(y[midx], probs = 0.75, na.rm = TRUE) - quantile(y[midx], probs = 0.25, na.rm = TRUE)
      
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

iqr_rising = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zgauss[zidx][rising_idx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zgauss[zidx][rising_idx] + epochs$sfr_10Myr_pipes_l1_zgauss[zidx][rising_idx]), 
  epochs$stellar_mass_pipes_zgauss[zidx][rising_idx], 
  0.5*(epochs$stellar_mass_pipes_u1_zgauss[zidx][rising_idx] + epochs$stellar_mass_pipes_l1_zgauss[zidx][rising_idx]),
  sm_bins, 
  100
)
iqr_intermediate = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zgauss[zidx][intermediate_idx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zgauss[zidx][intermediate_idx] + epochs$sfr_10Myr_pipes_l1_zgauss[zidx][intermediate_idx]), 
  epochs$stellar_mass_pipes_zgauss[zidx][intermediate_idx], 
  0.5*(epochs$stellar_mass_pipes_u1_zgauss[zidx][intermediate_idx] + epochs$stellar_mass_pipes_l1_zgauss[zidx][intermediate_idx]),
  sm_bins, 
  100
)
iqr_declining = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zgauss[zidx][declining_idx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zgauss[zidx][declining_idx] + epochs$sfr_10Myr_pipes_l1_zgauss[zidx][declining_idx]), 
  epochs$stellar_mass_pipes_zgauss[zidx][declining_idx], 
  0.5*(epochs$stellar_mass_pipes_u1_zgauss[zidx][declining_idx] + epochs$stellar_mass_pipes_l1_zgauss[zidx][declining_idx]),
  sm_bins, 
  100
)



# png("~/Documents/Epochs/plots/rough-plots/Mstar-sigmaSFR.png", width = 8, height = 6, units = "in", res = 240)
par(mar = c(3.5, 3.5, 2.5, 1.0), mfrow = c(1,1))
magplot(
  sm_mids, iqr_rising$Q50, col = "blue", 
  xlim = c(7.5, 10.5), ylim = c(0, 1.0), 
  xlab = "logMstar", ylab = "IQR SFR10",
  main = paste(zmin, "< z < ", zmax)
)
magerr(
  sm_mids, iqr_rising$Q50, 
  ylo = iqr_rising$Q50 - iqr_rising$Q16,
  yhi = iqr_rising$Q84 - iqr_rising$Q50, 
  col = alpha("blue", 0.6)
)
points(
  sm_mids, iqr_intermediate$Q50, col = "forestgreen"
)
magerr(
  sm_mids, iqr_intermediate$Q50, 
  ylo = iqr_intermediate$Q50 - iqr_intermediate$Q16,
  yhi = iqr_intermediate$Q84 - iqr_intermediate$Q50, 
  col = alpha("forestgreen", 0.6)
)
points(
  sm_mids, iqr_declining$Q50, col = "red"
)
magerr(
  sm_mids, iqr_declining$Q50, 
  ylo = iqr_declining$Q50 - iqr_declining$Q16,
  yhi = iqr_declining$Q84 - iqr_declining$Q50, 
  col = alpha("red", 0.6)
)
lines(
  sm_mids, iqr_total$Q50
)
lines(
  sm_mids, iqr_total$Q16, lty = 2
)
lines(
  sm_mids, iqr_total$Q84, lty = 2
)
legend(
  x = "topleft",
  col = c("cornflowerblue", "forestgreen", "red", "black"), 
  pch = 16,
  legend = c("SFR10/SFR100>2", "1<SFR10/SFR100<2", "SFR10/SFR100<1", "Total")
)
# dev.off()


magplot(
  epochs$stellar_mass_pipes_zfix, 
  epochs$stellar_mass_pipes_zgauss
)
abline(0,1)
magplot(
  epochs$zbest,
  epochs$redshift_pipes_zgauss
)

par(mfrow = c(2,1))
magplot(
  epochs$stellar_mass_pipes_zfix, 
  epochs$sfr_10Myr_pipes_zfix, 
  z = epochs$sfr_10Myr_pipes_zfix/epochs$sfr_100Myr_pipes_zfix,
  log = "y", 
  # zcol = alpha(c("cornflowerblue", "forestgreen", "red", "darkgrey"), 0.5),
  dobar = FALSE, 
  pch = 16,
  xlim = c(6.5, 11.0), 
  ylim = c(0.1, 300),
  xlab = "log10(Mstar[Msol])",
  ylab = "SFR10 [Msol yr^-1]", 
  main = paste(zmin, "< z < ", zmax)
)
magplot(
  epochs$stellar_mass_pipes_zgauss, 
  epochs$sfr_10Myr_pipes_zgauss, 
  z = epochs$sfr_10Myr_pipes_zgauss/epochs$sfr_100Myr_pipes_zgauss,
  log = "y", 
  # zcol = alpha(c("cornflowerblue", "forestgreen", "red", "darkgrey"), 0.5),
  dobar = FALSE, 
  pch = 16,
  xlim = c(6.5, 11.0), 
  ylim = c(0.1, 300),
  xlab = "log10(Mstar[Msol])",
  ylab = "SFR10 [Msol yr^-1]", 
  main = paste(zmin, "< z < ", zmax)
)


zmin = 7.5
zmax = 8.5
zidx = epochs$redshift_pipes_zgauss >= zmin & epochs$redshift_pipes_zgauss < zmax
iqr_total = calc_Mstar_IQRSFR(
  sfr = epochs$sfr_10Myr_pipes_zgauss[zidx], 
  sfr_err = 0.5*(epochs$sfr_10Myr_pipes_u1_zgauss[zidx] + epochs$sfr_10Myr_pipes_l1_zgauss[zidx]), 
  mstar = epochs$stellar_mass_pipes_zgauss[zidx], 
  mstar_err = 0.5*(epochs$stellar_mass_pipes_u1_zgauss[zidx] + epochs$stellar_mass_pipes_l1_zgauss[zidx]),
  sm_bins = sm_bins, 
  Nsamples = 100
)
zidx = epochs$zbest >= zmin & epochs$zbest < zmax
iqr_total_zfix = calc_Mstar_IQRSFR(
  epochs$sfr_10Myr_pipes_zfix[zidx], 
  0.5*(epochs$sfr_10Myr_pipes_u1_zfix[zidx] + epochs$sfr_10Myr_pipes_l1_zfix[zidx]), 
  epochs$stellar_mass_pipes_zfix[zidx], 
  0.5*(epochs$stellar_mass_pipes_u1_zfix[zidx] + epochs$stellar_mass_pipes_l1_zfix[zidx]),
  sm_bins, 
  100
)
magplot(
  sm_mids, iqr_total$Q50, col = "black", 
  xlim = c(7.5, 10.5), ylim = c(0, 1.0), 
  xlab = "logMstar", ylab = "IQR SFR10",
  main = paste(zmin, "< z < ", zmax)
)
magerr(
  sm_mids, iqr_total$Q50, 
  ylo = iqr_total$Q50 - iqr_total$Q16,
  yhi = iqr_total$Q84 - iqr_total$Q50, 
  col = alpha("black", 0.6)
)
points(
  sm_mids, iqr_total_zfix$Q50, col = "black", pch = 16
)
magerr(
  sm_mids, iqr_total_zfix$Q50, 
  ylo = iqr_total_zfix$Q50 - iqr_total_zfix$Q16,
  yhi = iqr_total_zfix$Q84 - iqr_total_zfix$Q50, 
  col = alpha("black", 0.6), 
  pch = 16
)
