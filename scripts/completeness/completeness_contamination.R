library(numDeriv)
library(Rfits)
library(data.table)
library(magicaxis)
library(foreach)
library(pracma)
library(Highlander)
library(rhdf5)
library(celestial)
library(fields)
library(MASS)

## Load in the 100 000 survey galaxies with magnitudes following the same distribution I observe (synthetic observation)
main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_survey"
fname_append = "_agn_blue_sfhz"

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv")
))
astro_cat_nn = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_nn_predictions.csv")
))

eazy_zPDFs = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_mock_photometry_eazy_zPDF",fname_append, ".csv")
))
eazy = Rfits_read_table(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot", fname_append, ".out.zout.fits")
)
eazy_ntemplates = h5read(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot", fname_append, ".out.h5"),
  'fit/coeffs_best'
)
eazy_ntemplates = as.numeric(colSums(
  eazy_ntemplates > 0
))

compute_sigmaNMAD <- function(zspec, zphot, zphot_lim = 0) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz   <- zphot[idx] - zspec[idx]
  
  sigmaNMAD <- 1.48 * median( abs( (deltaz - median(deltaz)) / (1 + zspec[idx]) ) )
  return(sigmaNMAD)
}
compute_outliers <- function(zspec, zphot, thresh = 0.15, zphot_lim = -2) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz     <- abs(zphot[idx] - zspec[idx])
  Noutlier   <- sum(deltaz/(1 + zspec[idx]) > thresh)
  outlier_frac <- Noutlier / length(zspec[idx])
  
  cat("Noutlier:", Noutlier, "\n")
  cat("N total: ", length(zspec[idx]), "\n")
  
  return(outlier_frac)
}
zspec_zphot_plot = function(zphot, zspec, zphot_lim = 0){
  magplot(
    zphot,
    zspec, 
    xlim = c(0, 40),
    ylim = c(0, 40),
    xlab = "zphot",
    ylab = "zspec",
    pch = "."
  )
  legend(
    x = "topright",
    legend = c(
      paste0("sigmaNMAD = ", round(compute_sigmaNMAD(zspec, zphot, zphot_lim = zphot_lim), 3)),
      paste0("Outliers = ", round(compute_outliers(zspec, zphot, zphot_lim = zphot_lim), 3))
    )
  )
}

zspec_zphot_plot(
  zphot = astro_cat_nn$z,
  zspec = astro_cat$z,
  zphot_lim = 10
)
zspec_zphot_plot(
  zphot = eazy$z_ml,
  zspec = astro_cat$z
)

magplot(
  log10(astro_cat$StellarMass),
  log10(astro_cat_nn$z), 
  pch = ".",
  col = adjustcolor("black", 0.3)
)

maghist(
  log10(astro_cat$StellarMass)-
  log10(astro_cat_nn$StellarMass)
)
maghist(
  ((astro_cat$z)-(astro_cat_nn$z))/(1+astro_cat$z)
)

compute_peak_integrals = function(eazy_zPDFs, eazy, kk){
  zgrid = eazy_zPDFs$zgrid
  zPDF = eazy_zPDFs[[paste0('gal',kk)]]
  
  zbest = eazy$z_ml[kk]
  
  idx = zgrid >= 0.9*zbest & zgrid <= 1.1*zbest
  trapz(
    zgrid[idx],
    zPDF[idx]
  )
}

zbestPzIntegral = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(eazy_zPDFs, eazy, x))
selection_idx = rep(TRUE, dim(astro_cat)[1])

# zbins = seq(0.0, 40.0, 1.0)
# zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

lbt_bins = c(seq(8.0, 13.0, 0.5), 13.2)
zbins = approx(
  x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
  y = seq(0,40,1),
  xout = lbt_bins
)$y
# zbins = seq(0.0, 12.0, 0.5)
zbins = c(0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10.0, 12.0)
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

Mstar_bins = seq(6.5, 12.5, 0.5)
Mstar_mids = Mstar_bins[1:(length(Mstar_bins)-1)] + diff(Mstar_bins)/2.0

MUV_bins = seq(-32, -10, 1.0)
MUV_mids = MUV_bins[1:(length(MUV_bins)-1)] + diff(MUV_bins)/2.0

calc_selection_z_M = function(
    zbins, zmids, ztest, ztrue, 
    xtest, xtrue, xbins, xmids,
    selection_idx
){
  Ndim = c(length(zmids), length(xmids))
  completeness_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  contamination_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  
  for(i in seq_along(zmids)){
    for(j in seq_along(xmids)){
      true_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (xtrue >= xbins[j]) & (xtrue < xbins[j+1])
      selected_idx = (ztest >= zbins[i]) & (ztest < zbins[i+1]) & (xtest >= xbins[j]) & (xtest < xbins[j+1]) & (selection_idx)
      
      # which_in = which(true_idx)
      which_in = astro_cat$ID[true_idx]
      which_select = astro_cat$ID[selected_idx]
      
      # which_select = which(selected_idx)
      which_exact = intersect(
        which_in, which_select
      )
      
      if(sum(true_idx) == 0){
        completeness_array[i,j] = NA
        contamination_array[i,j] = 1
      }else{
        # completeness_array[i,j] = sum(selected_idx) / sum(true_idx)
        completeness_array[i,j] = length(which_exact) / sum(true_idx)
        contamination_array[i,j] = 1 - (length(which_exact) / sum(selected_idx))
      }
      
      if(sum(selected_idx) == 0){
        contamination_array[i,j] = NA
      } else {
        contamination_array[i,j] = 1 - (length(which_exact) / sum(selected_idx))
      }
      
    }
  }
  
  out = list(
    "completeness" = completeness_array, 
    "contamination" = contamination_array, 
    "zbins" = zbins,
    "xbins" = xbins,
    "zmids" = zmids,
    "xmids" = xmids
  )
  return( out )
}

foo_Mstar = calc_selection_z_M(
  zbins = zbins,
  zmids = zmids,
  ztest = astro_cat_nn$z,
  ztrue = astro_cat_nn$z,
  xbins = Mstar_bins,
  xmids = Mstar_mids,
  xtrue = log10(astro_cat$StellarMass),
  xtest = log10(astro_cat_nn$StellarMass),
  # selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  selection_idx = selection_idx == 1L
)
foo_MUV = calc_selection_z_M(
  zbins = zbins,
  zmids = zmids,
  ztest = astro_cat_injected$z,
  ztrue = astro_cat_injected$z,
  xbins = MUV_bins,
  xmids = MUV_mids,
  xtrue = (astro_cat_injected$UV1500),
  xtest = (astro_cat_nn$UV1500),
  # selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  selection_idx = selection_idx == 1L
)

plot_matrices = function(foo){
  par(mfrow = c(1,2), mar = rep(3.0,4), oma = rep(2.0, 4))
  compl_mat = foo$completeness
  cont_mat = foo$contamination
  image(
    x = foo$zmids,
    y = foo$xmids,
    z = compl_mat,
    xlab = "z",
    ylab = "MUV",
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    xlim = range(foo$zmids),
    ylim = range(foo$xmids),
  )
  image.plot(
    legend.only = TRUE,
    zlim = range(compl_mat, na.rm = TRUE),
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    legend.lab = "Completeness"
  )
  
  image(
    x = foo$zmids,
    y = foo$xmids,
    z = cont_mat,
    xlab = "z",
    ylab = "MUV",
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    xlim = range(foo$zmids),
    ylim = range(foo$xmids),
  )
  image.plot(
    legend.only = TRUE,
    zlim = range(cont_mat, na.rm = TRUE),
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    legend.lab = "Completeness"
  )
}

compl_mat = foo_Mstar$completeness
cont_mat = foo_Mstar$contamination
saveRDS(1-cont_mat, "~/Desktop/cont_mat_test.rds")
saveRDS(compl_mat, "~/Desktop/comp_mat_test.rds")

plot_matrices(foo_MUV)
plot_matrices(foo_Mstar)

kk = 1
magplot(
  Mstar_mids,
  compl_mat[kk,],
  ylim = c(0,1),
  main = zmids[kk]
)
