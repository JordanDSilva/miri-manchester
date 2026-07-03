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
library(mgcv)

## Load in the 100 000 survey galaxies with magnitudes following the same distribution I observe (synthetic observation)
main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_uniform_detect"
fname_append = "_agn_blue_sfhz"

fittype = "photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth"

astro_cat_injected = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))
phot_cat_injected = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv")
))

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat_recovered.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry_recovered.csv")
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

compute_peak_integrals = function(eazy_zPDFs, eazy, kk){
  print(kk)
  zgrid = eazy_zPDFs$zgrid
  zPDF = eazy_zPDFs[[paste0('gal',kk)]]
  
  zbest = eazy$z_ml[kk]
  
  idx = zgrid >= 0.9*zbest & zgrid <= 1.1*zbest
  trapz(
    zgrid[idx],
    zPDF[idx]
  )
}

zbestPzIntegral = sapply(1:dim(astro_cat_injected)[1], function(x)compute_peak_integrals(eazy_zPDFs, eazy, x))

selection_idx = astro_cat_injected$ID %in% astro_cat$ID

# magplot(
#   (astro_cat_injected$z[selection_idx]),
#   (eazy$z_ml[selection_idx]),
#   z = eazy$z_ml[selection_idx],
#   xlim = c(0, 15),
#   ylim = c(0, 15),
#   pch = "."
# )
# magplot(
#   log10(astro_cat_injected$StellarMass[selection_idx]),
#   log10(astro_cat_nn$StellarMass[selection_idx]),
#   xlim = c(6, 12),
#   ylim = c(6, 12),
#   pch = "."
# )
# maghist(
#   log10(astro_cat_injected$StellarMass)[log10(astro_cat_nn$StellarMass)>9]-
#     log10(astro_cat_nn$StellarMass)[log10(astro_cat_nn$StellarMass)>9],
# )
# magplot(
#   (astro_cat_injected$UV1500[selection_idx]),
#   (astro_cat_nn$UV1500[selection_idx]),
#   xlim = c(-30, -10),
#   ylim = c(-30, -10),
#   pch = "."
# )
# maghist(
#   (astro_cat_injected$UV1500[selection_idx])-
#   (astro_cat_nn$UV1500[selection_idx])
# )

zbins = seq(0.0, 40.0, 1.0)
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

Mstar_bins = seq(3, 13, 0.75)
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
      which_in = astro_cat_injected$ID[true_idx]
      which_select = astro_cat_injected$ID[selected_idx]
      
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
  ztest = eazy$z_ml,
  ztrue = astro_cat_injected$z,
  xbins = Mstar_bins,
  xmids = Mstar_mids,
  xtrue = log10(astro_cat_injected$StellarMass),
  xtest = log10(astro_cat_nn$StellarMass),
  # selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  selection_idx = selection_idx == 1L
)
foo_MUV = calc_selection_z_M(
  zbins = zbins,
  zmids = zmids,
  ztest = eazy$z_ml,
  ztrue = astro_cat_injected$z,
  xbins = MUV_bins,
  xmids = MUV_mids,
  xtrue = (astro_cat_injected$UV1500),
  xtest = (astro_cat_nn$UV1500),
  # selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  selection_idx = selection_idx == 1L
)

plot_matrices = function(foo){
  par(mfrow = c(1,2), mar = rep(1.0,4), oma = rep(1.0, 4))
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

plot_matrices(foo_MUV)
plot_matrices(foo_Mstar)

kk = 7
magplot(
  Mstar_mids,
  oce::matrixSmooth(compl_mat)[kk,], 
  xlim = c(6, 12),
  ylim = c(0,1)
)
points(
  Mstar_mids,
  compl_mat[kk,],
  pch = 16
)

# compl_mat[cont_mat >= 0.5 | is.na(compl_mat)] = 0

interp_func <- function(z_in, M_in, zmids = zmids, Mstar_mids = Mstar_mids, mat = compl_mat) {
  compl = fields::interp.surface(
    obj = list(x = zmids, y = Mstar_mids, z = mat),
    loc = cbind(z_in, M_in)
  )
  return(compl)
}

saveRDS(
  list("fn" = interp_func, zmids = zmids, Mstar_mids = Mstar_mids, compl_mat = (compl_mat), cont_mat = (cont_mat)),
  "~/Desktop/zooweemama.rds"
)
# Single value
interp_func(z_in = 1.5, M_in = 10.5, zmids = zmids, Mstar_mids = Mstar_mids, mat = compl_mat)

# grid <- expand.grid(
#   z     = seq(0, 12, length.out = 200),
#   Mstar = seq(6, 12, length.out = 200)
# )
# 
# # Use your interpolation function
# grid$completeness <- interp_func(grid$z, grid$Mstar)
# 
# # Reshape for contour
# completeness_grid <- matrix(grid$completeness, 
#                             nrow = 200, ncol = 200)
# contour(seq(0, 12, length.out = 200),
#         seq(6, 12, length.out = 200),
#         completeness_grid,
#         levels  = c(0.3, 0.5, 0.8, 1.0),
#         add     = TRUE,
#         col     = "dodgerblue",
#         lwd     = 2,
#         labcex  = 0.8)


# compl_results <- array(NA, dim = c(length(zmids), length(Mstar_mids), 100))
# cont_results <- array(NA, dim = c(length(zmids), length(Mstar_mids), 100))
# for(i in 1:100){
#   Mstar_test = rnorm(
#     length(log10(astro_cat$StellarMass)), 
#     log10(astro_cat$StellarMass), 
#     0.1
#   )
#   z_test = rnorm(
#     length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]),
#     eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID], 
#     0.1
#   )
# 
#   foo_Mstar = calc_selection_z_M(
#     zbins = zbins,
#     zmids = zmids,
#     ztest = z_test,
#     ztrue = astro_cat_injected$z,
#     xbins = Mstar_bins,
#     xmids = Mstar_mids,
#     xtrue = log10(astro_cat_injected$StellarMass),
#     xtest = Mstar_test,
#     # selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
#     selection_idx = zbestPzIntegral >= 0.6
#   )
# 
#   compl_mat = foo_Mstar$completeness
#   compl_results[,,i] = compl_mat
#   
#   cont_mat = foo_Mstar$contamination
#   cont_results[,,i] = cont_mat
# }
# 
# compl_mat <- apply(compl_results, c(1,2), function(x){median(x, na.rm = TRUE)})
# cont_mat <- apply(cont_results, c(1,2), function(x){median(x, na.rm = TRUE)})
