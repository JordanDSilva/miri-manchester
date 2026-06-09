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

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_survey"
fname_append = "_larson"

fittype = "photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth"
astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
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

eazy_stats = data.frame(fread(paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_stats.csv")))

# selection_idx = eazy_stats$zbestPzIntegral >= 0.6 & eazy_stats$zbestChi2Pval > 0.005
selection_idx = eazy_stats$zbestPzIntegral >= 0.6

template_cols = names(
  readRDS(
    list.files(paste0(main_dir, fname, "/prospect/", fittype), full.names = TRUE)[1]
  )$astro
)

# registerDoParallel(cores = 8)
# df_combine = foreach(ii = 1:dim(astro_cat)[1], .combine = rbind) %dopar% {
#   message("Running galID: ", ii)
#   read_stub = paste0(
#     main_dir, fname, "/prospect/", fittype, "/galID_",ii, ".rds"
#   )
#   
#   if(file.exists(read_stub)){
#     out = readRDS(read_stub)
#     
#     astro = out$astro
#     astro$ID = out$galInfo["ID"]
#     astro$DIC = out$galInfo["DIC"]
#     
#     rm(out)
#     gc()
#   }
#   else{
#     astro = setNames(as.data.frame(matrix(NA, nrow = 1, ncol = length(template_cols))), template_cols)
#     astro$ID = ii
#     astro$DIC = NA
#   }
#   return(astro)
# }
# stopImplicitCluster()

# fwrite(
#   data.frame(df_combine), 
#   paste0(main_dir, fname, "/prospect/", fittype, "_astro_cat.csv")
# )
# prospect_catalogue = data.frame(fread(
#   paste0(main_dir, fname, "/prospect/", fittype, "_astro_cat.csv")
# ))

lbt_bins = c(
  c(seq(0.0, 13.0, 0.5), 13.1)
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

SFR_bins = seq(-5, 5, 0.75)
SFR_mids = SFR_bins[1:(length(SFR_bins)-1)] + diff(SFR_bins)/2.0

Mstar_bins = seq(1, 14, 0.75)
Mstar_mids = Mstar_bins[1:(length(Mstar_bins)-1)] + diff(Mstar_bins)/2.0

make_completeness_func = function(inter_obj, xmids, ymids) {
  force(inter_obj)
  force(xmids)
  force(ymids)
  function(x, y) {
    x_clamped = pmax(min(xmids), pmin(max(xmids), x))
    y_clamped = pmax(min(ymids), pmin(max(ymids), y))
    interp.surface(inter_obj, cbind(x_clamped, y_clamped))
  }
}
calc_selection_3d = function(
    zbins, zmids, ztest, ztrue, 
    xtest, xtrue, xbins, xmids,
    ytest, ytrue, ybins, ymids,
    selection_idx
){
  Ndim = c(length(zmids), length(xmids), length(ymids))
  completeness_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  contamination_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  completeness_func_array = vector("list", length(zmids))
  
  for(i in seq_along(zmids)){
    true_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1])
    true_selected_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (selection_idx)
    
    xlevels = seq_along(xmids)
    ylevels = seq_along(ymids)
    
    xbin_true = cut(xtrue[true_idx], breaks = xbins, labels = FALSE)
    ybin_true = cut(ytrue[true_idx], breaks = ybins, labels = FALSE)
    
    xbin_sel  = cut(xtrue[true_selected_idx], breaks = xbins, labels = FALSE)
    ybin_sel  = cut(ytrue[true_selected_idx], breaks = ybins, labels = FALSE)
    
    # Force same factor levels so tables always have same dimensions
    n_true = table(
      factor(xbin_true, levels = xlevels),
      factor(ybin_true, levels = ylevels)
    )
    n_sel = table(
      factor(xbin_sel, levels = xlevels),
      factor(ybin_sel, levels = ylevels)
    )
    
    compl_mat =  n_sel / n_true
    compl_mat[n_true > 0 & n_true < 5] = 0
    compl_mat[is.infinite(compl_mat) | is.na(compl_mat)] = 0
    
    inter_obj = list(
      "x" = xmids, 
      "y" = ymids,
      "z" = compl_mat
    )
    inter_obj = list(x = xmids, y = ymids, z = matrix(as.numeric(compl_mat), nrow = length(xmids)))
    
    completeness_func = make_completeness_func(inter_obj, xmids, ymids)
    
    fine_x = seq(min(xbins), max(xbins), length.out = 200)
    fine_y = seq(min(ybins), max(ybins), length.out = 200)
    
    # Evaluate interpolation on every combination
    fine_grid = expand.grid(x = fine_x, y = fine_y)
    fine_compl = completeness_func(fine_grid$x, fine_grid$y)
    fine_mat = matrix(fine_compl, nrow = length(fine_x), ncol = length(fine_y))
    
    magplot(
      log10(astro_cat$StellarMass[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1]]),
      log10(astro_cat$SFR100[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1]]),
      xlim = range(xbins),
      ylim = range(ybins), 
      col = adjustcolor("black", 0.4),
      pch = 16
    )
    points(
      log10(astro_cat$StellarMass[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1] & selection_idx == 1]),
      log10(astro_cat$SFR100[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1] & selection_idx == 1]),
      col = adjustcolor("magenta", 0.4), 
      pch = 16
    )
    contour(
      x = fine_x,
      y = fine_y,
      z = fine_mat,
      add = TRUE,
      levels = c(0.1, 0.3, 0.5, 0.7, 0.9),
      col = "red",
      lwd = 1.5,
      labcex = 1.0
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3))
      )
    )
    completeness_array[i,,] = compl_mat
    completeness_func_array[[i]] = completeness_func
  }
  return( list("completeness" = completeness_array, "contamination" = contamination_array, "zmids" = zmids, "xmids" = xmids, "ymids" = ymids, "compl_funcs" = completeness_func_array) )
}

foo = calc_selection_3d(
  zbins = zbins,
  zmids = zmids,
  ztest = astro_cat$z,
  ztrue = astro_cat$z,
  ybins = SFR_bins,
  ymids = SFR_bins,
  ytrue = log10(astro_cat$SFR100),
  ytest = log10(astro_cat$SFR100),
  xbins = Mstar_bins,
  xmids = Mstar_mids,
  xtrue = log10(astro_cat$StellarMass),
  xtest = log10(astro_cat$StellarMass),
  selection_idx = selection_idx == 1
)

# message("Compute Veff from integral of dz/[dVc dOmega]")
# OmegaM = 0.30
# OmegaL = 0.70
# H0 = 70
# Ez = function(z){
#   sqrt(
#     OmegaM * (1+z)^3 + OmegaL
#   )
# }
# DH = 3000 * (H0/100)^-1
# dVc = function(z){
#   DH * (1+z)^2 * (cosdistAngDist(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = 0))^2 / Ez(z)
# } ## differntial comoving volume element per dz and dOmega[steradian]
# 
# Veff_func = function(Mstar, SFR, zmin, zmax, area, compl_func){
#   pracma::integral(
#     fun = dVc, 
#     xmin = zmin,
#     xmax = zmax
#   ) * area * compl_func(Mstar, SFR)
# }
# 
# Veff_func_z = function(x){
#   Mstar = x[,1]
#   SFR = x[,2]
#   ff = foo$compl_funcs[[1]](Mstar, SFR)
#   return(ff)
# }

kk = 5
zidx = astro_cat$z > zbins[kk] & astro_cat$z < zbins[kk+1]
ztemp = log10(astro_cat[zidx, c("StellarMass", "SFR100")])

eff_vols = foo$compl_funcs[[kk]](ztemp$StellarMass, ztemp$SFR100)
Ndim = c(length(Mstar_mids), length(SFR_mids))
temp = array(data = rep(0, prod(Ndim)), dim = Ndim)
temp_err = array(data = rep(0, prod(Ndim)), dim = Ndim)
for(i in seq_along(Mstar_mids)){
  for(j in seq_along(SFR_mids)){
    idx = ztemp$StellarMass >= Mstar_bins[i] & ztemp$StellarMass < Mstar_bins[i+1] & ztemp$SFR100 >= SFR_bins[j] & ztemp$SFR100 < SFR_bins[j+1]
    # temp[i,j] = sum(1/eff_vols[idx])
    temp[i,j] = sum(idx / foo$completeness[kk, i, j])
    temp_err[i,j] = sqrt(sum((idx / foo$completeness[kk, i, j])^2))
  }
}

hh = maghist(log10(astro_cat$StellarMass[zidx]), breaks = Mstar_bins, log = "y", plot = FALSE)
hh_raw_select = maghist(log10(astro_cat$StellarMass[zidx & selection_idx]), breaks = Mstar_bins, log = "y", plot = FALSE)
magplot(
  hh$mids, 
  hh$counts, 
  type = "s", 
  log = "y",
  ylim = c(1, 1500)
)
lines(Mstar_mids, rowSums(temp, na.rm = TRUE) * diff(Mstar_bins), type = 's', col = "blue")
lines(Mstar_mids, hh_raw_select$counts, type = 's', col = "red")

for(kk in seq_along(zmids)){
  image(
    x = foo$xmids,
    y = foo$ymids,
    z = foo$completeness[kk,,],
    xlab = "Mstar",
    ylab = "SFR",
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    ylim = c(-5, 5),
    xlim = c(4, 13)
  )
  image.plot(
    legend.only = TRUE,
    zlim = c(0, 1),
    # zlim = c(min(foo$completeness[kk,,], na.rm = TRUE), max(foo$completeness[kk,,], na.rm = TRUE)),
    col = hcl.colors(100, "YlOrRd", rev = TRUE),
    legend.lab = "Completeness"
  )
  points(
    log10(astro_cat$StellarMass[astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk+1]]),
    log10(astro_cat$SFR100[astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk+1]])
  )
  points(
    log10(astro_cat$StellarMass[astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk+1] & selection_idx == 1]),
    log10(astro_cat$SFR100[astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk+1] & selection_idx == 1]),
  col = "magenta"
  )
  legend(
    x = "topleft",
    legend = c(
      paste0("z = ", round(zmids[kk],3))
    )
  )
}

## old 
# calc_selection = function(
#     zbins, zmids, ztest, ztrue, 
#     xtest, xtrue, xbins, xmids,
#     selection_idx
# ){
#   completeness_mat = matrix(
#     data = 0, 
#     ncol = length(xmids),
#     nrow = length(zmids)
#   )
#   contamination_mat = matrix(
#     data = 0, 
#     ncol = length(xmids),
#     nrow = length(zmids)
#   )
#   
#   for(i in seq_along(zmids)){
#     for(j in seq_along(xmids)){
#       
#       # true, known injected in bin
#       true_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (xtrue >= xbins[j]) & (xtrue < xbins[j+1])
#       # of the true injected, how many do I select ?
#       true_selected_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) & (selection_idx)
#       # how many recovered do I select
#       check_idx = (ztest >= zbins[i]) & (ztest < zbins[i+1]) & (xtest >= xbins[j]) & (xtest < xbins[j+1]) & (selection_idx)
#       
#       completeness = sum(true_selected_idx, na.rm = TRUE) / sum(true_idx, na.rm = TRUE)
#       purity = length( intersect(which(true_selected_idx), which(check_idx)) ) / sum(check_idx, na.rm = TRUE)
#       
#       if(is.na(completeness) | is.infinite(completeness) | sum(true_idx) <= 10){
#         completeness = NA
#       }
#       if(is.na(purity) | is.infinite(purity)){
#         purity = NA
#       }
#       
#       completeness_mat[i,j] = completeness
#       contamination_mat[i,j] = 1 - purity
#     }
#   }
#   return( list("completeness" = completeness_mat, "contamination" = contamination_mat, "zmids" = zmids, "xmids" = xmids) )
# }
# 
# compl_quantiles = function(compl_mat){
#   m = length(compl_mat[[1]]$zmids)
#   n = length(compl_mat[[1]]$xmids)
#   
#   compl_arr_mc <- array(unlist(lapply(compl_mat, function(x)x$completeness)), dim = c(m, n, length(compl_mat)))
#   # rowMedians works on a matrix, so reshape to (m*n) x 100
#   compl_mc <- matrix(rowMedians(matrix(compl_arr_mc, nrow = m*n)), nrow = m, ncol = n)
#   compl_mc_Q16 <- matrix(rowQuantiles(matrix(compl_arr_mc, nrow = m*n), probs = 0.16), nrow = m, ncol = n)
#   compl_mc_Q84 <- matrix(rowMedians(matrix(compl_arr_mc, nrow = m*n), probs = 0.84), nrow = m, ncol = n)
#   
#   cont_arr_mc <- array(unlist(lapply(compl_mat, function(x)x$contamination)), dim = c(m, n, length(compl_mat)))
#   # rowMedians works on a matrix, so reshape to (m*n) x 100
#   cont_mc <- matrix(rowMedians(matrix(cont_arr_mc, nrow = m*n)), nrow = m, ncol = n)
#   cont_mc_Q16 <- matrix(rowQuantiles(matrix(cont_arr_mc, nrow = m*n), probs = 0.16), nrow = m, ncol = n)
#   cont_mc_Q84 <- matrix(rowQuantiles(matrix(cont_arr_mc, nrow = m*n), probs = 0.84), nrow = m, ncol = n)
#   
#   out = list(
#     "z" = compl_mat[[1]]$zmids,
#     "x" = compl_mat[[1]]$xmids,
#     "compl" = compl_mc,
#     "compl_Q16" = compl_mc_Q16,
#     "compl_Q84" = compl_mc_Q84,
#     
#     "cont" = cont_mc,
#     "cont_Q16" = cont_mc_Q16,
#     "cont_Q84" = cont_mc_Q84
#   )
#   return(out)
# }
# MUV_completeness = compl_quantiles(lapply(
#   1:100,
#   function(x){
#     foo = calc_selection(
#       zbins = zbins,
#       zmids = zmids,
#       # ztest = eazy$z_ml,
#       ztest = rnorm(
#         n=length(prospect_catalogue$zQ50),
#         mean = prospect_catalogue$zQ50,
#         sd = 0.5*(prospect_catalogue$z_lo + prospect_catalogue$z_hi)
#       ),
#       ztrue = astro_cat$z,
#       xbins = MUV_bins,
#       xmids = MUV_mids,
#       xtrue = astro_cat$UV1500,
#       xtest = rnorm(
#         n = length(prospect_catalogue$StellarMassQ50),
#         mean = prospect_catalogue$UV1500Q50,
#         sd = 0.5*(prospect_catalogue$UV1500_lo + prospect_catalogue$UV1500_hi)
#       ),
#       selection_idx = selection_idx
#     )
#   }
# ))
# 
# Mstar_completeness = compl_quantiles(lapply(
#   1:100,
#   function(x){
#     foo = calc_selection(
#       zbins = zbins,
#       zmids = zmids,
#       # ztest = eazy$z_ml,
#       ztest = rnorm(
#         n=length(prospect_catalogue$zQ50),
#         mean = prospect_catalogue$zQ50,
#         sd = 0.5*(prospect_catalogue$z_lo + prospect_catalogue$z_hi)
#       ),
#       ztrue = astro_cat$z,
#       xbins = Mstar_bins,
#       xmids = Mstar_mids,
#       xtrue = log10(astro_cat$StellarMass),
#       xtest = (rnorm(
#         n = length(prospect_catalogue$StellarMassQ50),
#         mean = log10(prospect_catalogue$StellarMassQ50),
#         sd = 0.5 * (log10(prospect_catalogue$StellarMassQ84) - log10(prospect_catalogue$StellarMassQ16))
#       )),
#       selection_idx = selection_idx
#     )
#   }
# ))
# par(mfrow = c(2,2))
# image(
#   x = MUV_completeness$z,
#   y = MUV_completeness$x,
#   z = MUV_completeness$compl,
#   xlab = "z",
#   ylab = expression(M[UV]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   # xlim = c(0,16)
# )
# image.plot(
#   legend.only = TRUE,
#   zlim = c(min(MUV_completeness$compl, na.rm = TRUE), max(MUV_completeness$compl, na.rm = TRUE)),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   legend.lab = "Completeness"
# )
# 
# image(
#   x = MUV_completeness$z,
#   y = MUV_completeness$x,
#   z = MUV_completeness$cont,
#   xlab = "z",
#   ylab = expression(M[UV]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   # xlim = c(0,16)
# )
# image.plot(
#   legend.only = TRUE,
#   zlim = c(min(MUV_completeness$cont, na.rm = TRUE), max(MUV_completeness$cont, na.rm = TRUE)),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   legend.lab = "Completeness"
# )
# 
# image(
#   x = Mstar_completeness$z,
#   y = Mstar_completeness$x,
#   z = Mstar_completeness$compl,
#   xlab = "z",
#   ylab = expression(M[star]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   # xlim = c(0,16)
# )
# image.plot(
#   legend.only = TRUE,
#   zlim = c(min(Mstar_completeness$compl, na.rm = TRUE), max(Mstar_completeness$compl, na.rm = TRUE)),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   legend.lab = "Completeness"
# )
# 
# image(
#   x = Mstar_completeness$z,
#   y = Mstar_completeness$x,
#   z = Mstar_completeness$cont,
#   xlab = "z",
#   ylab = expression(M[star]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   # xlim = c(0,16)
# )
# image.plot(
#   legend.only = TRUE,
#   zlim = c(min(Mstar_completeness$cont, na.rm = TRUE), max(Mstar_completeness$cont, na.rm = TRUE)),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   legend.lab = "Completeness"
# )
# 
# 
# image(
#   x = Mstar_completeness$z,
#   y = Mstar_completeness$x,
#   z = Mstar_completeness$compl,
#   xlab = "z",
#   ylab = expression(M[star]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   xlim = c(8,12),
#   ylim = c(8,10.5)
# )
# image(
#   x = Mstar_completeness$z,
#   y = Mstar_completeness$x,
#   z = Mstar_completeness$cont,
#   xlab = "z",
#   ylab = expression(M[star]),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   xlim = c(8,12),
#   ylim = c(8,10.5)
# )
# image.plot(
#   legend.only = TRUE,
#   zlim = c(min(Mstar_completeness$compl, na.rm = TRUE), max(Mstar_completeness$compl, na.rm = TRUE)),
#   col = hcl.colors(100, "YlOrRd", rev = TRUE),
#   legend.lab = "Completeness"
# )
# points(
#   astro_cat$z,
#   log10(astro_cat$StellarMass),
#   col = "black", 
#   pch = "."
# )
# points(
#   astro_cat$z[selection_idx],
#   log10(astro_cat$StellarMass[selection_idx]),
#   col = "magenta", 
#   pch = "."
# )