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

selection_idx = zbestPzIntegral >= 0.6 & astro_cat_injected$ID %in% astro_cat$ID

# lbt_bins = c(seq(8.0, 13.0, 0.5), 13.1)
# zbins = approx(
#   x = cosdistTravelTime(z = seq(0,40,1), ref = "737"),
#   y = seq(0,40,1),
#   xout = lbt_bins
# )$y
zbins = seq(0.0, 12.0, 1.0)
zmids = zbins[1:(length(zbins)-1)] + diff(zbins)/2

# MUV_bins = seq(-30, -10, 1)
# MUV_mids = MUV_bins[1:(length(MUV_bins)-1)] + diff(MUV_bins)/2.0
# 
# MAB_bins = seq(10, 40, 0.5)
# MAB_mids = MAB_bins[1:(length(MAB_bins)-1)] + diff(MAB_bins)/2.0
# 
# SFR_bins = seq(-5, 5, 0.5)
# SFR_mids = SFR_bins[1:(length(SFR_bins)-1)] + diff(SFR_bins)/2.0

Mstar_bins = seq(6, 12, 0.5)
Mstar_mids = Mstar_bins[1:(length(Mstar_bins)-1)] + diff(Mstar_bins)/2.0

# det_F560W_flux_true = as.numeric(phot_cat$F560W_true_flux)
# det_F770W_flux_true = as.numeric(phot_cat$F770W_true_flux)
# 
# det_F560W_flux_test = as.numeric(phot_cat$F560W_flux)
# det_F770W_flux_test = as.numeric(phot_cat$F770W_flux)
# 
# det_F560W_flux_err = as.numeric(phot_cat$F560W_flux_err)
# det_F770W_flux_err = as.numeric(phot_cat$F770W_flux_err)
# 
# det_mag_true = -2.5*log10(colMeans(
#   rbind(det_F560W_flux_true, det_F770W_flux_true)
# ))+8.9
# det_mag_test = -2.5*log10(colMeans(
#   rbind(pmax(0,det_F560W_flux_test), pmax(0,det_F770W_flux_test)),
# ))+8.9

# Build a binary outcome dataframe
# Every injected source: was it recovered?
df_all <- data.frame(
  z    = astro_cat_injected$z,
  Mstar = log10(pmax(1, astro_cat_injected$StellarMass)),
  recovered = astro_cat_injected$ID %in% astro_cat$ID  # TRUE/FALSE
)

## Fit smooth 2D surface
fit_gam <- gam(recovered ~ s(z, Mstar, k=25),
               data = df_all,
               family = binomial(link="logit"))

# Predict completeness on a fine grid
grid <- expand.grid(
  z     = zbins,
  Mstar = Mstar_bins
)
grid$completeness <- predict(fit_gam, newdata = grid, type = "response")
completeness_matrix <- matrix(grid$completeness, 
                              nrow = length(zbins), 
                              ncol = length(Mstar_bins))

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
      which_select = astro_cat$ID[selected_idx]
      
      # which_select = which(selected_idx)
      which_exact = intersect(
        which_in, which_select
      )
      
      if(sum(true_idx) == 0){
        completeness_array[i,j] = NA
        contamination_array[i,j] = NA
      }else{
        completeness_array[i,j] = sum(selected_idx) / sum(true_idx)
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

results <- array(NA, dim = c(length(zmids), length(Mstar_mids), 100))
for(i in 1:100){
  Mstar_test = rnorm(
    length(log10(astro_cat$StellarMass)), 
    log10(astro_cat$StellarMass), 
    0.1
  )
  z_test = rnorm(
    length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]),
    eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID], 
    0.2
  )
  # foo_Mstar = calc_selection_z_M(
  #   zbins = zbins,
  #   zmids = zmids,
  #   ztest = astro_cat$z,
  #   ztrue = astro_cat_injected$z,
  #   xbins = Mstar_bins,
  #   xmids = Mstar_mids,
  #   xtrue = log10(astro_cat_injected$StellarMass),
  #   xtest = log10(astro_cat$StellarMass),
  #   # selection_idx = zbestPzIntegral[astro_cat_injected$ID %in% astro_cat$ID] >= 0.6
  #   selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  # )
  
  foo_Mstar = calc_selection_z_M(
    zbins = zbins,
    zmids = zmids,
    ztest = z_test,
    ztrue = astro_cat_injected$z,
    xbins = Mstar_bins,
    xmids = Mstar_mids,
    xtrue = log10(astro_cat_injected$StellarMass),
    xtest = Mstar_test,
    selection_idx = rep(1L, length(eazy$z_ml[astro_cat_injected$ID %in% astro_cat$ID]))
  )
  
  compl_mat = foo_Mstar$completeness
  results[,,i] = compl_mat
}

compl_mat <- apply(results, c(1,2), median)
# compl_mat[compl_mat > 1] = NA

library(mgcv)
library(reshape2)

# Convert your completeness matrix to long format
df_binned <- melt(compl_mat, varnames = c("z_idx", "Mstar_idx"), value.name = "recovered_frac")

# Replace indices with actual bin midpoint values
df_binned$z     <- zmids[df_binned$z_idx]
df_binned$Mstar <- Mstar_mids[df_binned$Mstar_idx]

# Remove NAs (empty bins)
df_binned <- df_binned[is.finite(df_binned$recovered_frac) | !is.na(df_binned$recovered_frac), ]

# Beta regression requires values strictly between 0 and 1 (not 0 or 1 exactly)
# Squish any 0s and 1s slightly
df_binned$recovered_frac <- pmin(pmax(df_binned$recovered_frac, 0.001), 0.999)

# Fit
fit_gam <- gam(recovered_frac ~ s(z, Mstar, k=25),
               family = betar(link="logit"),
               data = df_binned)

summary(fit_gam)

grid <- expand.grid(
  z     = zbins,
  Mstar = Mstar_bins
)

grid$completeness <- predict(fit_gam, newdata = grid, type = "response")
completeness_matrix <- matrix(grid$completeness, 
                              nrow = length(zbins), 
                              ncol = length(Mstar_bins))

kk = 1
magplot(
  Mstar_mids, compl_mat[kk,], 
  xlim = c(6,12),
  ylim = c(0,2)
)
lines(
  seq(6,12,0.01), predict(fit_gam, newdata = data.frame(z = zmids[kk], Mstar = seq(6,12,0.01)), type = "response")
)

image(
  x = foo_Mstar$zmids,
  y = foo_Mstar$xmids,
  z = compl_mat,
  xlab = "z",
  ylab = "MUV",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  xlim = range(foo_Mstar$zmids),
  ylim = range(foo_Mstar$xmids),
)
image.plot(
  legend.only = TRUE,
  zlim = range(compl_mat, na.rm = TRUE),
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  legend.lab = "Completeness"
)
image(
  x = zbins,
  y = Mstar_bins,
  z = completeness_matrix,
  xlab = "z",
  ylab = "MUV",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  xlim = range(foo_Mstar$zmids),
  ylim = range(foo_Mstar$xmids),
)
image.plot(
  legend.only = TRUE,
  zlim = range(completeness_matrix, na.rm = TRUE),
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  legend.lab = "Completeness"
)

image(
  x = foo_Mstar$zmids,
  y = foo_Mstar$xmids,
  # z = foo_Mstar$contamination,
  z = compl_mat_smooth,
  xlab = "z",
  ylab = "MUV",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  xlim = range(foo_Mstar$zmids),
  ylim = range(foo_Mstar$xmids),
)
image.plot(
  legend.only = TRUE,
  zlim = range(foo_Mstar$contamination, na.rm = TRUE),
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  legend.lab = "Completeness"
)

compl_mat[is.na(compl_mat)] = 0
comp_interp <- list(
  x = zmids,
  y = Mstar_mids,
  z = compl_mat
)
# look up completeness for a single galaxy
interp.surface(comp_interp, cbind(10, 8))

magplot(
  MUV_mids, 
  interp.surface(comp_interp, cbind(2, MUV_mids))
)

cont_interp <- list(
  x = zmids,
  y = MUV_mids,
  z = cont
)
# look up contamination for a single galaxy
interp.surface(cont_interp, cbind(z = 2.1, mag = -20))

completeness_fns = list(
  "compl_fn" = function(z, M, comp_interp){interp.surface(comp_interp, cbind(z = z, mag = M))},
  "cont_fn" = function(z, M, cont_interp){interp.surface(cont_interp, cbind(z = z, mag = M))}
)

saveRDS(
  list(
    "fns" = completeness_fns,
    "compl" = compl,
    "cont" = cont,
    "zbins" = zbins,
    "Mbins" = MUV_bins, 
    "compl_interp" = comp_interp,
    "cont_interp" = cont_interp
  ), 
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_selection_completeness.rds")
)

fit_logistic <- function(x, y) {
  
  # Remove NAs
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  # Better starting values
  A_start  <- max(y)
  x0_start <- x[which.min(abs(y - A_start/2))]  # x where y ~ 50% of max
  k_start  <- 2 / diff(range(x))                 # rough steepness estimate
  
  start <- list(A = A_start, x0 = x0_start, k = k_start)
  
  fit <- tryCatch({
    nls(
      y ~ A / (1 + exp(-k * (x - x0))),
      start   = start,
      control = nls.control(maxiter = 1000, tol = 1e-5),
      algorithm = "port",
      lower = list(A = 0,   x0 = min(x), k = 0),
      upper = list(A = 1,   x0 = max(x), k = 100)
    )
  }, error = function(e) {
    message("nls failed, trying grid search for starting values...")
    
    # Grid search over starting values
    best <- NULL
    best_ss <- Inf
    
    for(A_try in seq(0.5, 1, by=0.1)){
      for(x0_try in seq(min(x), max(x), length.out=10)){
        for(k_try in c(0.5, 1, 2, 5, 10)){
          try({
            fit_try <- nls(
              y ~ A / (1 + exp(-k * (x - x0))),
              start   = list(A=A_try, x0=x0_try, k=k_try),
              control = nls.control(maxiter=500),
              algorithm = "port",
              lower = list(A=0, x0=min(x), k=0),
              upper = list(A=1, x0=max(x), k=100)
            )
            ss <- sum(residuals(fit_try)^2)
            if(ss < best_ss){ best <- fit_try; best_ss <- ss }
          }, silent=TRUE)
        }
      }
    }
    return(best)
  })
  
  return(fit)
}
fit <- fit_logistic(x = Mstar_mids, y = compl_mat[kk, ])

# Predict on fine grid
x_fine <- seq(-100, 100, 0.01)
y_pred  <- predict(fit, newdata = list(x = x_fine))
lines(x_fine, y_pred, col="red", lwd=2)

