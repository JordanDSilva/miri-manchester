library(Rfits)
library(data.table)
library(magicaxis)
library(rhdf5)
library(pracma)
library(celestial)
library(fields)
library(matrixStats)

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_uniform_detect"

compute_sigmaNMAD <- function(zspec, zphot, zphot_lim = 0) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz   <- zphot[idx] - zspec[idx]
  
  sigmaNMAD <- 1.48 * median( abs( (deltaz - median(deltaz)) / (1 + zspec[idx]) ) )
  return(sigmaNMAD)
}
compute_outliers <- function(zspec, zphot, thresh = 0.15, zphot_lim = -2) {
  idx = zphot > zphot_lim & !is.na(zphot)
  deltaz     <- abs(zphot[idx] - zspec[idx])
  
  outlier_idx = deltaz/(1 + zspec[idx]) > thresh
  Noutlier   <- sum(outlier_idx)
  outlier_frac <- Noutlier / length(zspec[idx])
  
  cat("Noutlier:", Noutlier, "\n")
  cat("N total: ", length(zspec[idx]), "\n")
  
  return(list(outlier_frac=outlier_frac, outlier_idx=outlier_idx))
}

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
compute_chisq_GOF = function(eazy, eazy_ntemplates, kk){
  if(eazy$z_ml[kk] == -1 | eazy_ntemplates[kk] == 0){
    pval = 0.0
  }else{
    pval = pchisq(
      q = eazy$z_phot_chi2[kk],
      df = as.numeric(eazy$nusefilt[kk] - eazy_ntemplates[kk] - 1),
      lower.tail = FALSE
    )
  }
}

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))
# phot_cat = data.frame(fread(
#   paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv")
# ))

load_eazy = function(template_fname){
  eazy_zPDFs = data.frame(fread(
    paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_mock_photometry_eazy_zPDF", template_fname, ".csv")
  ))
  eazy = Rfits_read_table(
    paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot", template_fname, ".out.zout.fits")
  )
  eazy_ntemplates = h5read(
    paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot", template_fname, ".out.h5"),
    'fit/coeffs_best'
  )
  eazy_ntemplates = as.numeric(colSums(
    eazy_ntemplates > 0
  ))
  out = list(
    "eazy_zPDFs" = eazy_zPDFs,
    "eazy" = eazy,
    "eazy_ntemplates" = eazy_ntemplates
  )
  return(out)
}

# larson_default = load_eazy("_larson")
# sfhz_carnall = load_eazy("_sfhz_carnall")
agn_blue_sfhz = load_eazy("_agn_blue_sfhz")

# larson_eazy_chi2 = sapply(1:dim(astro_cat)[1], function(x)compute_chisq_GOF(larson_default$eazy, larson_default$eazy_ntemplates, x))
# larson_eazy_integrals = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(larson_default$eazy_zPDFs, larson_default$eazy, x))
# 
# sfhz_carnall_eazy_chi2 = sapply(1:dim(astro_cat)[1], function(x)compute_chisq_GOF(sfhz_carnall$eazy, sfhz_carnall$eazy_ntemplates, x))
# sfhz_carnall_eazy_integrals = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(sfhz_carnall$eazy_zPDFs, sfhz_carnall$eazy, x))

agn_blue_sfhz_eazy_chi2 = sapply(1:dim(astro_cat)[1], function(x)compute_chisq_GOF(agn_blue_sfhz$eazy, agn_blue_sfhz$eazy_ntemplates, x))
agn_blue_sfhz_eazy_integrals = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(agn_blue_sfhz$eazy_zPDFs, agn_blue_sfhz$eazy, x))

zspec_zphot_plot = function(zphot, zspec, zphot_lim = 0, select_idx = 1:length(zphot)){
  magplot(
    zphot,
    zspec, 
    z = select_idx,
    xlim = c(0,40),
    ylim = c(0,40),
    xlab = "zphot",
    ylab = "zspec"
  )
  legend(
    x = "topright",
    legend = c(
      paste0("sigmaNMAD = ", round(compute_sigmaNMAD(zspec, zphot, zphot_lim = zphot_lim), 3)),
      paste0("Outliers = ", round(compute_outliers(zspec, zphot, zphot_lim = zphot_lim)$outlier_frac, 3))
    )
  )
}

# zspec_zphot_plot(
#   zphot = larson_default$eazy$z_ml,
#   zspec = larson_default$eazy$z_spec,
#   # select_idx = log10(astro_cat$StellarMass)
# )
# outliers = compute_outliers(zspec = larson_default$eazy$z_spec,zphot = larson_default$eazy$z_ml)
# plot(density(log10(astro_cat$StellarMass[outliers$outlier_idx])))
# plot(density(log10(astro_cat$StellarMass)))
# 
# zspec_zphot_plot(
#   zphot = sfhz_carnall$eazy$z_ml,
#   zspec = sfhz_carnall$eazy$z_spec
# )
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml[selection_idx == 0],
  zspec = agn_blue_sfhz$eazy$z_spec[selection_idx == 0], 
  
)

# zspec_zphot_plot(
#   zphot = larson_default$eazy$z_ml[larson_eazy_integrals >= 0.6],
#   zspec = larson_default$eazy$z_spec[larson_eazy_integrals >= 0.6]
# )
# zspec_zphot_plot(
#   zphot = sfhz_carnall$eazy$z_ml[sfhz_carnall_eazy_integrals >= 0.6],
#   zspec = sfhz_carnall$eazy$z_spec[sfhz_carnall_eazy_integrals >= 0.6],
#   # select_idx = sfhz_carnall_eazy_integrals >= 0.6
# )
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml[agn_blue_sfhz_eazy_integrals >= 0.6],
  zspec = agn_blue_sfhz$eazy$z_spec[agn_blue_sfhz_eazy_integrals >= 0.6]
)

# magplot(
#   astro_cat$z,
#   log10(astro_cat$StellarMass),
#   z = larson_eazy_integrals >= 0.6, 
#   ylim = c(4, 15), 
#   pch = "."
# )
# magplot(
#   astro_cat$z,
#   log10(astro_cat$StellarMass),
#   z = sfhz_carnall_eazy_integrals >= 0.6
# )
magplot(
  # log(1+astro_cat$z),
  log(1 + agn_blue_sfhz$eazy$z_ml),
  log10(astro_cat$StellarMass),
  z = agn_blue_sfhz_eazy_integrals >= 0.6, 
  pch = ".",
  ylim = c(4, 15)
)

calc_selection_3d = function(
    zbins, zmids, ztrue, 
    xtrue, xbins, xmids,
    ytrue, ybins, ymids,
    selection_idx,
    ztest = NULL,
    xtest = NULL,
    ytest = NULL
){
  Ndim = c(length(zmids), length(xmids), length(ymids))
  completeness_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  purity_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  
  if(is.null(ztest)){
    ztest = rnorm(length(ztrue), ztrue, 0.06)
  }
  if(is.null(xtest)){
    xtest = rnorm(length(xtrue), xtrue, 0.3)
    # xtest = xtrue
  }
  if(is.null(ytest)){
    ytest = rnorm(length(ytrue), ytrue, 0.3)
    # ytest = ytrue
  }
  
  for(i in seq_along(zmids)){
    for(j in seq_along(xmids)){
      for(k in seq_along(ymids)){
        id_array_true = which(
          (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) & 
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1])
        )
        
        id_array_sel = which(
          (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) &
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1]) &
            (selection_idx)
        )
        
        id_array_det_sel = which(
          (ztest >= zbins[i]) & (ztest < zbins[i+1]) &
            (xtest >= xbins[j]) & (xtest < xbins[j+1]) &
            (ytest >= ybins[k]) & (ytest < ybins[k+1]) &
            (selection_idx)
        )
        
        id_array_det_sel_exact = intersect(
          id_array_sel, id_array_det_sel
        ) # Find exact matches in z-x bin
        
        # hh_true = maghist(
        #   x = ytrue[id_array_true],
        #   breaks = ybins[(k):(k+1)],
        #   plot = FALSE
        # )
        # hh_det_sel_exact = maghist(
        #   x = ytest[id_array_det_sel_exact], # the exact matches of ytest, purity will be how many are exact given how many we select
        #   breaks = ybins[(k):(k+1)],
        #   plot = FALSE
        # )
        # hh_det_sel = maghist(
        #   x = ytest[id_array_det_sel],
        #   breaks = ybins[(k):(k+1)],
        #   plot = FALSE
        # )
        
        det_sel_counts = length(id_array_det_sel)
        true_counts = length(id_array_true)
        true_sel_counts = length(id_array_sel)
        det_sel_exact_counts = length(id_array_det_sel_exact)
        
        compl = det_sel_counts/true_counts
        purity = det_sel_exact_counts/det_sel_counts
        
        compl[is.infinite(compl)] = NA
        purity[is.infinite(purity)] = NA
        
        completeness_array[i,j,k] = compl
        purity_array[i,j,k] = purity
      }
    }
    
    compl_mat = as.matrix(completeness_array[i,,])
    cont_mat = 1 - as.matrix(purity_array[i,,])
    
    # par(mfrow = c(1,2))
    # image(
    #   x = xmids,
    #   y = ymids,
    #   z = compl_mat,
    #   xlab = "Mstar",
    #   ylab = "SFR",
    #   col = hcl.colors(100, "YlOrRd", rev = TRUE),
    #   xlim = range(xbins),
    #   ylim = range(ybins), 
    # )
    # image.plot(
    #   legend.only = TRUE,
    #   zlim = c(0, 1),
    #   col = hcl.colors(100, "YlOrRd", rev = TRUE),
    #   legend.lab = "Completeness"
    # )
    # legend(
    #   x = "topleft",
    #   legend = c(
    #     paste0("z = ", round(zmids[i],3)), 
    #     paste0("Min = ", round(min(compl_mat, na.rm = TRUE), 3)),
    #     paste0("Max = ", round(max(compl_mat, na.rm = TRUE), 3))
    #   )
    # )
    # 
    # image(
    #   x = xmids,
    #   y = ymids,
    #   z = cont_mat,
    #   xlab = "Mstar",
    #   ylab = "SFR",
    #   col = hcl.colors(100, "YlOrRd", rev = TRUE),
    #   xlim = range(xbins),
    #   ylim = range(ybins), 
    # )
    # image.plot(
    #   legend.only = TRUE,
    #   zlim = c(0, 1),
    #   col = hcl.colors(100, "YlOrRd", rev = TRUE),
    #   legend.lab = "Contamination"
    # )
    # legend(
    #   x = "topleft",
    #   legend = c(
    #     paste0("z = ", round(zmids[i],3)), 
    #     paste0("Min = ", round(min(cont_mat, na.rm = TRUE), 3)),
    #     paste0("Max = ", round(max(cont_mat, na.rm = TRUE), 3))
    #   )
    # )
  }
  
  out = list(
    "compl" = completeness_array,
    "cont" = 1 - purity_array
  )
  return( out )
}

## Pretend that we have done SED fitting
# ytest = log10(rnorm(100000, astro_cat$SFR10, 0.0001*astro_cat$SFR10))
# xtest = log10(rnorm(100000, astro_cat$StellarMass, 0.0001*astro_cat$StellarMass))
# ztest = agn_blue_sfhz$eazy$z_ml

# ytest = log10(astro_cat$SFR10)
ytest = astro_cat$UV1500
xtest = log10(astro_cat$StellarMass)
ztest = astro_cat$z

lbt_bins = c(seq(8.0, 13.0, 0.5), 13.1)
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

selection_idx = agn_blue_sfhz_eazy_integrals >= 0.6 & agn_blue_sfhz$eazy$z_ml > 0
foo = calc_selection_3d(
  zbins = zbins,
  zmids = zmids,
  ztrue = astro_cat$z,
  ztest = agn_blue_sfhz$eazy$z_ml,
  
  xbins = Mstar_bins,
  xmids = Mstar_mids,
  xtrue = log10(astro_cat$StellarMass),
  xtest = xtest,
  
  ybins = MUV_bins,
  ymids = MUV_mids,
  # ytrue = log10(astro_cat$SFR10),
  ytrue = astro_cat$UV1500,
  ytest = ytest,
  
  selection_idx = selection_idx
)

saveRDS(
  foo, 
  "~/Desktop/mstar_sfr_compl_test.rds"
)

kk = 1
zidx_test = ztest >= zbins[kk] & ztest < zbins[kk + 1] & selection_idx == 1
ztemp = data.frame(
  "StellarMass" = xtest[zidx_test],
  "SFR10" = ytest[zidx_test]
)
zidx_true = astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk + 1]
zidx_true_sel = astro_cat$z >= zbins[kk] & astro_cat$z < zbins[kk + 1] & selection_idx == 1

Ndim = c(length(Mstar_mids), length(SFR_mids))
temp = array(data = rep(0, prod(Ndim)), dim = Ndim)
for(i in seq_along(Mstar_mids)){
  for(j in seq_along(SFR_mids)){
    idx = ztemp$StellarMass >= Mstar_bins[i] & ztemp$StellarMass < Mstar_bins[i+1] & ztemp$SFR10 >= SFR_bins[j] & ztemp$SFR10 < SFR_bins[j+1]
    
    contamination = foo$cont[kk, i, j]
    if((contamination > 0.6) | is.nan(contamination)){
      temp[i,j] = NA
    }else{
      temp[i,j] = sum(idx / foo$compl[kk, i, j])
    }
    
  }
}

par(mfrow = c(1,2))
hh = maghist(
  log10(astro_cat$StellarMass[zidx_true]), 
  breaks = Mstar_bins, 
  log = "y", 
  plot = FALSE
) # ground truth
hh_true_select = maghist(
  log10(astro_cat$StellarMass[zidx_true_sel]), 
  breaks = Mstar_bins, 
  log = "y", 
  plot = FALSE
) # ground truth
hh_raw_select = maghist(
  ztemp$StellarMass, 
  breaks = Mstar_bins, 
  log = "y", 
  plot = FALSE
)
magplot(
  hh$mids, 
  hh$counts, 
  type = "s", 
  log = "y",
  ylim = c(1, 1500), 
  lwd = 2
)
lines(Mstar_mids, rowSums(temp, na.rm = TRUE), type = "s", pch = 16, col = "blue", lwd = 3)
lines(Mstar_mids, hh_raw_select$counts, type = 's', col = "red")
lines(Mstar_mids, hh_true_select$counts, type = 's', col = "darkgreen")

hh = maghist(
  log10(astro_cat$SFR10[zidx_true]), 
  breaks = SFR_bins, 
  log = "y", 
  plot = FALSE
) # ground truth
hh_raw_select = maghist(
  ztemp$SFR10, 
  breaks = SFR_bins, 
  log = "y", 
  plot = FALSE
)
hh_true_select = maghist(
  log10(astro_cat$SFR10[zidx_true_sel]), 
  breaks = SFR_bins, 
  log = "y", 
  plot = FALSE
) # ground truth
magplot(
  hh$mids, 
  hh$counts, 
  type = "s", 
  log = "y",
  ylim = c(1, 1500), 
  lwd = 2
)
lines(SFR_mids, colSums(temp, na.rm = TRUE), type = 's', col = "blue")
lines(SFR_mids, hh_raw_select$counts, type = 's', col = "red")
lines(SFR_mids, hh_true_select$counts, type = 's', col = "darkgreen")

magplot(
  Mstar_mids,
  rowMedians(foo$compl[kk, , ], na.rm = TRUE)
)
magplot(
  SFR_mids,
  colMedians(foo$compl[kk, , ], na.rm = TRUE)
)

magplot(
  Mstar_mids,
  rowMedians(foo$cont[kk, , ], na.rm = TRUE)
)
magplot(
  SFR_mids,
  colMedians(foo$cont[kk, , ], na.rm = TRUE)
)


image(
  x = Mstar_mids,
  y = MUV_mids,
  # z = matrix(pmin(1,foo$compl[kk,,]), nrow = length(Mstar_mids), ncol = length(SFR_mids)),
  z = foo$compl[kk,,],
  xlab = "Mstar",
  ylab = "SFR",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  ylim = c(-30, -10),
  xlim = c(4, 13)
)
image.plot(
  legend.only = TRUE,
  # zlim = c(0, 1),
  zlim = c(
  min(foo$compl[kk,,], na.rm = TRUE),
  max(foo$compl[kk,,], na.rm = TRUE)
  ),
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  legend.lab = "Completeness"
)

image(
  x = Mstar_mids,
  y = SFR_mids,
  z = foo$cont[kk,,],
  xlab = "Mstar",
  ylab = "SFR",
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  ylim = c(-5, 5),
  xlim = c(4, 13)
)
image.plot(
  legend.only = TRUE,
  # zlim = c(0, 1),
  zlim = c(
    min(foo$cont[kk,,], na.rm = TRUE), 
    max(foo$cont[kk,,], na.rm = TRUE)
  ),
  col = hcl.colors(100, "YlOrRd", rev = TRUE),
  legend.lab = "Contamination"
)

inter_obj = list(
  "x" = Mstar_bins, 
  "y" = SFR_bins,
  "z" = foo$compl[kk,,]
)
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
completeness_func = make_completeness_func(inter_obj, Mstar_mids, SFR_mids)