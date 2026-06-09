library(Rfits)
library(data.table)
library(magicaxis)
library(rhdf5)
library(pracma)
library(celestial)
library(fields)

main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_survey"

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

larson_default = load_eazy("_larson")
sfhz_carnall = load_eazy("_sfhz_carnall")
agn_blue_sfhz = load_eazy("_agn_blue_sfhz")

larson_eazy_chi2 = sapply(1:dim(astro_cat)[1], function(x)compute_chisq_GOF(larson_default$eazy, larson_default$eazy_ntemplates, x))
larson_eazy_integrals = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(larson_default$eazy_zPDFs, larson_default$eazy, x))

sfhz_carnall_eazy_chi2 = sapply(1:dim(astro_cat)[1], function(x)compute_chisq_GOF(sfhz_carnall$eazy, sfhz_carnall$eazy_ntemplates, x))
sfhz_carnall_eazy_integrals = sapply(1:dim(astro_cat)[1], function(x)compute_peak_integrals(sfhz_carnall$eazy_zPDFs, sfhz_carnall$eazy, x))

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

zspec_zphot_plot(
  zphot = larson_default$eazy$z_ml,
  zspec = larson_default$eazy$z_spec,
  # select_idx = log10(astro_cat$StellarMass)
)
outliers = compute_outliers(zspec = larson_default$eazy$z_spec,zphot = larson_default$eazy$z_ml)
plot(density(log10(astro_cat$StellarMass[outliers$outlier_idx])))
plot(density(log10(astro_cat$StellarMass)))

zspec_zphot_plot(
  zphot = sfhz_carnall$eazy$z_ml,
  zspec = sfhz_carnall$eazy$z_spec
)
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml,
  zspec = agn_blue_sfhz$eazy$z_spec
)

zspec_zphot_plot(
  zphot = larson_default$eazy$z_ml[larson_eazy_integrals >= 0.6],
  zspec = larson_default$eazy$z_spec[larson_eazy_integrals >= 0.6]
)
zspec_zphot_plot(
  zphot = sfhz_carnall$eazy$z_ml[sfhz_carnall_eazy_integrals >= 0.6],
  zspec = sfhz_carnall$eazy$z_spec[sfhz_carnall_eazy_integrals >= 0.6],
  # select_idx = sfhz_carnall_eazy_integrals >= 0.6
)
zspec_zphot_plot(
  zphot = agn_blue_sfhz$eazy$z_ml[agn_blue_sfhz_eazy_integrals >= 0.6],
  zspec = agn_blue_sfhz$eazy$z_spec[agn_blue_sfhz_eazy_integrals >= 0.6],
  # select_idx = agn_sfhz_carnall_eazy_integrals >= 0.6
)

magplot(
  astro_cat$z,
  log10(astro_cat$StellarMass),
  z = larson_eazy_integrals >= 0.6, 
  ylim = c(4, 15)
)
magplot(
  astro_cat$z,
  log10(astro_cat$StellarMass),
  z = sfhz_carnall_eazy_integrals >= 0.6
)
magplot(
  log(1+astro_cat$z),
  log10(astro_cat$StellarMass),
  z = agn_blue_sfhz_eazy_integrals >= 0.6, 
  ylim = c(4, 15)
)

calc_selection_3d = function(
    zbins, zmids, ztrue, ztest, 
    xtrue, xbins, xmids,
    ytrue, ybins, ymids,
    selection_idx,
    xerr = 0.3,
    yerr = 0.3
){
  Ndim = c(length(zmids), length(xmids), length(ymids))
  completeness_array = array(data = rep(0, prod(Ndim)), dim = Ndim)
  
  for(i in seq_along(zmids)){
    true_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1])
    true_selected_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (selection_idx)
    
    # ztest = rnorm(length(ztrue), ztrue, ztrue * 0.1) # use eazy zml
    xtest = rnorm(length(xtrue), xtrue, xerr)
    ytest = rnorm(length(ytrue), ytrue, yerr)
    
    det_selected_idx = (ztest >= zbins[i]) & (ztest < zbins[i+1]) & (selection_idx)
    
    xlevels = seq_along(xmids)
    ylevels = seq_along(ymids)
    
    xbin_true = cut(xtrue[true_idx], breaks = xbins, labels = FALSE)
    ybin_true = cut(ytrue[true_idx], breaks = ybins, labels = FALSE)
    
    xbin_sel  = cut(xtrue[true_selected_idx], breaks = xbins, labels = FALSE)
    ybin_sel  = cut(ytrue[true_selected_idx], breaks = ybins, labels = FALSE)
    
    xbin_det_sel  = cut(xtest[det_selected_idx & true_selected_idx], breaks = xbins, labels = FALSE)
    ybin_det_sel  = cut(ytest[det_selected_idx & true_selected_idx], breaks = ybins, labels = FALSE)
    
    # Force same factor levels so tables always have same dimensions
    n_true = table(
      factor(xbin_true, levels = xlevels),
      factor(ybin_true, levels = ylevels)
    )
    n_sel = table(
      factor(xbin_sel, levels = xlevels),
      factor(ybin_sel, levels = ylevels)
    )
    n_det_sel = table(
      factor(xbin_det_sel, levels = xlevels),
      factor(ybin_det_sel, levels = ylevels)
    )
    
    compl_mat =  (n_det_sel) / n_sel
    compl_mat[n_true > 0 & n_true < 5] = NA
    compl_mat[is.infinite(compl_mat) | is.na(compl_mat)] = NA
    # compl_mat[compl_mat > 1] = 1
    
    image(
      x = xmids,
      y = ymids,
      z = compl_mat,
      xlab = "Mstar",
      ylab = "SFR",
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      xlim = range(xbins),
      ylim = range(ybins), 
    )
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      # zlim = c(min(foo$completeness[kk,,], na.rm = TRUE), max(foo$completeness[kk,,], na.rm = TRUE)),
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      legend.lab = "Completeness"
    )
    # points(
    #   log10(astro_cat$StellarMass[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1]]),
    #   log10(astro_cat$SFR100[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1]]),
    #   col = adjustcolor("black", 0.4),
    #   pch = 16
    # )
    # points(
    #   log10(astro_cat$StellarMass[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1] & selection_idx == 1]),
    #   log10(astro_cat$SFR100[astro_cat$z >= zbins[i] & astro_cat$z < zbins[i+1] & selection_idx == 1]),
    #   col = adjustcolor("magenta", 0.4), 
    #   pch = 16
    # )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3)), 
        paste0("Min = ", round(min(compl_mat, na.rm = TRUE), 3)),
        paste0("Max = ", round(max(compl_mat, na.rm = TRUE), 3))
      )
    )
    completeness_array[i,,] = compl_mat
  }
  return( completeness_array )
}

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
    # xtest = rnorm(length(xtrue), xtrue, 0.3)
    xtest = xtrue
  }
  if(is.null(ytest)){
    # ytest = rnorm(length(ytrue), ytrue, 0.3)
    ytest = ytrue
  }
  
  for(i in seq_along(zmids)){
    for(j in seq_along(xmids)){
      for(k in seq_along(ymids)){
        
        true_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1])
        true_selected_idx = (ztrue >= zbins[i]) & (ztrue < zbins[i+1]) & (selection_idx)
        det_selected_idx = (ztest >= zbins[i]) & (ztest < zbins[i+1]) & (selection_idx)
        
        id_array_true = which(
          true_idx &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) &
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1])
        )
        
        id_array_sel = which(
          true_selected_idx &
            (xtrue >= xbins[j]) & (xtrue < xbins[j+1]) &
            (ytrue >= ybins[k]) & (ytrue < ybins[k+1])
        )
        
        id_array_det_sel = which(
          det_selected_idx &
            (xtest >= xbins[j]) & (xtest < xbins[j+1]) &
            (ytest >= ybins[k]) & (ytest < ybins[k+1])
        )
        
        id_array_det_sel_exact = intersect(
          id_array_sel, id_array_det_sel
        )
        
        completeness_array[i,j,k] = length(id_array_sel)/length(id_array_true)
        purity_array[i,j,k] = length(id_array_det_sel_exact)/length(id_array_det_sel)
      }
    }
    
    cont_mat = 1 - as.matrix(purity_array[i,,])
    # cont_mat[cont_mat > 0.5] = 2
    compl_mat = as.matrix(completeness_array[i,,])

    par(mfrow = c(1,2))
    image(
      x = xmids,
      y = ymids,
      z = cont_mat,
      xlab = "Mstar",
      ylab = "SFR",
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      xlim = range(xbins),
      ylim = range(ybins), 
    )
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      legend.lab = "Contamination"
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3)), 
        paste0("Min = ", round(min(cont_mat, na.rm = TRUE), 3)),
        paste0("Max = ", round(max(cont_mat, na.rm = TRUE), 3))
      )
    )
    image(
      x = xmids,
      y = ymids,
      z = compl_mat,
      xlab = "Mstar",
      ylab = "SFR",
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      xlim = range(xbins),
      ylim = range(ybins), 
    )
    image.plot(
      legend.only = TRUE,
      zlim = c(0, 1),
      col = hcl.colors(100, "YlOrRd", rev = TRUE),
      legend.lab = "Completeness"
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("z = ", round(zmids[i],3)), 
        paste0("Min = ", round(min(compl_mat, na.rm = TRUE), 3)),
        paste0("Max = ", round(max(compl_mat, na.rm = TRUE), 3))
      )
    )
  }
  return( completeness_array )
}

lbt_bins = c(seq(9.0, 13.0, 0.5), 13.1)
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

Mstar_bins = seq(1, 15, 0.5)
Mstar_mids = Mstar_bins[1:(length(Mstar_bins)-1)] + diff(Mstar_bins)/2.0

foo = calc_selection_3d(
  zbins = zbins,
  zmids = zmids,
  ztrue = astro_cat$z,
  ybins = SFR_bins,
  ymids = SFR_mids,
  ytrue = log10(astro_cat$SFR10),
  xbins = Mstar_bins,
  xmids = Mstar_mids,
  xtrue = log10(astro_cat$StellarMass),
  selection_idx = larson_eazy_integrals >= 0.6, 
  ztest = larson_default$eazy$z_ml
)
