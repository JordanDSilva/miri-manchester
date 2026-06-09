library(numDeriv)
library(Rfits)
library(data.table)
library(magicaxis)
library(foreach)
library(pracma)
library(Highlander)
library(rhdf5)
library(celestial)

compute_sigmaNMAD <- function(zspec, zphot, zphot_lim = 0) {
  idx = zphot > zphot_lim
  deltaz   <- zphot[idx] - zspec[idx]
  
  sigmaNMAD <- 1.48 * median( abs( (deltaz - median(deltaz)) / (1 + zspec[idx]) ) )
  return(sigmaNMAD)
}
compute_outliers <- function(zspec, zphot, thresh = 0.15, zphot_lim = 0) {
  idx = zphot > zphot_lim
  deltaz     <- abs(zphot[idx] - zspec[idx])
  Noutlier   <- sum(deltaz > thresh * (1 + zspec[idx]))
  outlier_frac <- Noutlier / length(zspec[idx])
  
  cat("Noutlier:", Noutlier, "\n")
  cat("N total: ", length(zspec[idx]), "\n")
  
  return(list("frac" = outlier_frac, "idx" = deltaz > thresh * (1 + zspec[idx])))
}

## compute Hessian matrices about MLE solutions in EAZY
fname = "smiles"
eazy_zPDFs = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_mock_photometry_eazy_zPDF.csv")
))
eazy = Rfits_read_table(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot.out.zout.fits")
)
eazy_ntemplates = h5read(
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_zphot.out.h5"),
  'fit/coeffs_best'
)
eazy_ntemplates = as.numeric(colSums(
  eazy_ntemplates > 0
))

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))

hessian_Pz = function(kk){
  print(kk)
  doplot = FALSE
  zbest = eazy$z_ml[kk]
  
  zgrid = eazy_zPDFs$zgrid
  zPDF = pmax(1e-99, eazy_zPDFs[[paste0('gal',kk)]])
  
  window_idx <- zPDF > 0

  if(sum(window_idx) <= 1){
    sigma_z = 40
  }else if(zbest == -1){
    sigma_z = 40
  }else{
    log_zPDF = splinefun(
      x = c(zgrid[window_idx]),
      y = c(log(zPDF[window_idx]))
    )
    
    LL = function(parm, Data){
      log_zPDF(parm)
    }
    opt = Highlander(
      parm = zbest,
      likefunc = LL,
      Data = list(),
      lower = pmin(40, pmax(0.01, c(-0.5+zbest, 0.5+zbest)))[1],
      upper = pmin(40, pmax(0.01, c(-0.5+zbest, 0.5+zbest)))[2],
      liketype = "max", 
      optim_iters = 2,
      Niters = c(100,100),
      NfinalMCMC = 100
    )
    
    H = hessian(log_zPDF, as.numeric(opt$parm))
    if(H > 0){
      H = hessian(log_zPDF, as.numeric(zbest))
    }
    
    sigma_z <- as.numeric( sqrt(-1 / H) )
    if(is.nan(sigma_z) | is.na(sigma_z) | is.infinite(sigma_z)){
      sigma_z = 40
    }
  }
  
  if(doplot){
    xx = seq(0, 50, 0.01)
    yy = dnorm(
      xx, mean = zbest, sd = sigma_z, log = FALSE 
    )
    magplot(
      zgrid,
      zPDF,
      type = "l", 
      xlim = pmin(40, pmax(0.01, c(-1.5+zbest, 1.5+zbest))),
      ylim = c( min(c(yy,zPDF[zPDF>1e-99]))*0.5, max(yy,zPDF)*1.1 )
    )
    lines(
      xx, 
      yy, 
      col = "blue"
    )
    abline(
      v = zbest, 
      lty = 2
    )
  }
  
  return(sigma_z)
}
hessian_Pz_plot = function(kk){
  zgrid = eazy_zPDFs$zgrid
  zPDF = eazy_zPDFs[[paste0('gal',kk)]]
  
  # print(zgrid[zPDF > 0])
  
  xx = seq(0.01, 50, 0.001)
  fit = dnorm(
    xx, 
    mean = eazy$z_ml[kk],
    sd = zerrH[kk], 
    log = FALSE
  )
  magplot(
    xx,
    fit,
    type = "l",
    # log = "y",
    xlim = c(0.01,c(range(zgrid[zPDF > 0]))[2]),
    # xlim = c(range(zgrid[zPDF > 0])),
    ylim = c(pmax(1e-99, range(c(zPDF,fit)))),
    xlab = "Redshift",
    ylab = "P(z)", 
    grid = FALSE
  )
  lines(
    zgrid,
    (zPDF),
    col = "blue",
    lwd = 2
  )
  abline(
    v = eazy$z_ml[kk],
    lty = 2
  )
  abline(
    v = astro_cat$z[kk],
    lty = 3
  )
  legend(
    x = 'topleft',
    c(
      paste0("z true = ", round(astro_cat$z[kk], 3)),
      paste0("z phot = ", round(eazy$z_ml[kk], 3)),
      paste0("AGNlum = ", round(log10(astro_cat$AGNlum[kk]), 3)),
      paste0("Stellar mass = ", round(log10(astro_cat$StellarMass[kk]), 3)),
      paste0("Det mag = ", round(astro_cat$detmag[kk], 3)), 
      "Eazy P(z)",
      "Normal on z best"
    ), 
    pch = c(NA, NA, NA, NA, NA, NA),
    lty = c(NA, NA, NA, NA, 1, 1),
    col = c(NA, NA, NA, NA, "blue", "black")
  )
}
compute_peak_integrals = function(kk){
  zgrid = eazy_zPDFs$zgrid
  zPDF = eazy_zPDFs[[paste0('gal',kk)]]
  
  zbest = eazy$z_ml[kk]
  
  idx = zgrid >= 0.9*zbest & zgrid <= 1.1*zbest
  trapz(
    zgrid[idx],
    zPDF[idx]
  )
}
compute_chisq_GOF = function(kk){
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

zerrH = sapply(1:dim(astro_cat)[1], hessian_Pz)
zbest_integrals = sapply(1:dim(astro_cat)[1], compute_peak_integrals)
zbest_chisq_pval = sapply(1:dim(astro_cat)[1], compute_chisq_GOF)

# which(is.na(zerrH))
# which.max(zerrH)
# max(zerrH)

zPDFs_normals = data.frame(
  foreach(i = 1:dim(eazy)[1], .combine = cbind) %do% {
    if(zerrH[i] == 40){
      Pz = rep(1, length(eazy_zPDFs$zgrid)) ## uniform, no prior
    }else{
      Pz = dnorm(
        x = eazy_zPDFs$zgrid, mean = eazy$z_ml[i], sd = zerrH[i]
      )
    }
    return(Pz)
  }
)
names(zPDFs_normals) = paste0("gal", 1:dim(eazy)[1])
zPDFs_normals$zgrid = eazy_zPDFs$zgrid

eazy_stats = data.frame(
  # "zbestHessianSD" = zerrH,
  "zbestChi2" = eazy$z_ml_chi2,
  "zbestPzIntegral" = zbest_integrals,
  "zbestChi2Pval" = zbest_chisq_pval
)

fwrite(
  eazy_stats,
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_stats.csv")
)

fwrite(
  zPDFs_normals,
  paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_zPDF_normal.csv")
)

eazy_stats = data.frame(fread(paste0("~/Documents/Epochs/source_inject/", fname, "/eazy/", fname, "_eazy_stats.csv")))
idx = astro_cat$detected == 1

magplot(
  eazy$z_ml,
  astro_cat$StellarMass, 
  z = idx,
  log = "y"
)


pdf("~/Desktop/eazy_mock_galaxies_Pz.pdf", height = 6, width = 8)
hessian_Pz_plot(8)
hessian_Pz_plot(4)
hessian_Pz_plot(19)
hessian_Pz_plot(2995)
hessian_Pz_plot(4156)
dev.off()

pdf("~/Documents/Epochs/plots/completeness/zeazy_ztrue_with_checks.pdf", width = 12, height = 7)
par(mfrow = c(1,2), oma = rep(1.5,4), mar = rep(2.5, 4))
magplot(
  eazy$z_ml,
  astro_cat$z, 
  xlim = c(0, 40),
  ylim = c(0, 40),
  xlab = "z phot",
  ylab = "z true",
  pch = "."
)
legend(
  x = "topleft",
  legend = c(
    paste0("N = ", dim(astro_cat)[1]),
    paste0("sigmaNMAD = ", round(compute_sigmaNMAD(astro_cat$z, eazy$z_ml),3)),
    paste0("Outlier frac = ", round(compute_outliers(astro_cat$z, eazy$z_ml)$frac,3))
  )
)

magplot(
  eazy$z_ml[idx],
  astro_cat$z[idx], 
  pmax(6, log10(astro_cat$StellarMass[idx])),
  xlim = c(0, 40),
  ylim = c(0, 40),
  xlab = "z phot",
  ylab = "z true",
  pch = 16, 
  main = "With photoz checks", 
)
legend(
  x = "topleft",
  legend = c(
    paste0("N = ", sum(idx)),
    paste0("sigmaNMAD = ", round(compute_sigmaNMAD(astro_cat$z[idx], eazy$z_ml[idx]),3)),
    paste0("Outlier frac = ", round(compute_outliers(astro_cat$z[idx], eazy$z_ml[idx])$frac,3))
  )
)
dev.off()

plot(
  density(astro_cat$z, bw = 0.1), 
  xlim = c(0, 40), 
  ylim = c(0, 0.2)
)
lines(
  density(astro_cat$z[idx], bw = 0.1),     
  col = "blue"
)
plot(
  density(log10(astro_cat$StellarMass), bw = 0.1), 
  xlim = c(5, 12), ylim = c(0, 1.1)
)
lines(
  density(log10(astro_cat$StellarMass[idx]), bw = 0.1),    
  col = "blue"
)
plot(
  density(log10(astro_cat$SFR100), bw = 0.1),
  xlim = c(-5, 5), ylim = c(0, 0.8)
)
lines(
  density(log10(astro_cat$SFR100[idx]), bw = 0.1),     
  col = "blue"
)
plot(
  density(log10(astro_cat$AGNlum), bw = 0.1),
  xlim = c(35, 50), ylim = c(0, 0.15)
)
lines(
  density(
    log10(astro_cat$AGNlum[idx]), bw = 0.1
  ),     col = "blue"
)

