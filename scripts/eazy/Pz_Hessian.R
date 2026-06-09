library(numDeriv)
library(Rfits)
library(data.table)
library(Highlander)
library(magicaxis)
library(foreach)
library(pracma)
library(rhdf5)

## compute Hessian matrices about MLE solutions in EAZY
eazy_dir = "~/Documents/Epochs/data/eazy/"
fname = "smiles"
fname_append = "_agn_blue_sfhz"

compute_eazy_stuff = function(eazy_dir, fname, fname_append){
  eazy_zPDFs = data.frame(fread(
    paste0(eazy_dir, fname, fname_append, "/", fname, "_eazy_zPDF", fname_append, ".csv")
  ))
  eazy = Rfits_read_table(
    paste0(eazy_dir, fname, fname_append, "/", fname, "_zphot", fname_append, ".out.zout.fits")
  )
  eazy_ntemplates = h5read(
    paste0(eazy_dir, fname, fname_append, "/", fname, "_zphot", fname_append, ".out.h5"),
    'fit/coeffs_best'
  )
  eazy_ntemplates = colSums(
    eazy_ntemplates > 0
  )
  
  hessian_Pz = function(kk){
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
        xx, mean = zbest, sd = sigma_z, log = TRUE 
      )
      magplot(
        zgrid,
        log(zPDF),
        type = "l", 
        # xlim = c(1.5, 3.5)
        # xlim = pmin(40, pmax(0.01, c(-1.5+zbest, 1.5+zbest))),
        # ylim = c( min(c(yy,zPDF[zPDF>0]))*0.5, max(yy,zPDF)*1.01 )
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
  
  zerrH = sapply(1:dim(eazy)[1], hessian_Pz)
  zbest_integrals = sapply(1:dim(eazy)[1], compute_peak_integrals)
  zbest_chisq_pval = sapply(1:dim(eazy)[1], compute_chisq_GOF)
  
  # now make Normal distribution with  sigma from the Hessian in the Laplace Approximation
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
    "zbestHessianSD" = zerrH,
    "zbestChi2" = eazy$z_ml_chi2,
    "zbestPzIntegral" = zbest_integrals,
    "zbestChi2Pval" = zbest_chisq_pval
  )
  
  fwrite(
    eazy_stats,
    paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_stats", fname_append, ".csv")
  )
  
  fwrite(
    zPDFs_normals,
    paste0("~/Documents/Epochs/data/eazy/", fname, fname_append, "/", fname, "_eazy_zPDF_normal", fname_append, ".csv")
  )
  
  out = list(
    "zPDFs_normals" = zPDFs_normals,
    "eazy_stats" = eazy_stats
  )
  return(out)
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
}

agn_blue_sfhz = compute_eazy_stuff(
  eazy_dir, 
  "smiles",
  "_agn_blue_sfhz"
)

# idx = zbest_integrals >= 0.6
# maghist(
#   eazy$z_ml[idx], log = "y", breaks = 100
# )
# magplot(
#   eazy$z_ml[idx],
#   eazy$z500[idx]
# )
# 
# pdf("~/Desktop/eazy_smiles_galaxies_Pz.pdf", height = 6, width = 8)
# hessian_Pz_plot(1)
# hessian_Pz_plot(6)
# hessian_Pz_plot(711)
# hessian_Pz_plot(2001)
# hessian_Pz_plot(3003)
# dev.off()
# 
# hessian_Pz_plot(3048)
