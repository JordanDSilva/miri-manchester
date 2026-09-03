library(data.table)
library(magicaxis)

foo = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_lrd_astro_cat.csv"))
bar = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd_astro_cat.csv"))
poo = data.frame(fread("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))

lrd_idx = which(!is.na(foo$AGNlumQ50))

# ks test foo and poo
ks_samples = vector(length = 1000)
median_difference_samples = vector(length = 1000)
for(i in 1:1000){
  x = rnorm(
    n = length(lrd_idx),
    mean = log10(foo$AGNlumQ50[lrd_idx]),
    sd = 0.5 * ( log10(foo$AGNlumQ84[lrd_idx]) - log10(foo$AGNlumQ16[lrd_idx]) )
  )
  
  y = rnorm(
    n = length(lrd_idx),
    mean = log10(poo$AGNlumQ50[lrd_idx]),
    sd = 0.5 * ( log10(poo$AGNlumQ84[lrd_idx]) - log10(poo$AGNlumQ16[lrd_idx]) )
  )
  
  ks_ = ks.test(
    x,
    y
  )
  ks_samples[i] = ks_$p.value
  
  median_difference_samples[i] = median(
    x - y
  )
}
print(
  summary(ks_samples)
)
print(
  summary(median_difference_samples)
)

magplot(
  (foo$AGNlumQ50),
  (poo$AGNlumQ50),
  xlim = 10^c(40, 47),
  ylim = 10^c(35, 47),
  log = "xy"
)
magerr(
  (foo$AGNlumQ50),
  (poo$AGNlumQ50),
  xlo =   (foo$AGNlum_lo),
  ylo = (foo$AGNlum_lo),
  xhi =   (poo$AGNlum_hi),
  yhi = (poo$AGNlum_hi),
)
abline(0,1)

# ks_samples = vector(length = 1000)
# for(i in 1:1000){
#   x = rnorm(
#     n = length(lrd_idx),
#     mean = log10(foo$AGNlumQ50[lrd_idx]),
#     sd = 0.5 * ( log10(foo$AGNlumQ84[lrd_idx]) - log10(foo$AGNlumQ16[lrd_idx]) )
#   )
#   
#   y = rnorm(
#     n = length(lrd_idx),
#     mean = log10(bar$AGNlumQ50[lrd_idx]),
#     sd = 0.5 * ( log10(bar$AGNlumQ84[lrd_idx]) - log10(bar$AGNlumQ16[lrd_idx]) )
#   )
#   
#   ks_ = ks.test(
#     x,
#     y
#   )
#   ks_samples[i] = ks_$p.value
# }
# 
# ks.test(
#   log10(foo$AGNlumQ50[lrd_idx]),
#   log10(poo$AGNlumQ50[lrd_idx])
# )


magplot(
  (foo$AGNlumQ50)[lrd_idx],
  (bar$AGNlumQ50)[lrd_idx],
  xlim = 10^c(40, 47),
  ylim = 10^c(40, 47),
  log = "xy"
)
magerr(
  (foo$AGNlumQ50)[lrd_idx],
  (bar$AGNlumQ50)[lrd_idx],
  xlo =   (foo$AGNlum_lo)[lrd_idx],
  ylo = (bar$AGNlum_lo)[lrd_idx],
  xhi =   (foo$AGNlum_hi)[lrd_idx],
  yhi = (bar$AGNlum_hi)[lrd_idx],
)
abline(0,1)

magplot(
  (foo$StellarMassQ50),
  (bar$StellarMassQ50),
  xlim = 10^c(4, 11),
  ylim = 10^c(4, 11),
  log = "xy"
)
magerr(
  (foo$StellarMassQ50),
  (bar$StellarMassQ50),
  xlo =   (foo$StellarMass_lo),
  ylo = (bar$StellarMass_lo),
  xhi =   (foo$StellarMass_hi),
  yhi = (bar$StellarMass_hi)
)
abline(0,1)

akins25 = data.frame(fread(file = "~/Documents/Epochs/data/literature/AGN-seds/akins25_stack.csv"))
casey24 = read.table("~/Documents/Epochs/data/literature/AGN-seds/casey24_lrd_predicted_sed_agn.txt")

for(i in foo$CATID[!is.na(foo$zQ50)]){
  print(i)
  fritz = readRDS(
    paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_lrd/galID_", i, ".rds")
  )
  lrdmbb = readRDS(
    paste0("~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd/galID_", i, ".rds")
  )
  fritz_sed = readRDS(
    paste0("~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy_lrd/galID_", i, ".rds")
  )
  lrdmbb_sed = readRDS(
    paste0("~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd/galID_", i, ".rds")
  )
  
  pdf(paste0("~/Documents/Epochs/plots/rough-plots/lrd_seds/gal_", i, ".pdf"), width = 8, height = 8)
  magplot(
    fritz$bestfit$Data$flux$cenwave/1e4,
    fritz$bestfit$Data$flux$flux,
    log = "xy",
    xlim = c(0.1, 3e4),
    ylim = c(1e-10, 1e-1)
  )
  magerr(
    fritz$bestfit$Data$flux$cenwave/1e4,
    fritz$bestfit$Data$flux$flux,
    ylo = fritz$bestfit$Data$flux$fluxerr
  )
  lines(
    fritz_sed$obswave/1e4,
    fritz$bestfit$SEDout$FinalFlux$flux,
    col = "cornflowerblue"
  )
  magerr(
    fritz_sed$obswave/1e4,
    fritz_sed$fluxQ50,
    ylo = fritz_sed$fluxQ50 - fritz_sed$fluxQ16,
    yhi = fritz_sed$fluxQ84 - fritz_sed$fluxQ50,
    poly = TRUE,
    col = adjustcolor("cornflowerblue", 0.1),
    border = NA
  )
  
  lines(
    lrdmbb_sed$obswave/1e4,
    lrdmbb$bestfit$SEDout$FinalFlux$flux,
    col = "blue"
  )
  magerr(
    fritz_sed$obswave/1e4,
    lrdmbb_sed$fluxQ50,
    ylo = (lrdmbb_sed$fluxQ50 - lrdmbb_sed$fluxQ16) + 1e-323,
    yhi = (lrdmbb_sed$fluxQ84 - lrdmbb_sed$fluxQ50) + 1e-323,
    poly = TRUE,
    border = NA,
    col = adjustcolor("blue", 0.1)
  )
  
  # agn_rest = ProSpect::Lum2Flux(
  #   wave = fritz_sed$AGNWave,
  #   lum = fritz_sed$AGNQ50, 
  #   z = fritz$astro$zQ50, 
  #   ref = "737"
  # )
  # lines(
  #   agn_rest$wave/1e4,
  #   ProSpect::CGS2Jansky(
  #     ProSpect::convert_wave2freq(
  #       agn_rest$flux,
  #       agn_rest$wave
  #     )
  #   ),
  #   col = "cornflowerblue", lty = 2
  # )
  # 
  # agn_rest = ProSpect::Lum2Flux(
  #   wave = lrdmbb_sed$AGNWave,
  #   lum = lrdmbb_sed$AGNQ50, 
  #   z = lrdmbb$astro$zQ50,
  #   ref = "737"
  # )
  # lines(
  #   agn_rest$wave/1e4,
  #   ProSpect::CGS2Jansky(
  #     ProSpect::convert_wave2freq(
  #       agn_rest$flux,
  #       agn_rest$wave
  #     )
  #   ),
  #   col = "blue",
  #   lty = 2
  # )
  dev.off()
}

foo = readRDS(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_lrd/galID_688.rds"
)
foo_sed = readRDS(
  "~/Documents/Epochs/data/prospect_seds/photoZ_bc03_snorm_fritz_norm_eazy_lrd/galID_688.rds"
)
bar = readRDS(
  "~/Documents/Epochs/data/prospect/photoZ_bc03_snorm_LRDMBB_norm_eazy_lrd/galID_688.rds"
)
# 
#
