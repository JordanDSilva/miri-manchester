library(data.table)
library(magicaxis)

spec_match_agn_idx = fread("~/Documents/JWST-MIRI-AGN/data/prospect/match_better_fit_AGN.txt")[[2]]

prospect_agn = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_astro_cat.csv"))
prospect_noagn = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_astro_cat.csv"))

nircam_agn = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_fritz_norm_eazy_nircam_astro_cat.csv"))
nircam_noagn = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/prospect/photoZ_bc03_snorm_norm_eazy_nircam_astro_cat.csv"))

# magplot(
#   prospect_agn$zQ50[spec_match_agn_idx],
#   nircam_agn$zQ50
# )

magplot(
  prospect_agn$AGN_Dale_frac_Q50[spec_match_agn_idx],
  nircam_agn$AGN_Dale_frac_Q50,
  xlab = "HST+NIRCam+MIRI AGN fractional contribution to the SED",
  ylab = "HST+NIRCam AGN fractional contribution to the SED",
  xlim = c(0, 1),
  ylim = c(0, 1)
)
abline(0,1)

magplot(
  prospect_agn$StellarMassQ50[spec_match_agn_idx],
  nircam_agn$StellarMassQ50,
  xlab = "HST+NIRCam+MIRI stellar mass",
  ylab = "HST+NIRCam stellar mass",
  xlim = c(1e5, 1e12),
  ylim = c(1e5, 1e12),
  log = "xy"
)
magerr(
  prospect_agn$StellarMassQ50[spec_match_agn_idx],
  nircam_agn$StellarMassQ50,
  xlo = prospect_agn$StellarMass_lo[spec_match_agn_idx],
  xhi = prospect_agn$StellarMass_hi[spec_match_agn_idx],
  ylo = nircam_agn$StellarMass_lo,
  yhi = nircam_agn$StellarMass_hi,
)
abline(0,1)
print(
  median(
    log10(nircam_agn$StellarMassQ50/prospect_agn$StellarMassQ50[spec_match_agn_idx])
  )
)

magplot(
  prospect_agn$SFR10Q50[spec_match_agn_idx],
  nircam_agn$SFR10Q50,
  xlab = "HST+NIRCam+MIRI SFR",
  ylab = "HST+NIRCam SFR",
  xlim = c(1e-4, 1e2),
  ylim = c(1e-4, 1e2),
  log = "xy"
)
magerr(
  prospect_agn$SFR10Q50[spec_match_agn_idx],
  nircam_agn$SFR10Q50,
  xlo = prospect_agn$SFR10_lo[spec_match_agn_idx],
  xhi = prospect_agn$SFR10_hi[spec_match_agn_idx],
  ylo = nircam_agn$SFR10_lo,
  yhi = nircam_agn$SFR10_hi,
)
abline(0,1)
print(
  median(
    log10(nircam_agn$SFR10Q50/prospect_agn$SFR10Q50[spec_match_agn_idx])
  )
)


magplot(
  prospect_agn$AGNlumQ50[spec_match_agn_idx],
  nircam_agn$AGNlumQ50,
  xlab = "HST+NIRCam+MIRI LAGN",
  ylab = "HST+NIRCam LAGN",
  xlim = c(1e35, 1e47),
  ylim = c(1e35, 1e47),
  log = "xy"
)
magerr(
  prospect_agn$AGNlumQ50[spec_match_agn_idx],
  nircam_agn$AGNlumQ50,
  xlo = prospect_agn$AGNlum_lo[spec_match_agn_idx],
  xhi = prospect_agn$AGNlum_hi[spec_match_agn_idx],
  ylo = nircam_agn$AGNlum_lo,
  yhi = nircam_agn$AGNlum_hi,
)
abline(0,1)

median(mean(c(prospect_agn$SFR10_lo[spec_match_agn_idx],prospect_agn$SFR10_hi[spec_match_agn_idx]))) / median(mean(c(nircam_agn$SFR10_lo,nircam_agn$SFR10_hi)))

median(mean(c(prospect_agn$DustLumTotal_lo[spec_match_agn_idx],prospect_agn$DustLumTotal_hi[spec_match_agn_idx]))) / median(mean(c(nircam_agn$DustLumTotal_lo,nircam_agn$DustLumTotal_hi)))
# 
# median(mean(c(prospect_agn$SFR10_lo[spec_match_agn_idx],prospect_agn$SFR10_hi[spec_match_agn_idx])))
# median(mean(c(nircam_agn$SFR10_lo,nircam_agn$SFR10_hi)))
# 
# median(mean(c(prospect_agn$StellarMass_lo[spec_match_agn_idx],prospect_agn$StellarMass_hi[spec_match_agn_idx])))
# median(mean(c(nircam_agn$StellarMass_lo,nircam_agn$StellarMass_hi)))
# 
# median(mean(c(prospect_agn$AGNlum_lo[spec_match_agn_idx],prospect_agn$AGNlum_hi[spec_match_agn_idx])))
# median(mean(c(nircam_agn$AGNlum_lo,nircam_agn$AGNlum_hi)))
# 
