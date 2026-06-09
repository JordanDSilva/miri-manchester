library(ProFound)
library(ProSpect)
library(ProSpectUtils)
library(ProFit)
library(Rfits)
library(celestial)
library(magicaxis)
library(data.table)
library(pracma)
library(MASS)
library(foreach)
library(doParallel)
library(matrixStats)
library(nloptr)
library(Highlander)
library(furrr)
library(purrr)
library(stringr)
library(MASS)

fname = "smiles"
main_dir = "~/Documents/Epochs/"
load(
  paste0(
    main_dir, "data/prospect/prospect_mock_data_for_fitting.rda"
  )
)

pc_err = 0.1
Nsample = 10000

jaguar_mock = Rfits_read_table("~/Downloads/JADES_SF_mock_r1_v1.2.fits.gz")
# jaguar_Q_mock = Rfits_read_table("~/Downloads/JADES_Q_mock_r1_v1.2.fits.gz")

F444W_jaguar_distribution = density(
  -2.5*log10(jaguar_mock$NRC_F444W_fnu*1e-9)+8.9
)
F444W_jaguar_distribution_max = max(F444W_jaguar_distribution$y)
F444W_jaguar_rejection_sampler = function(F444W){
  check = approx(
    F444W_jaguar_distribution$x,
    F444W_jaguar_distribution$y,
    xout = F444W
  )$y / F444W_jaguar_distribution_max
  
  if(is.na(check)){
    return(0)
  }else{
    return(check)
  }
}

test_phot_cat = data.frame(fread(
  "~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_mock_eazy_photometry.csv"
))
test_cat = data.frame(fread(
  "~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_mock_cat.csv"
))
test_F444W_mag = -2.5*log10(as.numeric(test_phot_cat$F444W_true_flux))+8.9

down_sample = sample(
  x = test_F444W_mag,
  prob = approx(F444W_jaguar_distribution$x, F444W_jaguar_distribution$y, xout = test_F444W_mag, yleft = 0, yright = 0)$y,
  size = 10000
)
sample_idx = which(test_F444W_mag %in% down_sample)
control_idx = sample(1:100000, size = 10000)

magplot(
  test_cat$z[sample_idx],
  log10(test_cat$StellarMass[sample_idx]),
  pch = "."
  # ylim = c(-5, 5)
)
points(
  test_cat$z[control_idx],
  log10(test_cat$StellarMass[control_idx]),
  col = "red", 
  pch ="."
)

logged = c(
  FALSE,
  TRUE, FALSE, TRUE, FALSE,
  TRUE, TRUE, FALSE, FALSE,
  TRUE,
  TRUE, FALSE, TRUE
)
intervals = list(
  lo = c(
    0.01,
    -3, -2, log10(0.3), -0.5,
    -2.5, -5, 0, 0, 
    -4,
    35, 0.001, -1
  ), 
  hi = c(
    40, 
    4, 13.46728, 2, 1, 
    1.5, 1, 4, 4, 
    -1.3,
    49, 89.990, 1
  )
)

prospect_samples <- vector("list", Nsample)

N_accepted = 0
while(N_accepted < Nsample){
  
  parm <- runif(
    n   = length(intervals$lo),
    min = rep(intervals$lo),
    max = rep(intervals$hi)
  )
  
  parm[7] = runif(n = 1, min = intervals$lo[7], max = parm[6]) ## tau_screen < tau_birth
  parm[logged] = 10^parm[logged]
  
  Dale_M2L_func_updated = function(alpha_SF) {
    ProSpect::Dale_M2L_variableDTH_func(
      alpha_SF = alpha_SF, 
      qPAH_VSG = ProSpectUtils::shivaei24_qPAHZ(parm[10]) 
    )*(0.0073/RR14_BPL(parm[10]))
  }
  
  magemax = celestial::cosdistUniAgeAtz(z = parm[1], ref = '737')
  agemax = magemax * 1e9
  Zagemax = magemax
  
  mpeak_parm = runif(1, min = -2, max = magemax)
  parm[3] = mpeak_parm
  
  phot = ProSpectSED(
    z = parm[1],
    mSFR = parm[2],
    mpeak = parm[3],
    mperiod = parm[4],
    mskew = parm[5],
    tau_birth = parm[6],
    tau_screen = parm[7],
    alpha_SF_birth = parm[8],
    alpha_SF_screen = parm[9],
    Zfinal = parm[10],
    AGNlum = parm[11],
    AGNan = parm[12],
    AGNta = parm[13],
    
    agemax = agemax,
    magemax = magemax,
    Zagemax = Zagemax,
    
    AGN = Fritz,
    speclib = BC03hr,
    Dale = Dale_NormTot,
    Dale_M2L_func = Dale_M2L_func_updated,
    emission = TRUE,
    filtout = filter_transforms,
    filters = filters_info$filter_names,
    returnall = TRUE,
    
    SFH = SFHfunc,
    SMstar = TRUE,
    massfunc = massfunc_snorm_trunc,
    Z = Zfunc_massmap_lin,
    IGMabsorb = 'Inoue14',
    
    mtrunc = 2,
    
    AGNrm = 60,
    AGNbe = -0.5,
    AGNal = 4.0,
    AGNct = 100,
    
    ref='737',
    H0 = 70,
    OmegaM = 0.3,
    OmegaL = 0.7, 
  )
  
  mag_test = -2.5*log10(as.numeric(phot$Photom[filters_info$filters == "F444W"])) + 8.9
  mass_test = log10( as.numeric(phot$Stars$SMstar['TotSMstar']) )
  
  # Rejection step lives in the sampler, not the density function
  u <- runif(1)
  accept <- u < F444W_jaguar_rejection_sampler(mag_test) & mass_test > 6
  
  if(accept){
    print(N_accepted)
    N_accepted = N_accepted + 1
    
    phot_fesc = ProSpectSED(
      z = parm[1],
      mSFR = parm[2],
      mpeak = parm[3],
      mperiod = parm[4],
      mskew = parm[5],
      tau_birth = parm[6],
      tau_screen = parm[7],
      alpha_SF_birth = parm[8],
      alpha_SF_screen = parm[9],
      Zfinal = parm[10],
      AGNlum = parm[11],
      AGNan = parm[12],
      AGNta = parm[13],
      
      escape_frac = 1.0,
      
      agemax = agemax,
      magemax = magemax,
      Zagemax = Zagemax,
      
      AGN = Fritz,
      speclib = BC03hr,
      Dale = Dale_NormTot,
      Dale_M2L_func = Dale_M2L_func_updated,
      emission = TRUE,
      filtout = filter_transforms,
      filters = filters_info$filter_names,
      returnall = TRUE,
      
      SFH = SFHfunc,
      SMstar = TRUE,
      massfunc = massfunc_snorm_trunc,
      Z = Zfunc_massmap_lin,
      IGMabsorb = 'Inoue14',
      
      mtrunc = 2,
      
      AGNrm = 60,
      AGNbe = -0.5,
      AGNal = 4.0,
      AGNct = 100,
      
      ref='737',
      H0 = 70,
      OmegaM = 0.3,
      OmegaL = 0.7, 
    )
    
    UV1500 = UVLum1500(
      out = list("SEDout" = phot)
    )
    UV1500AGN = UVLumAGN1500(
      out = list("SEDout" = phot)
    )
    nions = NionEmissivity(
      out = list("SEDout" = phot_fesc)
    )
    xi_ion = ionisingEfficiency(
      out = list("SEDout" = phot_fesc)
    )
    nionsAGN = NionEmissivityAGN(
      out = list("SEDout" = phot_fesc)
    )
    xi_ionAGN = ionisingEfficiencyAGN(
      out = list("SEDout" = phot_fesc)
    )
    
    massfunc_args = c(
      mSFR = parm[2],
      mpeak = parm[3],
      mperiod = parm[4],
      mskew = parm[5],
      magemax = magemax,
      mtrunc = 2
    )
    SFR10 = do.call("integrate", c(list(f = massfunc_snorm_trunc, lower = 0, upper = 1e+07), massfunc_args))$value/1e+07
    
    df = list(
      'astro' = c(
        "z" = parm[1],
        "StellarMass" = as.numeric(phot$Stars$SMstar['TotSMstar']),
        "SFR10" = SFR10,
        "SFR100" = phot$Stars$SFRburst,
        "AGNlum" = parm[11],
        "UV1500" = UV1500,
        "UV1500AGN" = UV1500AGN,
        "DustMass" = as.numeric(phot$dustmass["total"]),
        "DustLum" = as.numeric(phot$dustlum["total"]),
        "Nions" = as.numeric(nions),
        "XiIon" = as.numeric(xi_ion),
        "NionsAGN" = as.numeric(nionsAGN),
        "XiIonAGN" = as.numeric(xi_ionAGN)
      ),
      parm = parm,
      phot = phot$Photom
    )
    
    prospect_samples[[N_accepted]] = df
  }
}

saveRDS(
  prospect_samples,
  paste0(main_dir, "/source_inject/", fname, "_jaguar/", fname, "_jaguar_mock_prospect_samples.RDS")
)
stop("Done")

prospect_samples = readRDS(
  paste0(main_dir, "/source_inject/", fname, "_jaguar/", fname, "_jaguar_mock_prospect_samples.RDS")
)

astro_cat = data.frame(do.call(rbind, lapply(prospect_samples, function(x){x$astro[1:length(x$astro)]})))

parm_cat = data.frame(do.call(rbind, lapply(prospect_samples, function(x){x$parm[1:length(x$parm)]})))
colnames(parm_cat) = c(
  'z',
  'mSFR','mpeak','mperiod','mskew',
  'tau_birth','tau_screen','alpha_SF_birth','alpha_SF_screen',
  'Zfinal',
  'AGNlum','AGNan','AGNta'
)

prospect_cat = cbind(astro_cat, parm_cat)
prospect_cat <- prospect_cat[, !duplicated(colnames(prospect_cat))]

frame_sigma = unlist(unname(frame_info[, paste0(filters_info$filters, "_2sigma_depth_32asec")]))/2.0

Nflux = length(filters_info$filters)

df_flux_raw = do.call(
  rbind, 
  lapply(prospect_samples, function(x)x$phot)
)
df_flux_err = data.frame(matrix(
  frame_sigma,
  ncol = Nflux,
  nrow = dim(df_flux_raw)[1],
  byrow = TRUE
))

df_flux_shuffle = matrix(rnorm(n = length(df_flux_raw), mean=as.matrix(df_flux_raw), sd=as.matrix(df_flux_err)), nrow(df_flux_raw))
df_flux = data.frame(df_flux_shuffle)
df_true_flux = data.frame(df_flux_raw)

flag_unrealistically_negative_fluxes = (df_flux + 4*df_flux_err) <= 0
df_flux[flag_unrealistically_negative_fluxes] = -999

df_flux_err = sqrt((df_flux_err)^2 + (pc_err * df_flux_err)^2)

names(df_flux) = paste0(filters_info$filters, "_flux")
names(df_true_flux) = paste0(filters_info$filters, "_true_flux")
names(df_flux_err) = paste0(filters_info$filters, "_flux_err")

# ## worse case scenario - only fit the bands with the most coverage...
bands_to_keep = paste0(c(
  "F435W", "F606W", "F775W", "F814W", "F850LP",
  "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W",
  "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W"
), "_flux")
bands_to_set_NA = names(df_flux) %notin% bands_to_keep
df_flux[,bands_to_set_NA] = -999

flux_df = cbind(
  ID = 1:Nsample, df_flux, df_flux_err, df_true_flux, zspec = prospect_cat$z
)

fwrite(
  flux_df, 
  paste0(main_dir, "/source_inject/", fname, "_jaguar/", fname, "_jaguar_mock_eazy_photometry.csv")
)
fwrite(
  prospect_cat, 
  paste0(main_dir, "/source_inject/", fname, "_jaguar/", fname, "_jaguar_mock_cat.csv")
)

eazy_translate = read.table(paste0("~/Documents/Epochs/data/eazy/", fname, "/eazy_", fname, "_zphot.translate"))
eazy_translate$V1 = str_remove_all(str_remove_all(eazy_translate$V1, "_app_0.32"), "_scaled")
fwrite(
  eazy_translate, 
  paste0(main_dir, "/source_inject/eazy_", fname, "_uniform_zphot_mock.translate"),
  sep = " ", 
  col.names = FALSE,
  row.names = FALSE
)
