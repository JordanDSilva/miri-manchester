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

jaguar_z_distribution = density(jaguar_mock$redshift)
z_samples = sample(jaguar_mock$redshift, Nsample)

kde2 <- kde2d(
  jaguar_mock$redshift, 
  jaguar_mock$mStar, 
  n = 200,
  lims = c( 
    range(jaguar_mock$redshift), 
    range(jaguar_mock$mStar) 
  )
)

# Precompute once
cond_max <- apply(kde2$z, 1, function(row) {
  r <- pmax(row, 0)
  max(r / sum(r))
})

# Function just returns the interpolated conditional probability
conditional_mstar_density <- function(z_val, mstar_val, kde2) {
  
  ix <- which.min(abs(kde2$x - z_val))
  
  cond_density <- kde2$z[ix, ]
  cond_density <- pmax(cond_density, 0)
  cond_density <- cond_density / sum(cond_density)
  
  prob <- approx(
    x    = kde2$y,
    y    = cond_density,
    xout = mstar_val
  )$y
  
  if (is.na(prob)) return(0)  # out of mStar range → always reject
  return(prob)
}

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

registerDoParallel(cores = 8)
prospect_samples = foreach(i = 1:Nsample) %dopar% {
  if(i %% 1000 == 0){
    cat("Accepted: ", i)
    flush.console()
  }

  count = 0
  rejection_sample_mstar_check = FALSE
  
  while(!rejection_sample_mstar_check){
    
    parm <- runif(
      n   = length(intervals$lo),
      min = rep(intervals$lo),
      max = rep(intervals$hi)
    )
    parm[1] = z_samples[i] # reset to z samples from JAGUAR catalogue
    
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
    
    mass_test = log10( as.numeric(phot$Stars$SMstar['TotSMstar']) )

    # Rejection step lives in the sampler, not the density function
    u <- runif(1)
    ix <- which.min(abs(kde2$x - parm[1]))
    accept <- u < conditional_mstar_density(parm[1], mass_test, kde2) / cond_max[ix]

    if(accept){
      rejection_sample_mstar_check = TRUE
      
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
    }else{
      count = count + 1
    }
  }
  
  return(df)
}
stopImplicitCluster()

saveRDS(
  prospect_samples,
  paste0(main_dir, "/source_inject/", fname, "_jaguar/", fname, "_jaguar_mock_prospect_samples.RDS")
)
# stop("Done")

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
