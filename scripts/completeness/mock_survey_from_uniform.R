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
    main_dir, "data/prospect_mock_data_for_fitting.rda"
  )
)

phot_cat = Rfits_read(
  "~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits"
)

# catalogues with 1e5 prospect sources sampled uniformly
mock_uniform_phot = data.frame(fread(
  "~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_mock_eazy_photometry.csv"
))
mock_uniform_astro = data.frame(fread(
  "~/Documents/Epochs/source_inject/smiles_uniform/smiles_uniform_mock_cat.csv"
))
mock_uniform_astro$ID = 1:dim(mock_uniform_astro)[1]

test_F770W_mag = -2.5*log10(as.numeric(mock_uniform_phot$F770W_true_flux))+8.9

# Make true F770W distribution for synthetic survey
F770W_distribution = density(
  -2.5*log10(phot_cat$Phot32$F770W_flux_app_0.32[phot_cat$Phot32$F770W_flux_app_0.32 > 0])+8.9, 
  bw = "SJ"
)

robust_idx = which(mock_uniform_astro$StellarMass > 1e3 & mock_uniform_astro$StellarMass < 1e12 & test_F770W_mag >= 10 & test_F770W_mag <= 40)

# use this for the main survey
sample_idx = sample(
  x = robust_idx,
  # prob = approx(F770W_distribution$x, F770W_distribution$y, xout = test_F770W_mag[robust_idx], yleft = 0, yright = 0)$y,
  size = 100000
)

# use this for training an emulator
training_idx = sample(
  x = robust_idx[robust_idx %notin% sample_idx],
  # prob = approx(F770W_distribution$x, F770W_distribution$y, xout = test_F770W_mag[robust_idx], yleft = 0, yright = 0)$y,
  size = 300000
)

fwrite(
  mock_uniform_phot[sample_idx, ], 
  paste0(main_dir, "/source_inject/", fname, "_survey/", fname, "_survey_mock_eazy_photometry.csv")
)
fwrite(
  mock_uniform_astro[sample_idx, ], 
  paste0(main_dir, "/source_inject/", fname, "_survey/", fname, "_survey_mock_cat.csv")
)

eazy_translate = read.table(paste0("~/Documents/Epochs/data/eazy/eazy_", fname, "_zphot.translate"))
eazy_translate$V1 = str_remove_all(str_remove_all(eazy_translate$V1, "_app_0.32"), "_scaled")
fwrite(
  eazy_translate, 
  paste0(main_dir, "/source_inject/eazy_", fname, "_survey_zphot_mock.translate"),
  sep = " ", 
  col.names = FALSE,
  row.names = FALSE
)

fwrite(
  mock_uniform_phot[training_idx, ], 
  paste0(main_dir, "/source_inject/", fname, "_survey/", fname, "_survey_mock_training_photometry.csv")
)
fwrite(
  mock_uniform_astro[training_idx, ], 
  paste0(main_dir, "/source_inject/", fname, "_survey/", fname, "_survey_mock_training_cat.csv")
)