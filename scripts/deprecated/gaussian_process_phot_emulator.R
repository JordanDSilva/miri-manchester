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
library(hetGP)

## Load in the 100 000 survey galaxies with magnitudes following the same distribution I observe (synthetic observation)
main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_survey"
fname_append = "_agn_blue_sfhz"

fittype = "photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth"

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv")
))

SED_matrix = as.matrix(
  phot_cat[, grep(".+rest_flux$", names(phot_cat), value = TRUE)]
)
n_wave <- ncol(SED_matrix)
gp_models <- vector("list", n_wave)

train_parms = astro_cat[, c("z", "mSFR", "mpeak", "mperiod", "mskew", "tau_birth", "tau_screen", "alpha_SF_birth", "alpha_SF_screen", "Zfinal", "AGNlum", "AGNan", "AGNta" )]
logged = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
train_parms[,logged] = log10(train_parms[,logged])

sample_idx = sample(1:100000, size = 10000)
for(w in 1:n_wave){
  cat("Training wavelength bin", w, "of", n_wave, "\r")
  parm_grid = train_parms[sample_idx, ]
  mag = -2.5*log10(SED_matrix[sample_idx, w])+8.9
  fit_idx = mag >= 5 & mag <= 40
  
  gp_models[[w]] <- mleHetGP(
    X = as.matrix(parm_grid[fit_idx, ]),
    Z = mag[fit_idx]
  )
}

predict_bands <- function(parm_new) {
  sapply(gp_models, function(gp) {
    predict(gp, x = matrix(as.numeric(parm_new), nrow=1))$mean
  })
}

test_idx = sample(1:100000, size = 1000)
foo = do.call(rbind, lapply(test_idx, function(x)predict_bands(train_parms[x,])))

magplot(
  foo[,15],
  -2.5*log10(SED_matrix[test_idx,15])+8.9, 
  # xlim = c(10, 40),
  # ylim = c(10, 40)
)
abline(0, 1, col ="red")


library(ranger)  # fast random forest

rf_models <- lapply(1:n_bands, function(b) {
  cat("Band", b, "\n")
  df <- as.data.frame(train_parms)
  df$y <- train_fluxes[, b]
  ranger(y ~ ., data = df, num.trees = 500)
})

predict_fluxes <- function(parm_new) {
  sapply(rf_models, function(rf) {
    predict(rf, data = as.data.frame(t(parm_new)))$predictions
  })
}