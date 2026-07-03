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
library(ranger)  # fast random forest

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

# sample_idx = sample(1:100000, size = 10000)
sample_idx = 1:100000

# Build training dataframe
# df_multi <- as.data.frame(train_parms[sample_idx, ])
df_multi = train_parms[sample_idx, ]
# Add all bands as output columns
mag_matrix <- -2.5*log10(SED_matrix[sample_idx, ]) + 8.9
mag_matrix[!is.finite(mag_matrix)] <- NA
colnames(mag_matrix) <- paste0("band_", 1:n_wave)

# Step 2: PCA on the SED matrix
pca <- prcomp(mag_matrix, scale. = TRUE, center = TRUE)

# How many components needed?
var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
n_comp <- which(var_explained > 0.999)[1]
cat("Need", n_comp, "components for 99% variance\n")

plot(var_explained, type="l", xlab="n components", ylab="variance explained")
abline(h=0.99, col="red", lty=2)

# Step 3: Train one RF per PCA component
scores <- pca$x[, 1:n_comp]  # n_samples x n_comp 

rf_pca <- lapply(1:n_comp, function(k) {
  df <- as.data.frame(train_parms)
  df$y <- scores[, k]
  ranger(y ~ ., data = df, num.trees = 1000)
})

# Step 4: Predict function
predict_fluxes_pca <- function(parm_df) {
  newdata <- as.data.frame(parm_df)
  
  # Predict each PCA component
  scores_new <- sapply(rf_pca, function(rf) {
    predict(rf, data = newdata)$predictions
  })
  # scores_new is n_samples x n_comp
  
  # Reconstruct mags from PCA
  reconstructed <- scores_new %*% t(pca$rotation[, 1:n_comp])
  
  # Unscale and uncenter
  reconstructed <- sweep(reconstructed, 2, pca$scale, "*")
  reconstructed <- sweep(reconstructed, 2, pca$center, "+")
  
  return(reconstructed)  # n_samples x n_bands
}

save(
  n_comp,
  pca,
  rf_pca, 
  predict_fluxes_pca,
  file = "~/Documents/Epochs/data/prospect_photom_emulator/prospect_photom_predictor_RF_PCA.rda"
)

# Much faster - one call per band instead of n_samples x n_bands calls
# Test
start_ = Sys.time()
test_idx <- sample(sample_idx, size = 10000)
foo <- predict_fluxes_pca(train_parms[test_idx, ])
print(
  Sys.time() - start_
)

magplot(
  foo[,1],
  mag_matrix[,1][test_idx], 
  pch = ".",
  xlim = c(10,40),
  ylim = c(10,40)
)
abline(0,1)

png("~/Documents/Epochs/data/prospect_photom_emulator/mag_mag.png", width = 10, height = 10, units = "in", res = 100)
par(
  mfrow = c(8,4), 
  oma = rep(1,4),
  mar = c(2.5, 2.5, 1.5, 1.5)
)
for(kk in 1:31){
  magplot(
    -2.5*log10(SED_matrix[test_idx,kk])+8.9, 
    foo[,kk],
    xlim = c(10, 40),
    ylim = c(10,40),
    xlab = "True mag",
    ylab = "Predicted mag", 
    pch = "."
  )
  abline(0,1,col = "red")
  legend(
    x = "topleft",
    legend = c(
      paste0(grep(".+true_flux$", names(phot_cat), value = TRUE)[kk])
    )
  )
}
dev.off()

png("~/Documents/Epochs/data/prospect_photom_emulator/mag_hist.png", width = 10, height = 10, units = "in", res = 100)
par(
  mfrow = c(8,4), 
  oma = rep(0.5,4),
  mar = c(2.5, 2.5, 1.5, 1.5)
)
breaks = seq(-10, 10, 0.5)
for(kk in 1:31){
  
  resid <- foo[,kk] - (-2.5*log10(SED_matrix[test_idx,kk])+8.9)
  
  maghist(
    resid,
    breaks = breaks,
    log = "y",
    xlim = c(-10, 10),
    ylim = c(1, 1000)
  )
  
  q16  <- quantile(resid, 0.16, na.rm=TRUE)
  q50  <- quantile(resid, 0.50, na.rm=TRUE)
  q84  <- quantile(resid, 0.84, na.rm=TRUE)
  
  legend("topright",
         legend = c(
           paste0(grep(".+true_flux$", names(phot_cat), value = TRUE)[kk]),
           sprintf("median = %.3f", q50),
           sprintf("q16 = %.3f",    q16),
           sprintf("q84 = %.3f",    q84)
         ),
         bty = "n"
  )
}
dev.off()

