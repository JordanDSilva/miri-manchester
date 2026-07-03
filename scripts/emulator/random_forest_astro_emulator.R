library(data.table)
library(magicaxis)
library(foreach)
library(pracma)
library(celestial)
library(ranger)  # fast random forest

## Load in the 100 000 survey galaxies with magnitudes following the same distribution I observe (synthetic observation)
main_dir = "~/Documents/Epochs/source_inject/"
fname = "smiles_survey"
fname_append = "_agn_blue_sfhz"

smiles_depths = data.frame(fread("~/Documents/Epochs/data/depths/smiles_profound_frame_stats.csv"))

fittype = "photoZ_bc03_snorm_fritz_norm_eazy_no_tau_birth"

astro_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_cat.csv")
))
phot_cat = data.frame(fread(
  paste0("~/Documents/Epochs/source_inject/", fname, "/", fname, "_mock_eazy_photometry.csv")
))

band_names = c(
  "F435W", "F606W", "F775W", "F814W", "F850LP",
  "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W",
  "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W"
)
bands_to_keep = paste0(c(
  "F435W", "F606W", "F775W", "F814W", "F850LP",
  "F090W", "F115W", "F150W", "F200W", "F277W", "F335M", "F356W", "F410M", "F444W",
  "F560W", "F770W", "F1000W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W"
), "_flux")

train_SED_matrix = as.matrix(
  phot_cat[, bands_to_keep]
) # can have negative fluxes 
n_wave <- ncol(train_SED_matrix)

band_noise = unname(as.numeric(smiles_depths[, paste0(band_names, "_5sigma_depth_32asec")]))

mag_matrix = -2.5*log10(exp(1)) * ( asinh(train_SED_matrix / (2*band_noise)) + log(band_noise/3631) )
colnames(mag_matrix) <- paste0("band_", 1:n_wave)

color_matrix <- sapply(2:n_wave, function(i) mag_matrix[, i-1] - mag_matrix[, i])
colnames(color_matrix) <- paste0("color_", 1:(n_wave-1))

feature_matrix <- cbind("z" = astro_cat$z, mag_matrix, color_matrix)

train_parms = astro_cat[, c("StellarMass", "SFR10", "AGNlum", "UV1500", "UV1500AGN")]
logged = c(TRUE, TRUE, TRUE, FALSE, FALSE)
train_parms[,logged] = log10(train_parms[,logged])

## --- assemble full training data.frame ---
full_data <- data.frame(feature_matrix, train_parms)

# Drop rows with non-finite values (bad fluxes/targets)
ok <- complete.cases(full_data) & apply(full_data, 1, function(x) all(is.finite(x)))
full_data <- full_data[ok, ]

cat(sprintf("Kept %d / %d rows after removing non-finite cases\n", nrow(full_data), length(ok)))

## --- train/test split ---
set.seed(42)
n <- nrow(full_data)
train_idx <- sample(seq_len(n), size = round(0.8 * n))
train_set <- full_data[train_idx, ]
test_set  <- full_data[-train_idx, ]

feature_names <- colnames(feature_matrix)
target_names  <- colnames(train_parms)

## --- train one ranger forest per target ---
rf_models <- list()
predictions <- data.frame(matrix(NA, nrow = nrow(test_set), ncol = length(target_names)))
colnames(predictions) <- target_names

for (targ in target_names) {
  
  form <- as.formula(paste(targ, "~", paste(feature_names, collapse = " + ")))
  
  rf_models[[targ]] <- ranger(
    formula = form,
    data = train_set,
    num.trees = 1000,
    mtry = max(1, floor(length(feature_names) / 3)),
    min.node.size = 5,
    importance = "permutation",
    num.threads = parallel::detectCores() - 1
  )
  
  pred <- predict(rf_models[[targ]], data = test_set)$predictions
  predictions[[targ]] <- pred
  
  resid <- pred - test_set[[targ]]
  rmse  <- sqrt(mean(resid^2))
  r2    <- 1 - sum(resid^2) / sum((test_set[[targ]] - mean(test_set[[targ]]))^2)
  
  cat(sprintf("%-12s  RMSE = %6.4f   R^2 = %6.4f\n", targ, rmse, r2))
}

saveRDS(
  list("models" = rf_models, "test_set" = test_set),
  "~/Documents/Epochs/data/prospect_photom_emulator/prospect_astro_predictor_RF.rda"
)

par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
for (targ in target_names) {
  magplot(
    test_set[[targ]], predictions[[targ]],
    xlab = paste("True", targ), ylab = paste("Predicted", targ),
    xlim = c(8, 12), ylim = c(8,12),
    pch = 16, cex = 0.3, col = adjustcolor("black", alpha.f = 0.3),
    main = targ
  )
  abline(0, 1, col = "red", lwd = 2)
}
