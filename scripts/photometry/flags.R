library(data.table)
library(Rfits)
library(rhdf5)

filter_info = data.frame(fread("~/Documents/Epochs/data/throughputs/wavelengths.txt"))
phot_cat = Rfits_read("~/Documents/Epochs/data/catalogues/smiles_profound_multiband_cat.fits", data.table = TRUE)

png("~/Documents/Epochs/plots/flags/hotpixels.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "Flag hot pixels"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_hotpixel == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_hotpixel == 1], 
  pch = 16, col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_hotpixel)), 
  pch = 16, 
  col = "red"
)
dev.off()

png("~/Documents/Epochs/plots/flags/outside_nircam.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "outside nircam"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_outside_nircam == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_outside_nircam == 1], 
  pch = 16, col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_outside_nircam)), 
  pch = 16, 
  col = "red"
)
dev.off()

png("~/Documents/Epochs/plots/flags/unmatched_epochs.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "Unmatched epochs"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_epochs_unmatch == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_epochs_unmatch == 1], 
  pch = 16, col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_epochs_unmatch)), 
  pch = 16, 
  col = "red"
)
dev.off()

png("~/Documents/Epochs/plots/flags/potential_star.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "Star with similar F560W-F1280W colour as GAIA sources"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_star == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_star == 1], 
  pch = 16, col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_star)), 
  pch = 16, 
  col = "red"
)
dev.off()

png("~/Documents/Epochs/plots/flags/detection.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "Sources with <4sigma detection"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_det == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_det == 1], 
  pch = '.', col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_det)), 
  pch = 16, 
  col = "red"
)
dev.off()

png("~/Documents/Epochs/plots/flags/masked.png")
magplot(
  phot_cat$ProFound$RAmax,
  phot_cat$ProFound$Decmax, 
  xlab = "RA", 
  ylab = "DEC",
  pch = '.', 
  main = "Source compromised by masked region"
)
points(
  phot_cat$ProFound$RAmax[phot_cat$ProFound$flag_mask == 1],
  phot_cat$ProFound$Decmax[phot_cat$ProFound$flag_mask == 1], 
  pch = 16, col = "red"
)
legend(
  x = "topleft", 
  legend = paste0("Nflag = ", sum(phot_cat$ProFound$flag_mask)), 
  pch = 16, 
  col = "red"
)
dev.off()