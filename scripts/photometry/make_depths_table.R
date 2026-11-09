filters_info = fread("~/Documents/JWST-MIRI-AGN/data/throughputs/wavelengths.txt")
smiles_frame_stats = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/depths/smiles_profound_frame_stats.csv"))
midis_frame_stats = data.frame(fread("~/Documents/JWST-MIRI-AGN/data/depths/midis_profound_frame_stats.csv"))

smiles_frame_stats[, paste0(filters_info$filters, "_5sigma_depth_32asec")]/1e-6
midis_frame_stats[, paste0(filters_info$filters, "_5sigma_depth_32asec")]/1e-6


safe_select <- function(df, cols) {
  df = data.frame(df)
  out <- data.frame(matrix(-999, nrow = nrow(df), ncol = length(cols)))
  names(out) <- cols
  existing <- cols[cols %in% names(df)]
  out[, existing] <- df[, existing]
  out
}

smiles_depths = safe_select(smiles_frame_stats, paste0(filters_info$filters, "_5sigma_depth_32asec"))
midis_depths = safe_select(midis_frame_stats, paste0(filters_info$filters, "_5sigma_depth_32asec"))

smiles_areas = safe_select(smiles_frame_stats, paste0(filters_info$filters, "_effarea_amin2"))
midis_areas = safe_select(midis_frame_stats, paste0(filters_info$filters, "_effarea_amin2"))

to_num = function(x) {
  if (is.list(x)) {
    # keeps length even if some elements are NULL or empty
    x = vapply(x, function(v) if (length(v) == 0) NA_real_ else as.numeric(v[1]), numeric(1))
  }
  as.numeric(x)
}

fmt = function(x, digits = 2) {
  x = to_num(x)
  x[x == -999] = NA
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "--")
}

mag = function(depth) -2.5 * log10(to_num(depth)) + 8.9

parts = strsplit(filters_info$filter_names, "_")
tel   = vapply(parts, `[`, character(1), 1)   # HST / JWST
inst  = vapply(parts, `[`, character(1), 2)   # ACSWFC / NIRCAM / MIRI
filt  = vapply(parts, `[`, character(1), 3)   # F435W, ...

inst_map = c(ACSWFC = "ACS", NIRCAM = "NIRCam", MIRI = "MIRI")

filter_labels = paste0("\\textit{", tel, "}/", unname(inst_map[inst]), " ", filt)

rows = paste0(
  filter_labels, " & ",
  fmt(mag(smiles_depths)), " & ",
  fmt(mag(midis_depths)), " & ",
  fmt(smiles_areas), " & ",
  fmt(midis_areas), " \\\\"
)

cat(rows, sep = "\n")

