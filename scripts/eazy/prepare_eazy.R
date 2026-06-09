# ============================================================
# Script: find_filter_entries.R
# Purpose: Read FILTER_RES_latest.info and extract row indices
#          for HST ACS, JWST NIRCam, and JWST MIRI filters
#          matching the flux column names supplied.
# ============================================================

library(stringr)
library(data.table)
library(Rfits)

main_dir = "~/Documents/Epochs/"
fname = "smiles"
appsize = "0.32"
min_pc_err = 0.1

# load(
#   paste0(main_dir, "/data/prospect/prospect_data_for_fitting.rda")
# )

phot_cat = Rfits_read(
  paste0(main_dir, "/data/catalogues/", fname, "_profound_multiband_cat.fits")
)
filters_info = data.frame(fread(paste0(main_dir, "data/throughputs/wavelengths.txt")))

flux_cat = data.frame(phot_cat$Phot32)
df_flux = flux_cat[, paste0(filters_info$filters, "_flux_app_", appsize)]
df_flux_err = flux_cat[, paste0(filters_info$filters, "_scaled_flux_err_app_", appsize)]
df_flux_err = sqrt((df_flux_err)^2 + (min_pc_err * df_flux)^2) # put minimum error

# --- 1. Target flux column names -------------------------------------------
flux_cols = names(df_flux)
flux_err_cols = names(df_flux_err)

# Extract bare filter names, e.g. "F435W", "F850LP", "F560W" ...
filter_names <- str_extract(flux_cols, "^[^_]+")   # everything before first "_"
cat("Filters to locate:\n")
print(filter_names)
cat("\n")

# --- 2. Instrument look-up table -------------------------------------------
instrument_map <- list(
  HST_ACS    = c("F435W", "F606W", "F775W", "F814W", "F850LP"),
  JWST_NIRCam = c("F070W", "F090W", "F115W", "F150W", "F162M", "F182M",
                  "F200W", "F210M", "F250M", "F277W", "F300M", "F335M",
                  "F356W", "F410M", "F430M", "F444W", "F460M", "F480M"),
  JWST_MIRI  = c("F560W", "F770W", "F1000W", "F1280W", "F1500W",
                 "F1800W", "F2100W", "F2550W")
)

# --- 3. Read the .info file -------------------------------------------------
info_file <- "~/miniconda3/envs/eazy/lib/python3.14/site-packages/eazy/data/eazy-photoz/filters/FILTER.RES.latest.info"   # adjust path if needed
lines <- readLines(info_file)

# --- 4. Search patterns per instrument -------------------------------------
# The file uses these prefixes in the filter-curve path:
#   HST ACS     ->  hst/ACS_update_sep07/wfc_f<band>...
#   JWST NIRCam ->  jwst_nircam_f<band>
#   JWST MIRI   ->  jwst_miri_f<band>  (note: MIRI uses f0560w, f0770w, etc.)

# Build a named list:  filter_name -> regex that uniquely identifies its line
make_pattern <- function(filter_name, instrument) {
  fn_lower <- tolower(filter_name)
  if (instrument == "HST_ACS") {
    # Two versions exist in the file: _t77 (obs_before_7-4-06) and _t81 (obs_AFTER_7-4-06).
    # The script picks the FIRST hit (t77, pre-2006 throughput curve).
    # Change grep() to hits[2] below if you need the post-2006 (t81) version.
    return(paste0("(?i)wfc_", fn_lower, "[^a-z0-9]"))
  } else if (instrument == "JWST_NIRCam") {
    return(paste0("(?i)jwst_nircam_", fn_lower, "\\s"))
  } else if (instrument == "JWST_MIRI") {
    # MIRI names in the file use a leading zero for sub-100µm filters
    # e.g. F560W -> f0560w, F1000W -> f1000w
    miri_fn <- fn_lower
    # Insert leading zero for wavelengths < 1000 (3-digit names after 'f')
    miri_fn <- sub("^f(\\d{3}[mw])$", "f0\\1", miri_fn)
    return(paste0("(?i)jwst_miri_", miri_fn, "\\s"))
  }
}

# --- 5. Match each filter to a line ----------------------------------------
results <- data.frame(
  flux_col   = flux_cols,
  filter     = filter_names,
  instrument = NA_character_,
  line_num   = NA_integer_,
  filter_idx = NA_integer_,   # first integer token on the matched line
  line_text  = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(filter_names)) {
  fn   <- filter_names[i]
  inst <- names(Filter(function(v) fn %in% v, instrument_map))
  if (length(inst) == 0) {
    warning(sprintf("No instrument assigned for filter: %s", fn))
    next
  }
  results$instrument[i] <- inst
  
  pat   <- make_pattern(fn, inst)
  hits  <- grep(pat, lines, perl = TRUE)
  
  if (length(hits) == 0) {
    warning(sprintf("No match found in .info file for: %s (%s)", fn, inst))
  } else {
    if (length(hits) > 1) {
      warning(sprintf("Multiple matches for %s; using first hit (lines: %s)",
                      fn, paste(hits, collapse = ", ")))
    }
    ln  <- if (results$instrument[i] == "HST_ACS" && length(hits) >= 2) hits[2] else hits[1]
    results$line_num[i]   <- ln
    results$line_text[i]  <- trimws(lines[ln])
    # The filter catalogue index is the first integer on the line
    results$filter_idx[i] <- as.integer(str_extract(trimws(lines[ln]), "^\\d+"))
  }
}

# --- 6. Print summary -------------------------------------------------------
cat("=== Matched filter entries ===\n\n")
print(results[, c("flux_col", "instrument", "filter_idx", "line_num", "line_text")],
      right = FALSE, row.names = FALSE)

cat("\n=== Filter indices (to pass to EAzY/LePhare/etc.) ===\n")
cat(paste(results$filter_idx, collapse = " "), "\n")

# --- 7. Save to CSV ---------------------------------------------------------
eazy_translate = data.frame(cbind(
  c(
    flux_cols, 
    flux_err_cols
  ), 
  c(
    paste0("F", results$filter_idx), 
    paste0("E", results$filter_idx)
  )
))

eazy_nircam_translate = data.frame(cbind(
  c(
    flux_cols[!grepl("MIRI", results$instrument)], 
    flux_err_cols[!grepl("MIRI", results$instrument)]
  ), 
  c(
    paste0("F", results$filter_idx[!grepl("MIRI", results$instrument)]), 
    paste0("E", results$filter_idx[!grepl("MIRI", results$instrument)])
  )
))

fwrite(
  eazy_translate, 
  paste0("~/Documents/Epochs/data/eazy/eazy_", fname, "_zphot.translate"), 
  sep = " ", 
  row.names = FALSE
)

fwrite(
  eazy_nircam_translate, 
  paste0("~/Documents/Epochs/data/eazy/eazy_", fname, "_zphot_nircam.translate"), 
  sep = " ", 
  row.names = FALSE
)

fwrite(
  cbind(
    ID = phot_cat$ProFound$ID, RA = phot_cat$ProFound$RAmax, Dec = phot_cat$ProFound$Decmax, df_flux, df_flux_err
  ), 
  paste0("~/Documents/Epochs/data/eazy/", fname, "_phot_eazy.csv")
)
