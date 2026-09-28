library(Rwcs)
library(Rfits)
library(ProFound)
library(data.table)
library(foreach)
library(stringr)
library(matrixStats)

field_name = 'midis'
# main_dir = '~/mnt/morgan/jdsilva/data/JP/'
main_dir = '/data/JP/'

data_dir = paste0(main_dir, '/ProFound/Data/', field_name, '/', field_name, '/')
profound_dir = paste0(main_dir, '/ProFound/Measurements/', field_name, '/', field_name, '/')
detect_dir = paste0(main_dir, '/ProFound/Detects/', field_name, '/', field_name, '/')

file_list = list.files(
  data_dir,
  pattern = glob2rx("*.fits"),
  full.names = TRUE
)
filt_names = toupper(unique( regmatches(file_list, regexpr("F\\d+[WML]P?", file_list, ignore.case = TRUE)) ))

calc_stats = function(pro_measure, filt_name){
  #effective observing area
  mask = pro_measure$mask
  # mask[pro_measure$] = 1L
  effarea = sum(pro_measure$mask == 0) * pixarea(pro_measure$keyvalues, unit = "amin2") ## unmasked area
  totarea = sum(!is.na(pro_measure$image)) * pixarea(pro_measure$keyvalues, unit = "amin2")
  
  ## depth
  random_aperture_radii = c(0.16, 0.32, 0.5, 0.6)
  random_apertures = profoundAperRan(
    image = pro_measure$image - pro_measure$sky,
    segim = pro_measure$segim,
    app_diam = random_aperture_radii * 2, 
    magzero = pro_measure$magzero,
    pixscale = pro_measure$pixscale,
    fluxtype = "Jansky", 
    Nran = 100000,
    depth = 5,
    correction = TRUE
  )
  df = data.frame(
    effarea,
    totarea,
    random_apertures$errors[1]*2,
    random_apertures$errors[2]*2,
    random_apertures$errors[3]*2,
    random_apertures$errors[1]*3,
    random_apertures$errors[2]*3,
    random_apertures$errors[3]*3,
    random_apertures$errors[1]*4,
    random_apertures$errors[2]*4,
    random_apertures$errors[3]*4,
    random_apertures$errors[1]*5,
    random_apertures$errors[2]*5,
    random_apertures$errors[3]*5, 
    pro_measure$magzero,
    pro_measure$pixscale
  )
  names(df) = c(
    paste0(filt_name, "_effarea_amin2"),
    paste0(filt_name, "_totarea_amin2"),
    paste0(filt_name, "_2sigma_depth_16asec"),
    paste0(filt_name, "_2sigma_depth_32asec"),
    paste0(filt_name, "_2sigma_depth_50asec"),
    paste0(filt_name, "_3sigma_depth_16asec"),
    paste0(filt_name, "_3sigma_depth_32asec"),
    paste0(filt_name, "_3sigma_depth_50asec"),
    paste0(filt_name, "_4sigma_depth_16asec"),
    paste0(filt_name, "_4sigma_depth_32asec"),
    paste0(filt_name, "_4sigma_depth_50asec"),
    paste0(filt_name, "_5sigma_depth_16asec"),
    paste0(filt_name, "_5sigma_depth_32asec"),
    paste0(filt_name, "_5sigma_depth_50asec"), 
    paste0(filt_name, "_magzero"),
    paste0(filt_name, "_pixscale")
  )
  return( df )
}

# Get statistics on the frames
frame_stats = foreach(i = 1:length(filt_names)) %do% {
  pro_measure_fstub = paste0(profound_dir, field_name, '_', field_name, '_', field_name, "_", filt_names[i], '_results.rds')
  message(pro_measure_fstub)
  
  pro_measure = readRDS(
    pro_measure_fstub
  )
  
  calc_stats(pro_measure, filt_names[i])
}
frame_stats = do.call(cbind, frame_stats)

pro_detect = readRDS(
  paste0(detect_dir, field_name, '_', field_name, '_', field_name, '_profound.rds')
)
det_stats = calc_stats(pro_detect, "detect")
# frame_stats$det_effarea_amin2 = sum(segim_stack$MASK[,]$imDat == 0L) * pixarea(segim_stack$MASK, unit = "amin2")
# frame_stats$det_2sigma_depth16asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.16)$y
# frame_stats$det_2sigma_depth32asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.32)$y
# frame_stats$det_2sigma_depth50asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.50)$y
# frame_stats$det_3sigma_depth16asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.16)$y
# frame_stats$det_3sigma_depth32asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.32)$y
# frame_stats$det_3sigma_depth50asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.50)$y
# frame_stats$det_4sigma_depth16asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.16)$y
# frame_stats$det_4sigma_depth32asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.32)$y
# frame_stats$det_4sigma_depth50asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.50)$y
# frame_stats$det_5sigma_depth16asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.16)$y
# frame_stats$det_5sigma_depth32asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.32)$y
# frame_stats$det_5sigma_depth50asec = approx(random_aperture_radii, y = random_apertures$errors, xout = 0.50)$y
# frame_stats$det_magzero = segim_stack$image$keyvalues$MAGZERO
# frame_stats$det_pixscale = pixscale(segim_stack$image$keyvalues)
fwrite(
  cbind(det_stats, frame_stats),
  paste0(main_dir, "/", field_name, "_profound_frame_stats.csv")
)
