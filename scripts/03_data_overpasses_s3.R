library(data.table); library(dplyr)

# ==============================================================================
# paths e directories
dir_data <- file.path("B:dati/icos_sites")
dir_out <- "outputs"
dir_icos <- file.path(dir_data, "Icos_archive")
file_meta <- file.path(dir_out, "icos_metadata.csv")
file_pass <- file.path(dir_data, "s3_overpasses.csv")
file_out <- file.path(dir_out, "s3_overpasses_tower.csv")


# ==============================================================================
# apertura dei file con i metadati e con i tempi
meta <- fread(file_meta)
pass <- fread(file_pass)
pass[, time_utc := as.POSIXct(time_utc, tz = "UTC")]

# eliminazione delle stazioni fuori dall'Europa
pass <- pass |> 
  filter(!site %in% c("CD-Ygb", "GF-Guy")) 
head(pass)
pass[, .(n = .N, first = as.Date(min(time_utc)), last = as.Date(max(time_utc))), by = site]

# nuova colonna con nome del satellite a partire da image_id
pass[, satellite := substr(image_id, 1, 3)] 
# unione di pass + utc_offset di ogni stazione
pass <- merge(pass, meta[, .(site, utc_offset)], by = "site")     # time zone of each station

# definizione di half-hour che contiene passaggio
pass[, hh_utc := as.POSIXct(floor(as.numeric(time_utc) / 1800) * 1800, origin = "1970-01-01", tz = "UTC")]
# conversione dell'ora nel formato del dataset icos
pass[, TIMESTAMP_START := format(hh_utc + utc_offset * 3600, "%Y%m%d%H%M", tz = "UTC")]
# minutes from the start of the half-hour (0-30)
pass[, min_in_hh := as.numeric(difftime(time_utc, hh_utc, units = "mins"))]
pass <- unique(pass, by = c("site", "TIMESTAMP_START"))            # one row per station and half-hour
nrow(pass); head(pass)


# ==============================================================================
# Variabili da estrarre da ICOS
# ==============================================================================
cols_hh <- c("LE_F_MDS", "LE_F_MDS_QC", "LE_CORR", "LE_CORR_25", "LE_CORR_75", "LE_RANDUNC",
             "H_F_MDS", "H_F_MDS_QC", "H_CORR", "NETRAD", "G_F_MDS", "G_F_MDS_QC",
             "NEE_VUT_REF", "NEE_VUT_REF_QC", "NEE_VUT_REF_JOINTUNC",
             "GPP_NT_VUT_REF", "GPP_NT_VUT_USTAR50", "GPP_NT_VUT_MEAN",
             "GPP_DT_VUT_REF", "GPP_DT_VUT_USTAR50", "GPP_DT_VUT_MEAN",
             "RECO_NT_VUT_REF", "RECO_DT_VUT_REF", "NIGHT", "SW_IN_POT",
             "TA_F", "VPD_F", "RH", "PA_F", "P_F", "SW_IN_F", "LW_IN_F", "SW_OUT", "LW_OUT",
             "WS", "WD", "USTAR")
cols_fp <- c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH", "ZL", "FETCH_70", "FETCH_80", "FETCH_90",
             "WSECT_FLAG", "LE_SSITC_TEST", "H_SSITC_TEST")

# SUPPORT FUNCTIONS
# read one file of a given type inside a zip
read_zip <- function(zip, site, type, ...) {
  files <- unzip(zip, list = TRUE)$Name
  file  <- grep(paste0(site, "_", type, "_"), files, value = TRUE, fixed = TRUE)
  tmp   <- unzip(zip, files = file, exdir = tempdir())
  on.exit(unlink(tmp))
  fread(tmp, na.strings = "-9999", ...)
}
# read requested columns that exist (+ soil sensors if soil = TRUE); -9999 -> NA
read_cols <- function(zip, site, type, cols, soil = FALSE) {
  header <- names(read_zip(zip, site, type, nrows = 0))
  extra  <- if (soil) grep("^(SWC|TS)_F_MDS_[0-9]+$", header, value = TRUE) else NULL
  d <- read_zip(zip, site, type, select = intersect(c("TIMESTAMP_START", cols, extra), header),
                colClasses = list(character = "TIMESTAMP_START"))
  for (v in names(d)) if (is.numeric(d[[v]])) set(d, which(d[[v]] == -9999), v, NA)
  d
}



# ==============================================================================
# Estrazione dei valori da ICOS
# ==============================================================================
zips <- list.files(dir_icos, pattern = "^ICOSETC_.*zip$", full.names = TRUE)
names(zips) <- substr(basename(zips), 9, 14) 

extract_site <- function(s) {
  ts <- pass[site == s, TIMESTAMP_START]
  hh <- read_cols(zips[s], s, "FLUXNET_HH", cols_hh, soil = TRUE)[TIMESTAMP_START %in% ts]
  fl <- read_cols(zips[s], s, "FLUXES", cols_fp)[TIMESTAMP_START %in% ts]
  setnames(hh, c("WS", "WD", "USTAR"), c("WS_HH", "WD_HH", "USTAR_HH"), skip_absent = TRUE)
  d  <- merge(hh, fl, by = "TIMESTAMP_START", all = TRUE)
  
  # se WS e WDnon sono in FLUXNET_HH allora vanno presi in FLUXES
  if (!"WS" %in% names(d) && "WS_HH" %in% names(d)) d$WS <- d$WS_HH   # wind fallback
  if (!"WD" %in% names(d) && "WD_HH" %in% names(d)) d$WD <- d$WD_HH
  d$site <- s
  d}

# creazione del dataset tower
tower <- rbindlist(lapply(unique(pass$site), extract_site), fill = TRUE)



# ==============================================================================
# Unione tower data con site info agli orari di passaggio
# ==============================================================================
out <- merge(pass, tower, by = c("site", "TIMESTAMP_START"), all.x = TRUE)       # keep every overpass
out <- merge(out, meta[, .(site, igbp, lat, lon, z_eddy, h_canopy)], by = "site", all.x = TRUE)
out[, zm := z_eddy - 0.67 * h_canopy]                                             # effective measurement height

# summary per station
fp_vars <- c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH")
out[, .(n           = .N,
        with_tower  = round(mean(!is.na(LE_F_MDS)), 2),
        fp_ok       = round(mean(complete.cases(.SD)), 2),
        le_measured = round(mean(LE_F_MDS_QC == 0, na.rm = TRUE), 2),
        zm          = round(zm[1], 1)),
    by = site, .SDcols = fp_vars]

fwrite(out, file_out)





