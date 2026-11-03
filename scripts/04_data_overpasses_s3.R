library(data.table)

# ==============================================================================
# paths e directories
dir_data <- file.path("B:/dati/icos_sites")
dir_out <- "outputs"
dir_icos <- file.path(dir_data, "Icos_archive")
file_meta <- file.path(dir_out, "icos_metadata.csv")
file_pass <- file.path(dir_data, "s3_overpasses.csv")
file_out <- file.path(dir_out, "s3_overpasses_tower.csv")



# ==============================================================================
# Orari a partire da Sentinel-3 e selezione
# ==============================================================================
meta <- fread(file_meta)
pass <- fread(file_pass)
pass[, time_utc := as.POSIXct(time_utc, tz = "UTC")]

# eliminazione delle stazioni fuori dall'Europa
pass <- pass[!site %in% c("CD-Ygb", "GF-Guy")]
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
# Variabili da estrarre da ICOS e support functions
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

# read one file of a given type inside a zip
read_zip <- function(zip, site, type, ...) {
  files <- unzip(zip, list = TRUE)$Name
  file  <- grep(paste0(site, "_", type, "_"), files, value = TRUE, fixed = TRUE)
  tmp   <- unzip(zip, files = file, exdir = tempdir())
  on.exit(unlink(tmp))
  fread(tmp, na.strings = "-9999", ...)
}

# read requested columns that exist; -9999 -> NA
read_cols <- function(zip, site, type, cols) {
  header <- names(read_zip(zip, site, type, nrows = 0))
  d <- read_zip(zip, site, type, select = intersect(c("TIMESTAMP_START", cols), header),
                colClasses = list(character = "TIMESTAMP_START"))
  for (v in names(d)) if (is.numeric(d[[v]])) set(d, which(d[[v]] == -9999), v, NA)
  d
}

# eddy system heights with their start date of validity
eddy_heights_dated <- function(zip, s) {
  vi <- read_zip(zip, s, "VARINFO_FLUXES")
  g  <- vi[VARIABLE == "VAR_INFO_VARNAME" & DATAVALUE %in% c("NEE", "FC", "USTAR"), GROUP_ID]
  if (length(g) == 0) return(NULL)
  h  <- dcast(vi[GROUP_ID %in% g & VARIABLE %in% c("VAR_INFO_HEIGHT", "VAR_INFO_DATE")],
              GROUP_ID ~ VARIABLE, value.var = "DATAVALUE")
  if (!"VAR_INFO_DATE" %in% names(h)) return(NULL)
  unique(data.table(site = s, from = as.Date(substr(h$VAR_INFO_DATE, 1, 8), "%Y%m%d"),
                    z = as.numeric(h$VAR_INFO_HEIGHT)))
}


# all canopy height groups of a station, one row per group
canopy_all <- function(zip, s) {
  an <- read_zip(zip, s, "ANCILLARY")[VARIABLE_GROUP == "GRP_HEIGHTC"]
  if (nrow(an) == 0) return(NULL)
  w <- dcast(an, GROUP_ID ~ VARIABLE, value.var = "DATAVALUE")
  col_or_na <- function(col) if (col %in% names(w)) w[[col]] else NA
  data.table(site = s, hc = as.numeric(col_or_na("HEIGHTC")), stat = col_or_na("HEIGHTC_STATISTIC"),
             species = col_or_na("HEIGHTC_SPP"), status = col_or_na("HEIGHTC_VEG_STATUS"),
             date = as.Date(substr(col_or_na("HEIGHTC_DATE_START"), 1, 8), "%Y%m%d"))
}

stat_order <- "Mean"     # same canopy statistic for every station

# canopy height values of one station: whole stand if available, otherwise tallest species
choose_canopy <- function(h) {
  h <- h[!is.na(hc) & (is.na(status) | status == "Alive") & stat %in% stat_order]
  if (nrow(h) == 0) return(NULL)
  level <- if (any(is.na(h$species))) "stand" else "species"
  if (level == "stand") h <- h[is.na(species)]
  best <- stat_order[stat_order %in% h$stat][1]
  h <- h[stat == best, .(hc = max(hc)), by = .(site, date)]
  h[, source := paste(level, best)]
}



# ==============================================================================
# Estrazione dei valori da ICOS
# ==============================================================================
zips <- list.files(dir_icos, pattern = "^ICOSETC_.*zip$", full.names = TRUE)
# assegnazione dei codici delle stazioni
names(zips) <- substr(basename(zips), 9, 14) 

extract_site <- function(s) {
  # prendere ts dal passaggio S3
  ts <- pass[site == s, TIMESTAMP_START]
  # estrazione deivalori per le variabili selezionate e per l'ora di passaggio
  hh <- read_cols(zips[s], s, "FLUXNET_HH", cols_hh)[TIMESTAMP_START %in% ts]
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
names(tower)
nrow(tower)


# ==============================================================================
# Unione tower data con site info agli orari di passaggio
# costruzione dell'altezza valida
# ==============================================================================
out <- merge(pass, tower, by = c("site", "TIMESTAMP_START"), all.x = TRUE)       # keep every overpass
out <- merge(out, meta[, .(site, igbp, lat, lon)], by = "site", all.x = TRUE)

# tabella delle altezze del sistema per stazione
z_tab  <- rbindlist(Map(eddy_heights_dated, zips, names(zips))); z_tab
# tabella delle altezze della chioma per stazione
hc_all <- rbindlist(Map(canopy_all, zips, names(zips)), fill = TRUE)
hc_tab <- rbindlist(lapply(split(hc_all, by = "site"), choose_canopy)); unique(hc_tab[, .(site, source)]); hc_tab
setorder(z_tab, site, from); setorder(hc_tab, site, date)

# aggiunta alla tabella out della colonna date con giorno di passaggio
out[, date := as.Date(time_utc)]

# Assegnazione a ogni passaggio dell’altezza del sistema in quel giorno
out[, z_eddy := {
  zt <- z_tab[site == .BY$site]
  if (nrow(zt) > 0) zt$z[pmax(1, findInterval(date, zt$from))] else NA_real_
}, by = site]

# Assegnazione a ogni passaggio dell’altezza della chioma in quel giorno
# canopy height: interpolated between dated surveys, otherwise median of undated values
out[, h_canopy := {
  ht <- hc_tab[site == .BY$site]
  hd <- ht[!is.na(date)]
  if (nrow(hd) > 1) approx(as.numeric(hd$date), hd$hc, xout = as.numeric(date), rule = 2, ties = mean)$y
  else if (nrow(ht) > 0) rep(if (nrow(hd) == 1) hd$hc else median(ht$hc), .N)
  else NA_real_
}, by = site]
out <- merge(out, unique(hc_tab[, .(site, hc_source = source)]), by = "site", all.x = TRUE)

# effective measurement height
out[, zm := z_eddy - 0.67 * h_canopy]                             

# summary per statione
fp_vars <- c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH")
out[, .(n           = .N,
        with_tower  = round(mean(!is.na(LE_F_MDS)), 2),
        fp_ok       = round(mean(complete.cases(.SD)), 2),
        le_measured = round(mean(LE_F_MDS_QC == 0, na.rm = TRUE), 2),
        zm          = round(zm[1], 1)),
    by = site, .SDcols = fp_vars]

out[, .(hc_ok = round(mean(!is.na(h_canopy)), 2)), by = site][hc_ok < 1]

fwrite(out, file_out)





