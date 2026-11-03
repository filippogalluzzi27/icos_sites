# Per ogni stazione ICOS => un file con semiore e solo colonne di interesse

# ==============================================================================
library(data.table) 
source("setup/parameters/vars.r")
dir_zip <- "B:/dati/icos_sites/Icos_archive"
dir_out <- "B:/dati/icos_sites/intermediate_icos" 


# ==============================================================================
# primi controlli
# ==============================================================================
# elenco dei file
zips  <- list.files(dir_zip, pattern = "^ICOSETC_.*zip$", full.names = TRUE)
# codice del sito a partire dal nome del file
sites <- substr(basename(zips), 9, 14)
# tabella
data.table(sites, file = basename(zips))          


# ==============================================================================
# funzioni di supporto
# ==============================================================================
# leggere dentro ai file zip
read_zip <- function(zip, site, type, ...) {
  names <- unzip(zip, list = TRUE)$Name 
  f    <- grep(paste0(site, "_", type, "_"), names, value = TRUE, fixed = TRUE)  
  tmp  <- unzip(zip, files = f, exdir = tempdir())
  on.exit(unlink(tmp))
  fread(tmp, na.strings = "-9999", ...)
}

si <- read_zip(zips[1], sites[1], "SITEINFO"); si

# carico la tabella con le info dei siti
stations <- fread(file.path("outputs", "sites_info.csv"))
setnames(stations, c("fonte", "sito", "inizio", "fine"), c("source", "site", "start", "end"))
stations[, .N, by = site][N > 1]   # stations with more than one row: should be empty

# ==============================================================================
# stazione di prova per estrarre le variabili
# ==============================================================================
# selezione
site <- "BE-Bra"; zip  <- zips[sites == site]
header <- names(read_zip(zip, site, "FLUXNET_HH", nrows = 0))
# colonne da tenere
keep   <- intersect(c(cols_hh), header) 
# lettura variabili semiorarie da
hh <- read_zip(zip, site, "FLUXNET_HH", select = keep,
               colClasses = list(character = c("TIMESTAMP_START", "TIMESTAMP_END")))
dim(hh); setdiff(cols_hh, header)
# lettura variabili dei flussi non contenute in "FLUXNET_HH"
header_fp <- names(read_zip(zip, site, "FLUXES", nrows = 0))
fp <- read_zip(zip, site, "FLUXES", select = intersect(cols_fp, header_fp),
               colClasses = list(character = "TIMESTAMP_START"))
setdiff(cols_fp, header_fp)

# Impostare i nomi uguali
setnames(hh, c("WS", "WD", "USTAR"), c("WS_HH", "WD_HH", "USTAR_HH"), skip_absent = TRUE)
d <- merge(hh, fp, by = "TIMESTAMP_START", all.x = TRUE)
dim(hh); dim(fp); dim(d)

# impostare -9999 come NA
for (col in names(d)) {
  if (is.numeric(d[[col]])) set(d, which(d[[col]] == -9999), col, NA)
}

# sostituire le variabili se non sono presenti nel primo dataset
if (!"WS" %in% names(d) && "WS_HH" %in% names(d)) d$WS <- d$WS_HH
if (!"WD" %in% names(d) && "WD_HH" %in% names(d)) d$WD <- d$WD_HH

# UTC coordinate
offset     <- stations$utc[stations$site == site]
d$time_utc <- as.POSIXct(d$TIMESTAMP_START, format = "%Y%m%d%H%M", tz = "UTC") - offset * 3600
d$site     <- site
head(d[, .(TIMESTAMP_START, time_utc, site)])



# ==============================================================================
# ciclo per tutte le stazioni
# ==============================================================================
read_station <- function(site) {
  zip <- zips[sites == site]
  # FLUXNET_HH: requested columns + soil sensors that exist
  header <- names(read_zip(zip, site, "FLUXNET_HH", nrows = 0))
  soil   <- grep("^(SWC|TS)_F_MDS_[0-9]+$", header, value = TRUE)
  hh <- read_zip(zip, site, "FLUXNET_HH", select = intersect(c(cols_hh, soil), header),
                 colClasses = list(character = c("TIMESTAMP_START", "TIMESTAMP_END")))
  setnames(hh, c("WS", "WD", "USTAR"), c("WS_HH", "WD_HH", "USTAR_HH"), skip_absent = TRUE)
  # FLUXES: footprint inputs that exist
  header_fp <- names(read_zip(zip, site, "FLUXES", nrows = 0))
  fp <- read_zip(zip, site, "FLUXES", select = intersect(cols_fp, header_fp),
                 colClasses = list(character = "TIMESTAMP_START"))
  # merge by half-hour, -9999 -> NA, wind fallback
  d <- merge(hh, fp, by = "TIMESTAMP_START", all.x = TRUE)
  for (col in names(d)) if (is.numeric(d[[col]])) set(d, which(d[[col]] == -9999), col, NA)
  if (!"WS" %in% names(d) && "WS_HH" %in% names(d)) d$WS <- d$WS_HH
  if (!"WD" %in% names(d) && "WD_HH" %in% names(d)) d$WD <- d$WD_HH
  # UTC time and station code
  offset     <- stations$utc[stations$site == site]
  d$time_utc <- as.POSIXct(d$TIMESTAMP_START, format = "%Y%m%d%H%M", tz = "UTC") - offset * 3600
  d$site     <- site
  d   # the last line is what the function returns
}

dir.create(dir_out, showWarnings = FALSE)
for (s in sites) {
  saveRDS(read_station(s), file.path(dir_out, paste0(s, "_HH.rds")))
  cat(s, "saved\n")   # progress
}

# cheks delle variabili (per completezza)
inp   <- c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH")      # footprint inputs
check <- rbindlist(lapply(sites, function(s) {
  d <- readRDS(file.path(dir_out, paste0(s, "_HH.rds")))
  ok <- all(inp %in% names(d)) 
  data.table(site    = s,
             start   = as.Date(min(d$time_utc[!is.na(d$LE_F_MDS)])),
             end     = as.Date(max(d$time_utc[!is.na(d$LE_F_MDS)])),
             fp_ok   = if (ok) round(mean(complete.cases(d[, inp, with = FALSE])), 2) else 0,
             missing = paste(setdiff(inp, names(d)), collapse = " "))
}))
check$group <- ifelse(check$fp_ok > 0, "ready", "approx")
check$group[check$site %in% c("CD-Ygb", "GF-Guy")] <- "outside"
fwrite(check, file.path(dir_out, "stations_check.csv"))
check
fwrite(siti[fonte == "ICOS", .(site = sito, lat, lon)], "outputs/icos_sites_gee.csv")


