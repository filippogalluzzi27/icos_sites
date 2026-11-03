# script per la creazione del dataset contenente le informazioni principali
# dei siti, le coordinate, e delle variabili di sito
# Sia per ICOS che per Fluxnet



# ==============================================================================
# ICOS
# ==============================================================================
library(data.table)

# Paths
PARAMETERS <- list(
  dir_icos = "B:/dati/icos_sites/Icos_archive",
  file_out = "outputs/sites_info_icos.csv")

# funzione che legge i file dentro a cartella zip
read_zip <- function(zip, site, type, ...) {
  files <- unzip(zip, list = TRUE)$Name
  file  <- grep(paste0(site, "_", type, "_"), files, value = TRUE, fixed = TRUE)
  tmp   <- unzip(zip, files = file, exdir = tempdir())
  on.exit(unlink(tmp))
  fread(tmp, na.strings = "-9999", ...)
}

# value of a variable in a long-format table 
get_value <- function(d, var) d$DATAVALUE[d$VARIABLE == var][1]

# ICOS archives and station codes
zips  <- list.files(PARAMETERS$dir_icos, pattern = "^ICOSETC_.*zip$", full.names = TRUE)
sites <- substr(basename(zips), 9, 14)

# ==============================================================================
# test on one station: columns from SITEINFO
si <- read_zip(zips[1], sites[1], "SITEINFO")
data.table(site = sites[1], 
           site_name = get_value(si, "SITE_NAME"), 
           igbp = get_value(si, "IGBP"),
           lat = as.numeric(get_value(si, "LOCATION_LAT")), 
           lon = as.numeric(get_value(si, "LOCATION_LONG")),
           elevation = as.numeric(get_value(si, "LOCATION_ELEV")), 
           utc_offset = as.numeric(get_value(si, "UTC_OFFSET")),
           homogeneity = get_value(si, "SURFACE_HOMOGENEITY"))



# ==============================================================================
# one row of site information per ICOS station
info_icos <- function(zip, site) {
  si <- read_zip(zip, site, "SITEINFO")
  yy <- read_zip(zip, site, "FLUXNET_YY", select = "TIMESTAMP")
  data.table(site = site, 
             site_name = get_value(si, "SITE_NAME"), 
             country = substr(site, 1, 2),
             source = "ICOS", 
             archive = basename(zip), 
             igbp = get_value(si, "IGBP"),
             lat = as.numeric(get_value(si, "LOCATION_LAT")), 
             lon = as.numeric(get_value(si, "LOCATION_LONG")),
             elevation = as.numeric(get_value(si, "LOCATION_ELEV")),
             utc_offset = as.numeric(get_value(si, "UTC_OFFSET")),
             mat = as.numeric(get_value(si, "MAT")), 
             map = as.numeric(get_value(si, "MAP")),
             homogeneity_m = as.numeric(get_value(si, "SURFACE_HOMOGENEITY")),
             start = min(yy$TIMESTAMP), 
             end = max(yy$TIMESTAMP),
)}

info_icos(zips[1], sites[1])                    # test on one station

icos <- rbindlist(Map(info_icos, zips, sites)); icos
sum(is.na(icos$h_canopy))
write.csv(icos, "outputs/icos_metadata.csv", row.names = FALSE)



# ==============================================================================
# FLUXNET
# ==============================================================================
PARAMETERS$dir_fx <- "B:/dati/icos_sites/Fluxnet_archive"

# Fluxnet archives and station codes
zips_fx  <- list.files(PARAMETERS$dir_fx, pattern = "_FLUXNET_[0-9]{4}-[0-9]{4}.*zip$", full.names = TRUE)
sites_fx <- regmatches(basename(zips_fx), regexpr("[A-Z]{2}-[A-Za-z0-9]{3}", basename(zips_fx)))
length(zips_fx); any(duplicated(sites_fx))

# test: files inside one zip and variables in its BIF
unzip(zips_fx[1], list = TRUE)$Name
bif <- read_zip(zips_fx[1], sites_fx[1], "FLUXNET_BIF")
unique(bif$VARIABLE)

# one row of site information per Fluxnet station
info_fx <- function(zip, site) {
  b  <- read_zip(zip, site, "FLUXNET_BIF")
  yr <- as.integer(strsplit(regmatches(basename(zip), regexpr("[0-9]{4}-[0-9]{4}", basename(zip))), "-")[[1]])
  data.table(site = site, site_name = get_value(b, "SITE_NAME"), country = substr(site, 1, 2),
             source = "FLUXNET", network = sub("_.*", "", basename(zip)), archive = basename(zip),
             igbp = get_value(b, "IGBP"),
             lat = as.numeric(get_value(b, "LOCATION_LAT")), lon = as.numeric(get_value(b, "LOCATION_LONG")),
             elevation = as.numeric(get_value(b, "LOCATION_ELEV")),
             utc_offset = as.numeric(get_value(b, "UTC_OFFSET")),
             mat = as.numeric(get_value(b, "MAT")), map = as.numeric(get_value(b, "MAP")),
             homogeneity_m = as.numeric(get_value(b, "SURFACE_HOMOGENEITY")),
             start = yr[1], end = yr[2],
)}

info_fx(zips_fx[1], sites_fx[1])                       # test on one station

fx <- rbindlist(Map(info_fx, zips_fx, sites_fx)); fx
colSums(is.na(fx))                                     # missing values per column
write.csv(fx, "outputs/fluxnet_metadata.csv", row.names = FALSE)




# ==============================================================================
# COMPARISON
# ==============================================================================
icos <- read.csv("outputs/icos_metadata.csv")
fluxn <- read.csv("outputs/fluxnet_metadata.csv")

# stations present in both archives, with the two periods side by side
both <- merge(icos[, c("site", "start", "end")],
              fluxn[, c("site", "start", "end")],
              by = "site", suffixes = c("_icos", "_fx"))
both <- both |> 
  relocate(site, start_fx, start_icos,end_fx, end_fx); both

setdiff(icos$site, fluxn$site)    # stations only in ICOS
setdiff(fluxn$site, icos$site)    # stations only in Fluxnet












