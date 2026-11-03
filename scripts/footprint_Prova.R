path_zip_dir <- file.path("B:/dati/icos_sites/Icos_archive/ICOSETC_IT-Cp2_ARCHIVE_INTERIM_L2.zip")
path_zip_file <- file.path("ICOSETC_IT-Cp2_FLUXNET_HH_INTERIM_L2.csv")

df <- read.csv(unz(path_zip_dir, path_zip_file), na.strings = "-9999")


f_flux <- "ICOSETC_IT-Cp2_FLUXES_INTERIM_L2.csv"
fl <- read.csv(unz(path_zip_dir, f_flux), na.strings = "-9999",
               colClasses = c(TIMESTAMP_START = "character",
                              TIMESTAMP_END   = "character"))
dim(fl)

grep("USTAR|WD|WS|SIGMA|MO_LENGTH|ZL|FETCH", names(fl), value = TRUE)


vars <- c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH", "ZL",
          "FETCH_70", "FETCH_80", "FETCH_90", "FETCH_MAX", "WSECT_FLAG")

ff <- fl[, vars]
ff$t_start <- as.POSIXct(fl$TIMESTAMP_START, format = "%Y%m%d%H%M", tz = "Etc/GMT-1")
ff$t_end   <- as.POSIXct(fl$TIMESTAMP_END,   format = "%Y%m%d%H%M", tz = "Etc/GMT-1")
attr(ff$t_start, "tzone") <- "UTC"
attr(ff$t_end,   "tzone") <- "UTC"

summary(ff[, vars])

h <- as.numeric(format(ff$t_start, "%H"))
mattina <- ff[h >= 9 & h < 11, ]
round(colMeans(is.na(mattina[, vars])), 2)   # quota di dati mancanti per variabile

vi  <- read.csv(unz(path_zip_dir, "ICOSETC_IT-Cp2_VARINFO_FLUXES_INTERIM_L2.csv"))
anc <- read.csv(unz(path_zip_dir, "ICOSETC_IT-Cp2_ANCILLARY_INTERIM_L2.csv"))

vi[apply(vi, 1, function(r) any(grepl("HEIGHT", r))), ]
anc[apply(anc, 1, function(r) any(grepl("HEIGHTC", r))), ]



hc <- anc[anc$VARIABLE_GROUP == "GRP_HEIGHTC", c("GROUP_ID", "VARIABLE", "DATAVALUE")]
hc <- reshape(hc, idvar = "GROUP_ID", timevar = "VARIABLE", direction = "wide")
names(hc) <- sub("DATAVALUE.", "", names(hc), fixed = TRUE)
hc[hc$HEIGHTC_VEG_STATUS %in% "Alive" & hc$HEIGHTC_STATISTIC %in% "Mean",
   c("HEIGHTC", "HEIGHTC_SPP", "HEIGHTC_DATE_START")]

# PARAMETRI (provvisori)
z_mis <- 21       # altezza di misura [m], da VARINFO_FLUXES
h_c   <- 10       # altezza della chioma [m]: SOSTITUISCI con il valore scelto
blh   <- 1000     # altezza dello strato limite [m], provvisoria
d     <- 0.67 * h_c    # altezza di spostamento, regola empirica
zm    <- z_mis - d     # altezza efficace usata dal modello


ok <- format(ff$t_start, "%H:%M") == "10:00" &
  complete.cases(ff[, c("USTAR", "WD", "WS", "V_SIGMA", "MO_LENGTH")]) &
  ff$USTAR >= 0.1
cand <- ff[ok, ]
x <- cand[format(cand$t_start, "%m") == "07", ][1, ]   # primo caso utile di luglio
x

source("scripts/FFP_R/calc_footprint_FFP.R")   # nome del file scaricato, adattalo se diverso
args(calc_footprint_FFP)

fp <- calc_footprint_FFP(zm = zm, z0 = NA, umean = x$WS, h = blh,
                         ol = x$MO_LENGTH, sigmav = x$V_SIGMA, ustar = x$USTAR,
                         wind_dir = x$WD, r = c(50, 80))
str(fp, max.level = 1)

plot(fp$xr[[2]], fp$yr[[2]], type = "l", asp = 1, xlab = "x [m]", ylab = "y [m]")
lines(fp$xr[[1]], fp$yr[[1]], lty = 2)
points(0, 0, pch = 17)




library(terra)
si  <- read.csv(unz(path_zip_dir, "ICOSETC_IT-Cp2_SITEINFO_INTERIM_L2.csv"))
lat <- as.numeric(si$DATAVALUE[si$VARIABLE == "LOCATION_LAT"][1])
lon <- as.numeric(si$DATAVALUE[si$VARIABLE == "LOCATION_LONG"][1])

torre <- project(vect(cbind(lon, lat), crs = "EPSG:4326"), "EPSG:32633")
xy0   <- crds(torre)

p <- vect(cbind(as.vector(fp$x_2d) + xy0[1], as.vector(fp$y_2d) + xy0[2]),
          atts = data.frame(f = as.vector(fp$f_2d)), crs = "EPSG:32633")

e <- ext(p)
griglia <- rast(xmin = floor(e$xmin / 10) * 10, xmax = ceiling(e$xmax / 10) * 10,
                ymin = floor(e$ymin / 10) * 10, ymax = ceiling(e$ymax / 10) * 10,
                res = 10, crs = "EPSG:32633")
w <- rasterize(p, griglia, field = "f", fun = mean) * 100   # frazione di flusso per pixel

global(w, "sum", na.rm = TRUE)   # vicino a 1, un po' meno se il dominio taglia la coda
plot(w); points(torre, pch = 17)



dist80 <- function(h) {
  f <- calc_footprint_FFP(zm = zm, z0 = NA, umean = x$WS, h = h, ol = x$MO_LENGTH,
                          sigmav = x$V_SIGMA, ustar = x$USTAR, wind_dir = x$WD, r = 80)
  max(sqrt(f$xr[[1]]^2 + f$yr[[1]]^2), na.rm = TRUE)
}
sapply(c(500, 1000, 2000), dist80)

