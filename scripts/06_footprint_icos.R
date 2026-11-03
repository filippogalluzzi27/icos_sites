library(data.table); library(sf); library(parallel); library(terra)


# ==============================================================================
# paths and directories
dir_out   <- "B:/dati/icos_sites/footprints_icos"
file_in   <- "outputs/s3_overpasses_selected.csv"   # output of the filter script
file_ffp  <- "scripts/FFP_R/calc_footprint_FFP.R"              # model code
blh_fixed <- 1000                                              # boundary layer height [m], placeholder until ERA5


# ==============================================================================
# dataset con valori semiorari di icos al passaggio di S3
d <- fread(file_in)[ok_all == TRUE]
nrow(d); d[, .N, by = site]
# chiamare la funzione del fooptinr
source(file_ffp)
args(calc_footprint_FFP)         # arguments of the function calc_footprint_FFP



# ==============================================================================
# Calcolo footprint per ogni stazione
# ==============================================================================
# calcolo UTM dalla latitudine
utm_epsg <- function(lon) 32600 + floor((lon + 180) / 6) + 1

# calcolo dei pixel del footprint di un passaggio
# pix = risoluzione dei pixels; frac = % per peso FP
footprint_pixels <- function(r, pix = 10, frac = 0.9) {
  # funzione del pacchetto
  fp <- tryCatch(calc_footprint_FFP(zm = r$zm, z0 = NaN, umean = r$WS, h = blh_fixed,
                                    ol = r$MO_LENGTH, sigmav = r$V_SIGMA, ustar = r$USTAR,
                                    wind_dir = r$WD),
                 error = function(e) NULL)
  if (is.null(fp) || fp$flag_err != 0) return(NULL)
  # posizione della torre in UTM
  epsg <- utm_epsg(r$lon)
  xy0  <- st_coordinates(st_transform(st_sfc(st_point(c(r$lon, r$lat)), crs = 4326), epsg))
  # per ogni punto del modello calcola su quale pixel cade
  pts <- data.table(px = floor((xy0[1] + as.vector(fp$x_2d)) / pix),
                    py = floor((xy0[2] + as.vector(fp$y_2d)) / pix),
                    w  = as.vector(fp$f_2d))
  # peso per ogni pixel
  pixels <- pts[w > 0, .(w = sum(w)), by = .(px, py)]
  pixels[, w := w / sum(w)]
  # tenere solo i pixels con valori entro il 90%
  setorder(pixels, -w)
  pixels <- pixels[cumsum(w) <= frac]
  pixels[, `:=`(id = r$id, epsg = epsg)]
  pixels
}

# ==============================================================================
# CHECK per una stazione
d[, id := paste0(site, "_", TIMESTAMP_START)]            # unique identifier of each overpass
x        <- d[site == "BE-Bra"][1]
pix_test <- footprint_pixels(x)
nrow(pix_test); sum(pix_test$w)                           # number of pixels, total weight (~0.9)

# pixel centres relative to the tower: the footprint must lie upwind
xy0 <- st_coordinates(st_transform(st_sfc(st_point(c(x$lon, x$lat)), crs = 4326), utm_epsg(x$lon)))
plot((pix_test$px + 0.5) * 10 - xy0[1], (pix_test$py + 0.5) * 10 - xy0[2], pch = 15, cex = 0.4, asp = 1,
     xlab = "x east of tower [m]", ylab = "y north of tower [m]")
points(0, 0, pch = 17, col = "red")                       # the tower


# ==============================================================================
# Run in parallelo per tutte le stazioni
n_cores <- detectCores(logical = FALSE) - 2              # physical cores, two left free
dir_fp <- normalizePath(file.path(dir_out, "fp_pixels"), mustWork = FALSE)
dir.create(dir_fp, showWarnings = FALSE)

# funzione di run => i dati sono salvati automaticamente
run_station <- function(s) {
  f <- file.path(dir_fp, paste0(s, "_fp.rds"))
  if (file.exists(f)) return(paste(s, "already done"))
  ds  <- d[site == s]
  res <- rbindlist(lapply(seq_len(nrow(ds)), function(i) footprint_pixels(ds[i])))
  saveRDS(res, f)
  paste(s, nrow(ds), "overpasses,", uniqueN(res$id), "footprints")
}

# impostazioni di run
ffp_path <- normalizePath(file_ffp)
t0 <- Sys.time()
cl <- makeCluster(n_cores)
clusterExport(cl, c("d", "dir_fp", "blh_fixed", "ffp_path", "utm_epsg", "footprint_pixels", "run_station"))
clusterEvalQ(cl, { library(data.table); library(sf); setDTthreads(1); source(ffp_path) })

# # RUN
# log <- parLapplyLB(cl, d[, .N, by = site][order(-N), site], run_station)
# stopCluster(cl)
# Sys.time() - t0
# unlist(log)

# ==============================================================================
# Risultati
check <- rbindlist(lapply(list.files(dir_fp, full.names = TRUE), function(f) {
  r <- readRDS(f)
  data.table(file = basename(f), footprints = uniqueN(r$id),
             pixels_per_fp = round(nrow(r) / uniqueN(r$id)))
}))
check
sum(check$footprints); nrow(d)                           # footprints obtained vs overpasses







# ==============================================================================
# random footprint from a random station
f   <- sample(list.files(dir_fp, full.names = TRUE), 1)
fps <- readRDS(f)
one <- fps[id == sample(unique(fps$id), 1)]

# pixels -> raster on the UTM grid (pixel centres = (index + 0.5) * 10 m)
r <- rast(one[, .(x = (px + 0.5) * 10, y = (py + 0.5) * 10, w)], type = "xyz",
          crs = paste0("EPSG:", one$epsg[1]))

# tower position in the same projection
tw <- d[id == one$id[1]]
tower <- project(vect(cbind(tw$lon, tw$lat), crs = "EPSG:4326"), crs(r))

plot(r, main = paste(one$id[1], " wind from", round(tw$WD), "deg"))
points(tower, pch = 17, col = "red")



