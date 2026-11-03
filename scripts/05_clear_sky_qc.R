library(data.table)

# ==============================================================================
# paths and directories
dir_out  <- "outputs"
file_in  <- file.path(dir_out, "s3_overpasses_tower.csv")
file_out <- file.path(dir_out, "s3_overpasses_selected.csv")


# ==============================================================================
# selezione di dati con "clear-sky" in base alla radiazione
d <- fread(file_in); nrow(d)
d[, kt := SW_IN_F / SW_IN_POT]
quantile(d$kt, probs = seq(0, 1, 0.1), na.rm = TRUE)

# soglie di "cielo sereno"
kt_min  <- 0.65      # clear-sky index threshold
pot_min <- 100       # minimum potential radiation [W m-2] (excludes very low sun)
d[, ok_sky := !is.na(kt) & SW_IN_POT > pot_min & kt >= kt_min & kt <= 1]
d[, .N, by = ok_sky]                                        # overall
d[, .(n = .N, clear = round(mean(ok_sky), 2)), by = site]   # per station


# ==============================================================================
# Quality check per le torri
d[, ok_le   := !is.na(LE_F_MDS_QC) & LE_F_MDS_QC == 0 &
    (is.na(LE_SSITC_TEST) | LE_SSITC_TEST <= 1)]
d[, ok_rain := !is.na(P_F) & P_F == 0]
d[, .(le = sum(ok_le), test_missing = sum(ok_le & is.na(LE_SSITC_TEST)), no_rain = sum(ok_rain))]
d[, .(n = .N, le = round(mean(ok_le), 2), no_rain = round(mean(ok_rain), 2)), by = site]


# ==============================================================================
# Filtri su FP inputs
# esclusione di valori fisicamente impossibili
ustar_min <- 0.1; ustar_max <- 3      # friction velocity [m/s]
vsig_max  <- 5                        # lateral wind standard deviation [m/s]

d[, ok_fp := complete.cases(USTAR, WD, WS, V_SIGMA, MO_LENGTH, zm) &
    USTAR >= ustar_min & USTAR <= ustar_max &
    V_SIGMA > 0 & V_SIGMA <= vsig_max &
    WS > 0 & WD >= 0 & WD <= 360 & zm > 0]
d[, .N, by = ok_fp]
d[, .(n = .N, fp = round(mean(ok_fp), 2)), by = site]


# ==============================================================================
# merge
d[, ok_all := ok_sky & ok_le & ok_rain & ok_fp]

# overpasses left after each filter, applied one after the other
d[, .(total        = .N,
      sky          = sum(ok_sky),
      sky_le       = sum(ok_sky & ok_le),
      sky_le_rain  = sum(ok_sky & ok_le & ok_rain),
      all          = sum(ok_all))]

# per station, sorted by number of selected overpasses
d[, .(n = .N, selected = sum(ok_all), share = round(mean(ok_all), 2)), by = site][order(-selected)]

fwrite(d, file_out)


