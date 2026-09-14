# ==============================================================================
# Does the fCover floor explain the model's flat response in open canopy?
# JO thread, 2026-08-21. Re-runs the 8 P1 loggers with their TRUE fractional
# cover instead of the 0.5 floor, everything else identical to the baseline
# (v3.2.3, abl_flag = 'iter', per-plot wind-corrected forcing, scan-angle LAI,
# native-20 m traits, 10-layer grid, read at 1 m).
#
# WHY NOT REUSE out_files/musica_hobo_nofloor/1111. That directory holds 4
# loggers dated 2026-07-10, a week BEFORE the native-20 baseline (2026-07-17),
# so it is not the same setup and a difference against the baseline would mix
# the floor with whatever else changed in between. Re-running is two minutes.
#
# The 4 non-floored P1 plots are re-run too: their cover is unchanged, so they
# must reproduce the baseline, which is the consistency check on this script.
# Writes: out_files/musica_hobo_native20_nofloor/1111/*.nc
#         out_files/Chapter1/tables/tab_nofloor_p1.csv
#   CH1_RUN_MUSICA=TRUE Rscript scripts/c1_nofloor_p1.R
# ==============================================================================
if (!identical(Sys.getenv("CH1_RUN_MUSICA"), "TRUE"))
  stop("gated stage: re-run with CH1_RUN_MUSICA=TRUE")
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(sf)
  library(terra); library(parallel); library(rmusica); library(musica.tools) })
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source)); source("pipeline/00_config.R")
FORC <- "in_files/FR-Blo_2021_v2.nc"
MB   <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
OUT  <- "out_files/musica_hobo_native20_nofloor/1111"
BASE <- "out_files/musica_hobo_native20/1111"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("R/wind_correction.R")

rasters <- load_lidar_rasters("in_files_native20", 1)
hob <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
                                       hobo_buffer_mode = "buffer25"))
# NO FLOOR HERE. This single missing line is the whole experiment.
sac <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")[, .(id_plot, sec_theta)]
hob <- merge(hob, sac, by = "id_plot", all.x = TRUE)
hob[, sec_theta := mean(hob$sec_theta, na.rm = TRUE)]
P1 <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[P == "P1", id_plot]
df <- as.data.frame(hob[id_plot %in% P1])
cat(sprintf("%d P1 plots, true fCover %.3f to %.3f (%d below the 0.5 floor)\n",
            nrow(df), min(df$fCover), max(df$fCover), sum(df$fCover < 0.5)))
mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
    leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
    LAI_max_per_cohort = lai1))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366; rbind(ph, d) }
invisible(mclapply(lapply(seq_len(nrow(df)), function(i) df[i, , drop = FALSE]), function(pr) {
  nc <- file.path(OUT, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
  if (file.exists(nc) && file.size(nc) > 1e6) return(NULL)
  lai <- pr$LAI / pr$sec_theta
  sc  <- list(lai_fn = function(p) lai, hmax_fn = function(p) pr$Hmax,
              fcover_fn = function(p) pr$fCover, lad_fn = make_lad_real,
              phenology_fn = mk_phen(lai))
  prc <- pr; prc$LAI <- lai
  tryCatch(run_musica_one(prc, sc, nc, windcorr_forcing(pr$Hmax), MB,
             extra_setup = list("abl_flag" = '"iter"')),
           error = function(e) cat(sprintf("ERR %s: %s\n", pr$id_plot, e$message)))
}, mc.cores = 4, mc.preschedule = FALSE))
cat(sprintf("nc written: %d / %d\n", length(list.files(OUT, "\\.nc$")), nrow(df)))

DS <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
MREF <- macro_ref(FORC, DS)
macH <- {nc <- nc_open(FORC); tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
  t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tmac = as.numeric(ncvar_get(nc, "Tair")) - 273.15); nc_close(nc)
  d[as.Date(time) %in% DS][, .(Tmac = mean(Tmac, na.rm = TRUE)), by = time]}
sc1 <- function(dir, id) { f <- file.path(dir, sprintf("musica_out_HOBO_%s.nc", id))
  if (!file.exists(f)) return(c(NA_real_, NA_real_))
  m <- micro_hourly_at(f, 1.0)
  c(delta_tmax_mean(m, MREF, DS),
    as.numeric(coef(lm(Tmic ~ Tmac, merge(m, macH, by = "time")))[2])) }
V <- fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")[id_plot %in% P1]
D <- rbindlist(lapply(V$id_plot, function(id) {
  a <- sc1(BASE, id); b <- sc1(OUT, id)
  data.table(id_plot = id, fl_dt = a[1], fl_sl = a[2], nf_dt = b[1], nf_sl = b[2]) }))
D <- merge(V, D, by = "id_plot")
D <- merge(D, as.data.table(df)[, .(id_plot, fCover_true = fCover, Hmax, LAI)], by = "id_plot")
D[, floored := fCover_true < 0.5]
fwrite(D, "out_files/Chapter1/tables/tab_nofloor_p1.csv")
cat("\n=== P1 loggers, floored baseline vs true cover, read at 1 m ===\n")
print(D[order(-obs_sl), .(id_plot, fCover = round(fCover_true, 3), planch = floored,
  Hmax = round(Hmax, 1), LAI = round(LAI, 2), obs = round(log(obs_sl), 3),
  sim_floor = round(log(fl_sl), 3), sim_true = round(log(nf_sl), 3),
  dTmax_floor = round(fl_dt, 2), dTmax_true = round(nf_dt, 2), dTmax_obs = round(obs_dt, 2))])
cat(sprintf("\nfloored plots (n=%d): log(slope) %+.3f -> %+.3f  (observed %+.3f)\n",
  sum(D$floored), mean(log(D[floored == TRUE, fl_sl])), mean(log(D[floored == TRUE, nf_sl])),
  mean(log(D[floored == TRUE, obs_sl]))))
cat(sprintf("control, unfloored plots (n=%d): %+.4f -> %+.4f (must be identical)\n",
  sum(!D$floored), mean(log(D[floored == FALSE, fl_sl])), mean(log(D[floored == FALSE, nf_sl]))))
cat("DONE\n")
