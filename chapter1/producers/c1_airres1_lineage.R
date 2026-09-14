# ==============================================================================
# Parallel lineage at AIR_RESOLUTION_LEVEL = 1 (20 vegetation / 30 air layers),
# the 53 logger plots. Everything else identical to the chapter's baseline
# (c1_hobo_native20_validation.R): v3.2.3, abl_flag = 'iter', per-plot
# wind-corrected forcing, scan-angle-corrected LAI, native-20 m traits,
# fCover floored at 0.5.
#
# WHY. JO (2026-08-20): read the model near the ground, the loggers being against
# trunks. With the standard grid the lowest air node is 0.0285 x Hmax, so 1.083 m
# on the tallest plot, and nothing can be read below that in closed canopy, which
# is exactly where the +1.6 degC warm bias lives. Doubling the grid puts the
# lowest node at 0.00694 x Hmax = 0.264 m there (measured, not extrapolated, by
# scripts/c1_air_resolution_test.R).
#
# READ THE SCORING IN THIS ORDER. The finer grid also re-discretizes the canopy,
# so this is a DIFFERENT simulation, not a lower reading of the same one. The
# 1 m column isolates what the grid alone does; only then does the lowest-node
# column say anything about sampling height.
#
# This is a SIDE lineage. It writes its own directory and its own table and
# touches nothing the chapter reads.
# Writes: out_files/musica_hobo_native20_airres1/1111/*.nc
#         out_files/Chapter1/tables/tab_airres1_compare.csv
#   CH1_RUN_MUSICA=TRUE Rscript scripts/c1_airres1_lineage.R
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
OUT  <- "out_files/musica_hobo_native20_airres1/1111"
BASE <- "out_files/musica_hobo_native20/1111"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("R/wind_correction.R")

rasters <- load_lidar_rasters("in_files_native20", 1)
hob <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
                                       hobo_buffer_mode = "buffer25"))
hob[fCover < 0.5, fCover := 0.5]
sac <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")[, .(id_plot, sec_theta)]
hob <- merge(hob, sac, by = "id_plot", all.x = TRUE)
hob[, sec_theta := mean(hob$sec_theta, na.rm = TRUE)]
df <- as.data.frame(hob)
mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
    leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
    LAI_max_per_cohort = lai1))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366; rbind(ph, d) }

cat(sprintf("running %d plots at AIR_RESOLUTION_LEVEL = 1 ...\n", nrow(df)))
invisible(mclapply(lapply(seq_len(nrow(df)), function(i) df[i, , drop = FALSE]), function(pr) {
  nc <- file.path(OUT, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
  if (file.exists(nc) && file.size(nc) > 1e6) return(NULL)
  lai <- pr$LAI / pr$sec_theta
  sc  <- list(lai_fn = function(p) lai, hmax_fn = function(p) pr$Hmax,
              fcover_fn = function(p) pr$fCover, lad_fn = make_lad_real,
              phenology_fn = mk_phen(lai))
  prc <- pr; prc$LAI <- lai
  tryCatch(run_musica_one(prc, sc, nc, windcorr_forcing(pr$Hmax), MB,
             extra_setup = list("abl_flag" = '"iter"', "air_resolution_level" = 1L)),
           error = function(e) cat(sprintf("ERR %s: %s\n", pr$id_plot, e$message)))
}, mc.cores = 4, mc.preschedule = FALSE))
cat(sprintf("nc written: %d / %d\n", length(list.files(OUT, "\\.nc$")), nrow(df)))

# ---- scoring: grid effect at 1 m, then sampling-height effect ----------------
DS <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
MREF <- macro_ref(FORC, DS); MOBS <- macro_ref_obs(FORC, DS)
macH <- {nc <- nc_open(FORC); tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
  t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tmac = as.numeric(ncvar_get(nc, "Tair")) - 273.15); nc_close(nc)
  d[as.Date(time) %in% DS][, .(Tmac = mean(Tmac, na.rm = TRUE)), by = time]}
hb <- as.data.table(read.csv(CFG$hobo_temp_csv))
hb[, datetime := as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
hb <- hb[position_sensor == "a" & as.Date(datetime) %in% DS & !id_plot %in% CFG$ids_to_remove]
hb[, time := floor_date(datetime, "hour")]
HO <- hb[, .(Tmic = mean(t_hobo, na.rm = TRUE)), by = .(id_plot, time)]
OBS <- HO[, {m <- .SD[, .(time, Tmic)]; mm <- merge(m, macH, by = "time")
  .(obs_dt = delta_tmax_mean(m, MOBS, DS),
    obs_sl = as.numeric(coef(lm(Tmic ~ Tmac, mm))[2]))}, by = id_plot]
lowest <- function(f) { nc <- nc_open(f); on.exit(nc_close(nc))
  t0 <- as.POSIXct(sub("hours since ", "", ncatt_get(nc, "time", "units")$value), tz = "UTC")
  data.table(time = floor_date(t0 + dhours(ncvar_get(nc, "time")), "hour"),
             Tmic = ncvar_get(nc, "Tair_z")[1, ] - 273.15) }
sc_one <- function(dir, id, tag) {
  f <- file.path(dir, sprintf("musica_out_HOBO_%s.nc", id))
  if (!file.exists(f)) return(NULL)
  nc <- nc_open(f); z1 <- ncvar_get(nc, "relative_height")[1] *
    median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE); nc_close(nc)
  m1 <- micro_hourly_at(f, 1.0); m0 <- lowest(f)
  setNames(data.table(z1, delta_tmax_mean(m1, MREF, DS), delta_tmax_mean(m0, MREF, DS),
    as.numeric(coef(lm(Tmic ~ Tmac, merge(m1, macH, by = "time")))[2]),
    as.numeric(coef(lm(Tmic ~ Tmac, merge(m0, macH, by = "time")))[2])),
    paste0(tag, c("_z1", "_dt1m", "_dtlow", "_sl1m", "_sllow"))) }
D <- rbindlist(lapply(OBS$id_plot, function(id) {
  a <- sc_one(BASE, id, "g10"); b <- sc_one(OUT, id, "g20")
  if (is.null(a) || is.null(b)) return(NULL); cbind(data.table(id_plot = id), a, b) }))
D <- merge(OBS, D, by = "id_plot")
D <- merge(D, fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, Hmax)], by = "id_plot")
fwrite(D, "out_files/Chapter1/tables/tab_airres1_compare.csv")

rep <- function(lab, dt, sl) cat(sprintf("  %-26s dTmax r=%.3f bias=%+.2f amp=%3.0f%% | slope r=%.3f mean %.3f\n",
  lab, cor(D$obs_dt, D[[dt]]), mean(D[[dt]] - D$obs_dt), 100 * coef(lm(D[[dt]] ~ D$obs_dt))[2],
  cor(D$obs_sl, D[[sl]]), mean(D[[sl]])))
cat(sprintf("\n=== n = %d | lowest node: 10-layer %.2f-%.2f m, 20-layer %.2f-%.2f m ===\n",
  nrow(D), min(D$g10_z1), max(D$g10_z1), min(D$g20_z1), max(D$g20_z1)))
cat("\nSTEP 1, what the grid alone does (both read at 1 m):\n")
rep("10 layers @ 1 m (baseline)", "g10_dt1m", "g10_sl1m")
rep("20 layers @ 1 m",           "g20_dt1m", "g20_sl1m")
cat("\nSTEP 2, then what the sampling height does (20-layer grid):\n")
rep("20 layers @ 1 m",           "g20_dt1m", "g20_sl1m")
rep("20 layers @ lowest node",   "g20_dtlow", "g20_sllow")
D[, terc := cut(Hmax, quantile(Hmax, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                labels = c("short", "mid", "tall"))]
cat("\nby Hmax tercile (bias against the loggers):\n")
print(D[, .(n = .N, Hmax = round(mean(Hmax), 1), znode = round(mean(g20_z1), 2),
            g10_1m = round(mean(g10_dt1m - obs_dt), 2), g20_1m = round(mean(g20_dt1m - obs_dt), 2),
            g20_low = round(mean(g20_dtlow - obs_dt), 2)), by = terc][order(terc)])
cat("DONE\n")
