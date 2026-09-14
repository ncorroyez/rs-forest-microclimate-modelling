# ==============================================================================
# What changes if the model is read NEAR THE GROUND instead of at 1 m?
# Asked by JO on 2026-08-20: the HOBOs sit against trunks, sheltered from air
# masses, so they might behave like a sensor closer to the soil.
#
# WHAT IS AND IS NOT POSSIBLE. MuSICA's air grid has no node below
# 0.0285 x Hmax: 0.085 m in the shortest stand here, 1.083 m in the tallest.
# A fixed 0.1 m is therefore unavailable for any plot above ~3.5 m, and
# micro_hourly_at() CLAMPS to level 1 rather than extrapolating, so asking for
# 0.1 m returns level 1 on 52 of the 53 loggers. The only well-defined
# experiment is the model's lowest node itself, nair == 1.
#
# CAVEAT THAT MUST TRAVEL WITH THE RESULT. nair == 1 is NOT a uniform downward
# shift. It reads 0.085 m in an open plot and 1.083 m in a tall stand, so the
# manipulation is confounded with Hmax, which is one of the figure's axes.
# Everything below is therefore also split by Hmax tercile.
#
# Reads : out_files/musica_hobo_native20/1111/*.nc, in_files/FR-Blo_2021_v2.nc,
#         CFG$hobo_temp_csv, tab_hobo_perplot_cluster.csv
# Writes: out_files/Chapter1/tables/tab_sampling_height_test.csv
#   Rscript scripts/c1_sampling_height_test.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4); library(lubridate); library(data.table)
  library(dplyr); library(musica.tools)})
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("pipeline/00_config.R")
DS <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
FORC <- "in_files/FR-Blo_2021_v2.nc"; SIMD <- "out_files/musica_hobo_native20/1111"
MREF <- macro_ref(FORC, DS); MOBS <- macro_ref_obs(FORC, DS)
macH <- {nc <- nc_open(FORC); tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
  t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tmac = as.numeric(ncvar_get(nc, "Tair")) - 273.15); nc_close(nc)
  d[as.Date(time) %in% DS][, .(Tmac = mean(Tmac, na.rm = TRUE)), by = time]}

hob <- as.data.table(read.csv(CFG$hobo_temp_csv))
hob[, datetime := as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
hob <- hob[position_sensor == "a" & as.Date(datetime) %in% DS & !id_plot %in% CFG$ids_to_remove]
hob[, time := floor_date(datetime, "hour")]
HO <- hob[, .(Tmic = mean(t_hobo, na.rm = TRUE)), by = .(id_plot, time)]
OBS <- HO[, {m <- .SD[, .(time, Tmic)]; mm <- merge(m, macH, by = "time")
  .(obs_dt = delta_tmax_mean(m, MOBS, DS),
    obs_sl = as.numeric(coef(lm(Tmic ~ Tmac, mm))[2]))}, by = id_plot]

lowest <- function(f) {                       # nair == 1, the model's lowest air node
  nc <- nc_open(f); on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  Tk <- ncvar_get(nc, "Tair_z")[1, ]
  data.table(time = floor_date(t0 + dhours(ncvar_get(nc, "time")), "hour"), Tmic = Tk - 273.15)
}
SIM <- rbindlist(lapply(OBS$id_plot, function(id) {
  f <- file.path(SIMD, sprintf("musica_out_HOBO_%s.nc", id))
  nc <- nc_open(f); z1 <- ncvar_get(nc, "relative_height")[1] *
    median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE); nc_close(nc)
  m1 <- micro_hourly_at(f, 1.0); m0 <- lowest(f)
  data.table(id_plot = id, z1 = z1,
             dt_1m  = delta_tmax_mean(m1, MREF, DS), dt_low = delta_tmax_mean(m0, MREF, DS),
             sl_1m  = as.numeric(coef(lm(Tmic ~ Tmac, merge(m1, macH, by = "time")))[2]),
             sl_low = as.numeric(coef(lm(Tmic ~ Tmac, merge(m0, macH, by = "time")))[2])) }))
D <- merge(OBS, SIM, by = "id_plot")
D <- merge(D, fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, Hmax)], by = "id_plot")
fwrite(D, "out_files/Chapter1/tables/tab_sampling_height_test.csv")

cat(sprintf("\nlowest air node: %.3f - %.3f m (median %.3f). Z = 1 m clamps to it on %d/53 plots\n",
            min(D$z1), max(D$z1), median(D$z1), sum(D$z1 >= 1)))
rep <- function(lab, s) cat(sprintf("  %-10s dTmax: r = %.3f | bias = %+.2f degC | amplitude = %3.0f%%   slope: r = %.3f | mean %.3f\n",
  lab, cor(D$obs_dt, D[[s[1]]]), mean(D[[s[1]]] - D$obs_dt), 100 * coef(lm(D[[s[1]]] ~ D$obs_dt))[2],
  cor(D$obs_sl, D[[s[2]]]), mean(D[[s[2]]])))
cat(sprintf("\n=== whole network (n = %d), observed bias reference: obs mean dTmax %+.2f degC ===\n", nrow(D), mean(D$obs_dt)))
rep("1 m",     c("dt_1m",  "sl_1m"))
rep("nair==1", c("dt_low", "sl_low"))

# The manipulation is confounded with height: split it so the confound is visible.
D[, terc := cut(Hmax, quantile(Hmax, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                labels = c("short", "mid", "tall"))]
cat("\n=== by Hmax tercile: how far down the reading actually moved, and what it did ===\n")
print(D[, .(n = .N, Hmax = round(mean(Hmax), 1), z_low = round(mean(z1), 2),
            bias_1m = round(mean(dt_1m - obs_dt), 2), bias_low = round(mean(dt_low - obs_dt), 2),
            dTmax_shift = round(mean(dt_low - dt_1m), 2),
            slope_1m = round(mean(sl_1m), 3), slope_low = round(mean(sl_low), 3)), by = terc][order(terc)])
cat(sprintf("\nlog(slope) range   1 m : %+.3f..%+.3f | nair==1 : %+.3f..%+.3f\n",
            log(min(D$sl_1m)), log(max(D$sl_1m)), log(min(D$sl_low)), log(max(D$sl_low))))
cat(sprintf("plots with slope > 1   1 m : %d | nair==1 : %d | observed : %d\n",
            sum(D$sl_1m > 1), sum(D$sl_low > 1), sum(D$obs_sl > 1)))
cat("DONE\n")
