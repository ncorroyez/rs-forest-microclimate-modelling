# ==============================================================================
# c3_rebuild_annual_notmasked.R
#
# Rebuild the annual Sentinel-2 leaf-area series on the unmasked product, so
# that every one of the 53 Blois plots carries a series.
#
# The series shipped so far comes from `03_RESULTS/Blois/Metrics/Deciduous_Only`,
# whose rasters are masked at fCover >= 0.90 (the Chapter 2 filter that removes
# gaps and edges). Eight of the 53 plots have no valid pixel at any date under
# that mask, seven of them in the open archetype P1, so `smooth_s2_ts(min_obs =
# 4)` drops them and every dynamic scenario is undefined there.
#
# `Not_Masked` carries the same retrieval: on 2021-06-14 the two rasters agree
# to within 0 on all 266 034 shared pixels, and the unmasked one defines 25 961
# more. It also covers 19 dates over 2021 instead of 11. Rebuilding on it costs
# no change of product, only the removal of a mask.
#
# The previous `annual` series is not reproducible: its producer is absent from
# `build_all_s2_ts()`, the cached rds dates from 10 June 2026, and the series
# matches neither `rescaled` nor `atbd`. It is therefore redefined here rather
# than reverse-engineered: the shape comes from the smoothed unmasked ATBD
# series, and the magnitude is anchored so the summer peak equals LAI_ALS, the
# convention `build_dopt_pure_ts()` already uses.
#
# Input : 03_RESULTS/Blois/Metrics/Not_Masked/s2lai_<date>_atbd_res_10_m.tif
#         out_files/Chapter3/lai_prep/df_plots_real53.rds
# Output: out_files/Chapter3/lai_prep/ts_annual_notmasked.rds
#         out_files/Chapter3/tables/annual_notmasked_report.csv
#   Rscript c3_rebuild_annual_notmasked.R
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(data.table); library(mgcv)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config_CHS41.R")

NM       <- file.path("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics",
                      "Not_Masked")
DATE_ALL <- as.Date(c("2021-01-01", "2021-12-31"))
OUT_TS   <- "out_files/Chapter3/lai_prep/ts_annual_notmasked.rds"
OUT_REP  <- "out_files/Chapter3/tables/annual_notmasked_report.csv"

df  <- as.data.table(readRDS("out_files/Chapter3/lai_prep/df_plots_real53.rds"))
df[, pid := sprintf("X%d_Y%d", as.integer(x), as.integer(y))]
pts <- vect(as.matrix(df[, .(x, y)]), crs = "EPSG:32631")

stk <- assemble_s2_ts_stack(NM, "atbd", date_range = DATE_ALL)
dts <- as.Date(names(stk))
cat(sprintf("pile non masquee : %d dates, %s -> %s\n",
            length(dts), min(dts), max(dts)))

vals <- as.matrix(as.data.table(terra::extract(stk, pts))[, -1])
df[, n_obs := apply(vals, 1, function(z) sum(is.finite(z)))]
cat(sprintf("observations valides par placette : min %d, mediane %d\n",
            min(df$n_obs), median(df$n_obs)))

# Summer window used to set the anchor, matching the LiDAR acquisition season.
is_summer <- format(dts, "%m") %in% c("06", "07", "08", "09")

# ---- one plot -> smoothed daily series, peak anchored to LAI_ALS -------------
# `smooth_s2_ts()` is reused verbatim so the smoother stays the one the chapter
# already documents: cyclic spline on day-of-year, k capped by the number of
# observations. The anchor divides by the series' OWN summer peak, so a plot
# whose Sentinel-2 magnitude is wrong still ends up on the LiDAR scale.
build_one <- function(i) {
  lai <- vals[i, ]
  ok  <- is.finite(lai)
  if (sum(ok) < 4L) return(NULL)
  long <- data.frame(plot_id = df$pid[i],
                     doy = as.integer(format(dts[ok], "%j")),
                     lai = as.numeric(lai[ok]))
  s <- smooth_s2_ts(long, k = min(10L, floor(sum(ok) * 0.8)))
  if (!length(s)) return(NULL)
  d    <- s[[1]]
  peak <- max(d$lai[d$doy %in% as.integer(format(dts[is_summer & ok], "%j"))],
              na.rm = TRUE)
  if (!is.finite(peak) || peak <= 0) peak <- max(d$lai, na.rm = TRUE)
  if (!is.finite(peak) || peak <= 0) return(NULL)
  d$lai <- pmax(d$lai * (df$LAI_ALS[i] / peak), 0)
  d
}

ts_new <- setNames(lapply(seq_len(nrow(df)), build_one), df$pid)
ts_new <- Filter(Negate(is.null), ts_new)
cat(sprintf("\nseries construites depuis leurs propres pixels : %d / %d\n",
            length(ts_new), nrow(df)))
manque <- df[!pid %in% names(ts_new)]$id_plot
if (length(manque))
  cat("placettes sous le seuil de 4 observations :",
      paste(manque, collapse = ", "), "\n")

# ---- fallback: archetype-mean shape, same LAI_ALS anchor ---------------------
# A plot with fewer than four valid dates cannot carry its own phenology. It
# takes the mean shape of its archetype instead, anchored on its own LAI_ALS by
# the same rule as every other plot. The shape is borrowed, the magnitude is
# not. At 41_39 the anchor is 0.001 m2/m2, so the choice of shape moves the
# series by less than a thousandth of a unit of leaf area.
cl <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
df <- merge(df, cl[, .(id_plot, P)], by = "id_plot", all.x = TRUE)
borrowed <- character(0)
for (id in manque) {
  i    <- which(df$id_plot == id)
  peers <- df[P == df$P[i] & pid %in% names(ts_new)]$pid
  if (!length(peers)) next
  shape <- rowMeans(vapply(peers, function(k) {
    v <- ts_new[[k]]$lai; v / max(v)          # each peer on a 0-1 scale
  }, numeric(365)))
  ts_new[[df$pid[i]]] <- data.frame(doy = 1:365,
                                    lai = pmax(shape * df$LAI_ALS[i], 0))
  borrowed <- c(borrowed, id)
  cat(sprintf("  %s : forme moyenne de %s (%d placettes), ancree sur LAI_ALS = %.4f\n",
              id, df$P[i], length(peers), df$LAI_ALS[i]))
}
cat(sprintf("total : %d / %d series\n", length(ts_new), nrow(df)))

saveRDS(ts_new, OUT_TS)

# ---- report: what the new series looks like next to the old one -------------
old <- readRDS("out_files/Chapter3/lai_prep/ts_by_plot.rds")$annual
rep <- rbindlist(lapply(df$pid, function(k) {
  n <- ts_new[[k]]; o <- old[[k]]
  data.table(pid = k,
             n_obs      = df$n_obs[match(k, df$pid)],
             pic_new    = if (is.null(n)) NA_real_ else max(n$lai),
             ampl_new   = if (is.null(n)) NA_real_ else max(n$lai) - min(n$lai),
             pic_old    = if (is.null(o)) NA_real_ else max(o$lai),
             ampl_old   = if (is.null(o)) NA_real_ else max(o$lai) - min(o$lai))
}))
rep <- merge(rep, df[, .(pid, id_plot, LAI_ALS, fCover, P)], by = "pid")
setorder(rep, -n_obs)
fwrite(rep, OUT_REP)

cat("\n--- controle d'ancrage sur 3 placettes qui avaient deja une serie ---\n")
print(rep[!is.na(pic_old)][1:3, .(id_plot, LAI_ALS = round(LAI_ALS, 2),
                                  pic_new = round(pic_new, 2),
                                  pic_old = round(pic_old, 2),
                                  ampl_new = round(ampl_new, 2),
                                  ampl_old = round(ampl_old, 2))])
cat("\n--- les 8 placettes qui n'avaient rien ---\n")
bad <- c("41_17","41_18","41_19","41_27","41_30","41_39","41_47","41_49")
print(rep[id_plot %in% bad, .(id_plot, fCover = round(fCover, 3), n_obs,
                              LAI_ALS = round(LAI_ALS, 3),
                              pic_new = round(pic_new, 3))][order(n_obs)])
message("\nwrote ", OUT_TS, "\nwrote ", OUT_REP)
