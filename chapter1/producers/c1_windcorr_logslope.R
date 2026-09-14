# ==============================================================================
# Does the canopy wind correction explain why the model no longer amplifies in
# open canopy? JO, 2026-08-21, after the ABL iteration turned out to explain only
# the dense end.
#
# NO NEW SIMULATION. The correction levels already exist on disk (the ledger of
# c1_correction_levels_rescore.R). The clean wind contrast is uncorr vs windonly:
# both carry the raw LiDAR leaf area, so only the in-canopy wind profile moves.
# The baseline adds the scan-angle leaf-area correction on top, so it is shown
# but is not the wind contrast.
# Reads : out_files/hobo_native20_{uncorr,windonly}/*.nc
#         out_files/musica_hobo_native20/1111/*.nc
# Writes: out_files/Chapter1/tables/tab_windcorr_logslope.csv
#         out_files/Chapter1/figures/Fig_logslope_windcorr.{png,pdf}
#   Rscript scripts/c1_windcorr_logslope.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4); library(lubridate); library(data.table)
  library(ggplot2); library(grid)})
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source)); source("pipeline/00_config.R")
source("scripts/_article_style.R")
DS <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
FORC <- "in_files/FR-Blo_2021_v2.nc"
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

VAR <- data.table(
  vkey  = c("uncorr", "windonly", "full"),
  dir   = c("out_files/hobo_native20_uncorr", "out_files/hobo_native20_windonly",
            "out_files/musica_hobo_native20/1111"),
  label = c("MuSICA, no wind correction", "MuSICA, wind correction",
            "MuSICA, wind + scan angle (baseline)"))
sc <- function(dir, id) {
  f <- file.path(dir, sprintf("musica_out_HOBO_%s.nc", id))
  if (!file.exists(f)) return(c(NA_real_, NA_real_))
  m <- micro_hourly_at(f, 1.0); if (is.null(m)) return(c(NA_real_, NA_real_))
  c(delta_tmax_mean(m, MREF, DS),
    as.numeric(coef(lm(Tmic ~ Tmac, merge(m, macH, by = "time")))[2])) }
D <- copy(OBS)
for (i in seq_len(nrow(VAR))) {
  v <- t(sapply(D$id_plot, function(id) sc(VAR$dir[i], id)))
  D[[paste0(VAR$vkey[i], "_dt")]] <- v[, 1]; D[[paste0(VAR$vkey[i], "_sl")]] <- v[, 2] }
D <- merge(D, fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[
  , .(id_plot, LAI, VCI, Hmax, P)], by = "id_plot")
fwrite(D, "out_files/Chapter1/tables/tab_windcorr_logslope.csv")

rep <- function(k, lab) { dt <- D[[paste0(k, "_dt")]]; sl <- D[[paste0(k, "_sl")]]
  cat(sprintf("%-38s dTmax r=%.3f bias=%+.2f amp=%3.0f%% | log(slope) OLS %3.0f%% span %3.0f%% | r=%.3f | range %+.3f..%+.3f\n",
    lab, cor(D$obs_dt, dt), mean(dt - D$obs_dt), 100 * coef(lm(dt ~ D$obs_dt))[2],
    100 * coef(lm(log(sl) ~ log(D$obs_sl)))[2],
    100 * diff(range(log(sl))) / diff(range(log(D$obs_sl))),
    cor(D$obs_sl, sl), log(min(sl)), log(max(sl)))) }
cat(sprintf("\nobserved                               log(slope) range %+.3f..%+.3f\n",
            log(min(D$obs_sl)), log(max(D$obs_sl))))
for (i in seq_len(nrow(VAR))) rep(VAR$vkey[i], VAR$label[i])
cat("\nmean log(slope) per archetype:\n")
print(D[, .(n = .N, obs = round(mean(log(obs_sl)), 3),
            uncorr = round(mean(log(uncorr_sl)), 3), wind = round(mean(log(windonly_sl)), 3),
            full = round(mean(log(full_sl)), 3),
            obs_max = round(max(log(obs_sl)), 3), uncorr_max = round(max(log(uncorr_sl)), 3),
            wind_max = round(max(log(windonly_sl)), 3)), by = P][order(P)])

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c("Field sensors (HOBO, 1 m)", VAR$label[1:2])
L <- rbindlist(lapply(names(MET), function(m) rbindlist(list(
  data.table(metric = MET[[m]], x = D[[m]], y = log(D$obs_sl),      src = SRC[1]),
  data.table(metric = MET[[m]], x = D[[m]], y = log(D$uncorr_sl),   src = SRC[2]),
  data.table(metric = MET[[m]], x = D[[m]], y = log(D$windonly_sl), src = SRC[3])))))
L <- L[is.finite(x) & is.finite(y)]
L[, `:=`(metric = factor(metric, levels = unname(MET)), src = factor(src, levels = SRC))]
PAL <- setNames(c("#128d84", "#D55E00", "#E69F00"), SRC)
p <- ggplot(L, aes(x, y, colour = src)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.9) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = "log(slope)") + theme_article(12) +
  theme(legend.position = "bottom", strip.placement = "outside", strip.background = element_blank(),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(4, "pt")) +
  guides(colour = guide_legend(nrow = 1))
ggsave_article("out_files/Chapter1/figures/Fig_logslope_windcorr", p, 10.5, 4.2)
cat("DONE\n")
