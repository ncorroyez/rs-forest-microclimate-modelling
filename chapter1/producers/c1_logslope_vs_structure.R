# ==============================================================================
# log(micro-macro slope) against the three structural metrics, observed and
# simulated. Regeneration, on the current native-20 m lineage, of panel (c) of
# the SilviLaser poster (Plots/Silvilaser/Figure_Silvilaser.png, and its spring
# 2026 refresh Plots/Silvilaser/compfull.png), built from main_compare_era5_ign.R.
# Asked for by JO on 2026-08-20 to see how the panel moved across the summer.
#
# Slope = OLS coefficient of hourly T_micro (1 m) on hourly T_macro, JJAS 2021,
# one fit per logger, exactly as in the chapter's validation
# (c1_hobo_native20_validation_regen.R). log(slope) < 0 = buffering, > 0 =
# amplification. All 53 slopes are positive, so the log is defined everywhere.
#
# Two figures:
#   _native20   two series (field sensors, MuSICA now), the direct analogue of
#               the poster panel.
#   _evolution  three series, adding the pre-summer simulations. Each lineage is
#               regressed on ITS OWN forcing macro (musica_in_Blois.nc for the
#               pre-summer runs, FR-Blo_2021_v2.nc for the current ones), which
#               is what makes each series internally consistent; the x axis is
#               the CURRENT trait table in both cases, so the two model curves
#               sit on the same abscissa and only the simulation lineage moves.
#
# WHY THE POSTER'S AXES ARE NOT SUPERIMPOSABLE WITH THESE:
#   - traits recomputed natively at 20 m; VCI recomputed on ground-normalized
#     heights (the poster's VCI came from the un-normalized raster, hence its
#     0-0.73 range against 0.05-0.93 here).
#   - the third panel is the ONE-SIDED LAI used in the chapter, not the poster's
#     raw PAI (which ran past 12), so only the shape is comparable.
#
# Reads : out_files/Chapter1/tables/tab_hobo_native20_validation.csv  (stage A1)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv      (stage A2)
#         out_files/musica_hobo_ERA5_results/*.nc   (pre-summer runs, 2026-03-31)
# Writes: out_files/Chapter1/figures/Fig_logslope_vs_structure_native20.{png,pdf}
#         out_files/Chapter1/figures/Fig_logslope_vs_structure_evolution.{png,pdf}
#         out_files/Chapter1/tables/tab_logslope_vs_structure.csv
#   Rscript scripts/c1_logslope_vs_structure.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4); library(lubridate)
  library(data.table); library(ggplot2); library(grid)})
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("scripts/_article_style.R")
DS <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day"); Z <- 1
PRE_DIR <- "out_files/musica_hobo_ERA5_results"          # pre-summer lineage
CACHE   <- "out_files/Chapter1/tables/tab_logslope_presummer.csv"

macro_hourly <- function(f) {
  nc <- nc_open(f); tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
  t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tmac = as.numeric(ncvar_get(nc, "Tair")) - 273.15); nc_close(nc)
  d[as.Date(time) %in% DS][, .(Tmac = mean(Tmac, na.rm = TRUE)), by = time]
}
ols_slope <- function(m, M) {
  mm <- merge(m, M, by = "time")
  if (nrow(mm) > 50) as.numeric(coef(lm(Tmic ~ Tmac, mm))[2]) else NA_real_
}

V <- fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")
C <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, LAI, Hmax, VCI)]
D <- merge(V, C, by = "id_plot")[is.finite(obs_sl) & is.finite(sim_sl)]

# ---- pre-summer model slopes, cached (53 NetCDFs, ~2 min on a cold cache) -----
if (file.exists(CACHE)) {
  D <- merge(D, fread(CACHE), by = "id_plot")
} else {
  MOLD <- macro_hourly("in_files/musica_in_Blois.nc")
  P <- rbindlist(lapply(D$id_plot, function(id) {
    m <- micro_hourly_at(file.path(PRE_DIR, sprintf("musica_out_HOBO_%s.nc", id)), Z)
    data.table(id_plot = id,
               sim_sl_pre = if (is.null(m)) NA_real_ else ols_slope(m, MOLD)) }))
  fwrite(P, CACHE); D <- merge(D, P, by = "id_plot")
}
stopifnot(all(D$obs_sl > 0), all(D$sim_sl > 0), all(D$sim_sl_pre > 0, na.rm = TRUE))
D[, `:=`(obs = log(obs_sl), sim = log(sim_sl), pre = log(sim_sl_pre))]
fwrite(D, "out_files/Chapter1/tables/tab_logslope_vs_structure.csv")

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c(obs = "Field sensors (HOBO, n = 53)",
         pre = "MuSICA v3.2.0, before summer",
         sim = "MuSICA v3.2.3 iter, after summer")

long <- function(keys) {
  L <- rbindlist(lapply(names(MET), function(m) rbindlist(lapply(keys, function(k)
    data.table(metric = MET[[m]], x = D[[m]], y = D[[k]], src = SRC[[k]])))))
  L[is.finite(y)][, `:=`(metric = factor(metric, levels = unname(MET)),
                         src = factor(src, levels = unname(SRC[keys])))]
}
PAL <- setNames(c("#128d84", "#8c6bb1", "#E69F00"),        # poster colors kept
                SRC[c("obs", "pre", "sim")])
mk <- function(L) ggplot(L, aes(x, y, colour = src)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.9) +
  geom_point(alpha = 0.75, size = 1.9) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = "log(slope)") +
  theme_article(12) +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.background = element_blank(),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(4, "pt"))
ggsave_article("out_files/Chapter1/figures/Fig_logslope_vs_structure_native20",
               mk(long(c("obs", "sim"))), 10, 4.2)
ggsave_article("out_files/Chapter1/figures/Fig_logslope_vs_structure_evolution",
               mk(long(c("obs", "pre", "sim"))), 10, 4.2)

# ---- what moved across the summer -------------------------------------------
st <- long(c("obs", "pre", "sim"))[, {f <- lm(y ~ poly(x, 2))
  .(R2 = round(summary(f)$r.squared, 3),
    spearman = round(cor(x, y, method = "spearman"), 3))}, by = .(metric, src)]
cat("\n=== log(slope) vs structure, quadratic fit (n = ", nrow(D), ") ===\n", sep = "")
print(st)
sp <- function(v) sprintf("%+.3f..%+.3f (span %.3f)", min(v), max(v), max(v) - min(v))
cat(sprintf("\nobs                  %s\nMuSICA before summer %s\nMuSICA after summer  %s\n",
            sp(D$obs), sp(D$pre), sp(D$sim)))
cat(sprintf("share of the observed log(slope) span covered: before %.0f%%, after %.0f%%\n",
            100 * diff(range(D$pre)) / diff(range(D$obs)),
            100 * diff(range(D$sim)) / diff(range(D$obs))))
cat(sprintf("r(obs, model): before %.3f, after %.3f\n",
            cor(D$obs_sl, D$sim_sl_pre), cor(D$obs_sl, D$sim_sl)))
cat(sprintf("loggers with slope > 1: obs %d | before %d | after %d\n",
            sum(D$obs_sl > 1), sum(D$sim_sl_pre > 1), sum(D$sim_sl > 1)))
cat("DONE\n")
