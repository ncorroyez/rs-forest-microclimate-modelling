# ==============================================================================
# Chapter 3 — Full-period (Apr–Nov) LAI time-series per scenario, ONE-SIDED.
# Shows the exact Leaf_area the six scenarios feed MuSICA: the four S2-driven
# (dynamic) scenarios now carry REAL leaf-out / senescence from the newly inverted
# shoulder dates (Feb 24, Apr 23, Nov 9, Dec 21), not cyclic-spline extrapolation.
# The two LiDAR scenarios carry the parametric phenology (budburst 115) scaled to a
# single summer ALS magnitude. Real S2 acquisition dates marked; observation window
# shaded. Run from z_Example root:  Rscript fig_lai_timeseries_fullperiod.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(ggplot2); library(stringr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE)
df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt]
rA  <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD)
bod <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT)
ropt_p <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
fp     <- setNames(df$LAI_ALS/df$LAI_S2_ATBD, df$pid)   # S2 shape -> LiDAR full magnitude (ATBD-cal)
fA     <- setNames(df$LAI_ALS/df$LAI_S2_DOPT, df$pid)   # S2 shape -> opt-calibrated (~4.1)

# ---- rebuild ATBD daily S2 series from Not_Masked rasters (now incl. shoulders) ----
nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
acq_doy <- sort(as.integer(format(dates, "%j")))
stk <- rast(files)
pts <- vect(as.data.frame(df[, .(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[, -1, drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
cat(sprintf("series rebuilt: %d/53 plots, %d acquisition dates\n", length(atbd), length(dates)))

# ---- build the exact Leaf_area_1yr (one-sided) for each scenario, year 2021 ----
YR <- 2021L
pheno_s2  <- function(ts) make_phenology_from_s2(ts, YR)$Leaf_area_1yr[1:365]
pheno_als <- function(lai_max) calc_phenology(list.year=YR, nleafage=1, budburst_date=115,
              leaf_age_max_in=0.56, relative_age_firstmax=0.10, relative_age_lastmax=0.75,
              LAI_max_per_cohort=lai_max)$Leaf_area_1yr[1:365]

rows <- list()
for (i in seq_len(nrow(df))) {
  p <- df$pid[i]; id <- as.character(df$id_plot[i]); a <- atbd[[p]]
  als <- pheno_als(df$LAI_ALS[i]); alsd <- pheno_als(df$LAI_ALS_DOPT[i])
  rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="LiDAR full (static)",        LAI=als)
  rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="LiDAR d_opt (static)",       LAI=alsd)
  if (!is.null(a)) {
    ts_atbd  <- a
    ts_opt   <- data.frame(doy=a$doy, lai=a$lai*ropt_p[[p]])
    ts_rfull <- data.frame(doy=a$doy, lai=a$lai*rA)
    ts_ordpt <- data.frame(doy=a$doy, lai=a$lai*ropt_p[[p]]*bod)
    ts_fus   <- data.frame(doy=a$doy, lai=a$lai*fp[[p]])
    ts_fopt  <- data.frame(doy=a$doy, lai=a$lai*fA[[p]])
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2 ATBD (dyn)",            LAI=pheno_s2(ts_atbd))
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2 opt (dyn)",             LAI=pheno_s2(ts_opt))
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2 ATBD ×ratio→full (dyn)", LAI=pheno_s2(ts_rfull))
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2 opt ×ratio→d_opt (dyn)", LAI=pheno_s2(ts_ordpt))
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2t·ALS/ATBD (→LiDAR full)", LAI=pheno_s2(ts_fus))
    rows[[length(rows)+1]] <- data.table(id_plot=id, doy=1:365, scenario="S2t·ALS/opt (→~4.1)",        LAI=pheno_s2(ts_fopt))
  }
}
m <- rbindlist(rows, fill=TRUE)
lev <- c("S2 ATBD (dyn)","S2 opt (dyn)","S2 ATBD ×ratio→full (dyn)","S2 opt ×ratio→d_opt (dyn)",
         "S2t·ALS/ATBD (→LiDAR full)","S2t·ALS/opt (→~4.1)","LiDAR full (static)","LiDAR d_opt (static)")
m[, scenario := factor(scenario, levels=lev)]
m[, typ := fcase(grepl("S2t·", scenario), "Fusion",
                 grepl("LiDAR", scenario), "LiDAR", default="Sentinel-2")]
fwrite(m, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4s_lai_timeseries_fullperiod.csv")

# window: Apr 15 - Nov 15 (doy 105-319); per-plot median for the heavy line
WIN <- c(105, 319)
med <- m[doy %between% WIN, .(LAI=median(LAI)), by=.(scenario, typ, doy)]
acq_in <- acq_doy[acq_doy %between% WIN]
obs_win <- as.integer(format(as.Date(c("2021-05-28","2021-09-30")), "%j"))  # original summer window

g <- ggplot(m[doy %between% WIN], aes(doy, LAI, group=id_plot, colour=typ)) +
  annotate("rect", xmin=obs_win[1], xmax=obs_win[2], ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.5) +
  geom_line(alpha=0.10, linewidth=0.25) +
  geom_line(data=med, aes(doy, LAI, group=scenario), colour="black", linewidth=0.6, inherit.aes=FALSE) +
  geom_rug(data=data.frame(doy=acq_in), aes(x=doy), inherit.aes=FALSE, sides="b", colour="grey25", length=unit(0.03,"npc")) +
  facet_wrap(~scenario, ncol=3) +
  scale_colour_manual(values=PAL_SENSOR, name=NULL) +
  scale_x_continuous(breaks=as.integer(format(as.Date(paste0("2021-",sprintf("%02d",c(5,7,9,11)),"-01")),"%j")),
                     labels=c("May","Jul","Sep","Nov")) +
  labs(x="2021", y="LAI (one-sided, m² m⁻²)",
       subtitle="Per-plot LAI fed to MuSICA over Apr 15–Nov 15. Grey band = original summer S2 window (28 May–30 Sep); ticks = real S2 acquisitions; black = median.") +
  theme_article(10) + theme(legend.position="top")
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_lai_timeseries_fullperiod"), g, 9.5, 5.2)
ggsave_article("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/Fig_lai_timeseries_fullperiod", g, 9.5, 5.2)
cat("DONE -> Fig_lai_timeseries_fullperiod (+ Table4s)\n")
