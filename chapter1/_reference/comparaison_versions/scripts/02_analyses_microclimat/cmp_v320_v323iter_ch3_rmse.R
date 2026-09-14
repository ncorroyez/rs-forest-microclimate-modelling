# ==============================================================================
# Chapter 3 — v3.2.0 vs v3.2.3+iter validation RMSE grid (summer 2021), for the
# 13 scenarios with both binaries. Two metrics x two periods (c1 design):
#   metric : ΔTmax (daily, °C)        | slope micro/macro (per-plot, unitless)
#   period : all summer days (JJAS)   | hottest days (macro daily Tmax ≥ p90)
# RMSE is sim-vs-obs: ΔTmax over (plot,day); slope over plots (one slope/plot).
# IDENTICAL pipeline for both binaries; only the nc directory differs.
#   Rscript cmp_v320_v323iter_ch3_rmse.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr); library(data.table)
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

ROOTS <- c("v3.2.0"=file.path(CFG_C3$out_dir,"nc"), "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter"))
SCN <- c("STATIC_ALS","STATIC_ALS_DOPT","DYN_S2_ATBD_NM","DYN_S2_OPT_NM","DYN_ATBD_rfull_NM",
         "DYN_OPT_rdopt_NM","DYN_ALS_S2TIMING_NM","DYN_ALSoOPT_S2TIMING_NM",
         "FUSION_H","DYN_RF","DYN_S2_ANNUAL","CONST_ALS","NAIVE_S2_FORMSH")
LAB <- c(STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt", DYN_S2_ATBD_NM="S2 ATBD",
         DYN_S2_OPT_NM="S2 opt", DYN_ATBD_rfull_NM="S2 ATBD x->full", DYN_OPT_rdopt_NM="S2 opt x->d_opt",
         DYN_ALS_S2TIMING_NM="S2t.ALS/ATBD", DYN_ALSoOPT_S2TIMING_NM="S2t.ALS/opt",
         FUSION_H="Fusion-H (layered)", DYN_RF="S2 RF (dyn)", DYN_S2_ANNUAL="S2 annual (dyn)",
         CONST_ALS="LiDAR const", NAIVE_S2_FORMSH="S2 naive")
fam_of <- function(s) fifelse(grepl("STATIC|CONST_ALS",s),"LiDAR",
                       fifelse(grepl("S2TIMING|FUSION_H",s),"Fusion","S2"))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")    # JJAS

# ---- macro: hourly (slope) + daily Tmax + hottest-day set --------------------
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds][, `:=`(hr=floor_date(time,"hour"), date=as.Date(time))]
mday <- macro_h[, .(Tmax_macro=max(Tmacro)), by=date]
thr <- quantile(mday$Tmax_macro, 0.90, na.rm=TRUE); hot_days <- mday[Tmax_macro>=thr, date]
cat(sprintf("JJAS days=%d | hottest(p90>=%.1f°C)=%d\n", nrow(mday), thr, length(hot_days)))
dm <- as.data.frame(mday)   # for extract_deltatmax_*

# ---- observed: ΔTmax (daily) + slope (hourly), per plot ----------------------
obs_dt <- as.data.table(read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove))[
  , .(id_plot=as.character(id_plot), date, Delta_obs)]
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"), date=as.Date(time))]
hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
slp <- function(y,x) if (sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
obs_slope <- rbind(
  hb[, .(period="all", so=slp(Tobs,Tmacro)), by=id_plot],
  hb[date %in% hot_days, .(period="hot", so=slp(Tobs,Tmacro)), by=id_plot])

read_hourly <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(r)||!nrow(r))return(NULL)
  as.data.table(r)[, .(Tsim=Tair_sim, hr=floor_date(time,"hour"), date=as.Date(time))][date %in% ds] }

rmse <- function(a,b){ d<-a-b; sqrt(mean(d*d,na.rm=TRUE)) }
out <- list()
for (v in names(ROOTS)) for (scn in SCN) {
  d <- file.path(ROOTS[v], scn); fs <- list.files(d, pattern="\\.nc$", full.names=TRUE)
  if (!length(fs)) { cat(sprintf("  [%s] %s: MISSING\n", v, scn)); next }
  # ---- ΔTmax (daily) via canonical extractor ----
  dts <- as.data.table(extract_deltatmax_hobo_scenario(d, dm, ds))
  dts <- dts[, .(id_plot=as.character(id_plot), date, Delta_sim=Delta_Tmax)]
  DT <- merge(dts, obs_dt, by=c("id_plot","date"))
  # ---- slope (hourly) per plot ----
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)")
    s<-read_hourly(f); if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
  sim <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
  sim_slope <- rbind(
    sim[, .(period="all", ss=slp(Tsim,Tmacro)), by=id_plot],
    sim[date %in% hot_days, .(period="hot", ss=slp(Tsim,Tmacro)), by=id_plot])
  SL <- merge(sim_slope, obs_slope, by=c("id_plot","period"))
  for (per in c("all","hot")) {
    dd <- if (per=="all") DT else DT[date %in% hot_days]
    out[[length(out)+1]] <- data.table(version=v, scenario=scn, metric="dTmax", period=per,
      rmse=rmse(dd$Delta_sim, dd$Delta_obs), r=suppressWarnings(cor(dd$Delta_sim,dd$Delta_obs,use="complete.obs")),
      bias=mean(dd$Delta_sim-dd$Delta_obs,na.rm=TRUE), n=nrow(dd))
    ss <- SL[period==per]
    out[[length(out)+1]] <- data.table(version=v, scenario=scn, metric="slope", period=per,
      rmse=rmse(ss$ss, ss$so), r=suppressWarnings(cor(ss$ss,ss$so,use="complete.obs")),
      bias=mean(ss$ss-ss$so,na.rm=TRUE), n=nrow(ss))
  }
  cat(sprintf("  [%s] %s ok\n", v, scn))
}
R <- rbindlist(out)
R[, `:=`(label=factor(LAB[scenario], levels=unname(LAB)), fam=fam_of(scenario),
         version=factor(version, levels=names(ROOTS)),
         metric=factor(metric, levels=c("dTmax","slope")),
         period=factor(period, levels=c("all","hot")))]
fwrite(R, "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables/Table_v320_v323iter_ch3_rmse.csv")
cat("\n=== RMSE grid (sim-vs-obs) ===\n")
print(dcast(R, label+fam ~ version+metric+period, value.var="rmse")[order(fam,label)])

# ---- figure: RMSE, scenario x version, facet metric x period -----------------
plab <- c(all="All summer days", hot="Hottest days (p90)")
mlab <- c(dTmax="ΔTmax RMSE (°C)", slope="slope micro/macro RMSE")
R[, facet := factor(paste0(mlab[as.character(metric)], " — ", plab[as.character(period)]),
     levels=as.vector(outer(plab, mlab, function(p,m) paste0(m," — ",p))))]
g <- ggplot(R, aes(label, rmse, fill=version)) +
  geom_col(position=position_dodge(0.8), width=0.74) +
  facet_wrap(~facet, scales="free_y", ncol=2) +
  scale_fill_manual(values=c("v3.2.0"="#0072B2","v3.2.3 iter"="#D55E00"), name=NULL) +
  labs(x=NULL, y="RMSE sim vs observed (HOBO)",
       subtitle="Validation RMSE per scenario, summer 2021. Lower = closer to ground truth. v3.2.0 = article binary (validated); v3.2.3 = ABL-coupled control.") +
  theme_article(9) + theme(legend.position="top", axis.text.x=element_text(angle=40, hjust=1))
outdir <- "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/figures"
ggsave_article(file.path(outdir,"Fig_v320_v323iter_rmse_grid"), g, 10, 7)
cat("\nDONE -> Table_v320_v323iter_ch3_rmse.csv + Fig_v320_v323iter_rmse_grid\n")
