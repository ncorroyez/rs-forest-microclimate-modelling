# ==============================================================================
# Chapter 3 — continuous error-vs-structure: instead of a median split, plot the
# per-plot summer ΔTmax error (sim - obs) and absolute error against the canopy
# structure (LiDAR LAI, fractional cover, height) for the head-to-head scenarios,
# with a LOESS smooth, to locate the structure threshold where Sentinel-2 optical
# saturation overtakes LiDAR. 1 m, no shift, summer. Run from z_Example root.
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(STATIC_ALS="LiDAR full", DYN_S2_ATBD_NM="S2 ATBD (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full")   # head-to-head: structure vs raw vs rescaled S2
ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")

cov <- as.data.table(readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds"))[
  , .(id_plot=as.character(id_plot), LAI_ALS, fCover, Hmax)]
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, date=as.Date(time))]
OBS <- merge(hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date)], mday, by="date")[
  , .(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)), by=id_plot]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% ds_all, .(date=as.Date(time), Tsim=Tair_sim)]; r }
pp <- rbindlist(lapply(names(scen), function(scn){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL)
    sd <- merge(s[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=date], mday, by="date")
    data.table(scenario=scen[[scn]], id_plot=id, ds=mean(sd$Tmax_sim-sd$Tmax_macro,na.rm=TRUE)) }), fill=TRUE)
}), fill=TRUE)
pp <- merge(merge(pp, OBS, by="id_plot"), cov, by="id_plot")
pp[, `:=`(err=ds-do, aerr=abs(ds-do))]
pp[, scenario:=factor(scenario, levels=c("LiDAR full","S2 ATBD (dyn)","S2 ATBD ×ratio→full"))]
fwrite(pp, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4n_perplot_error_structure.csv")

# long over the 3 structure variables
pl <- melt(pp, id.vars=c("scenario","id_plot","err","aerr"),
           measure.vars=c("LAI_ALS","fCover","Hmax"), variable.name="svar", value.name="sval")
pl[, svar:=factor(svar, levels=c("LAI_ALS","fCover","Hmax"),
                  labels=c("(a) LiDAR LAI (one-sided)","(b) fractional cover","(c) canopy height (m)"))]

mk <- function(yv, ylab){
  ggplot(pl, aes(sval, get(yv), colour=scenario, fill=scenario)) +
    geom_hline(yintercept=if(yv=="err") 0 else NA, linetype=3, colour="grey50") +
    geom_point(size=1, alpha=0.5) +
    geom_smooth(method="loess", se=TRUE, span=1, linewidth=0.8, alpha=0.12) +
    facet_wrap(~svar, scales="free_x") +
    scale_colour_manual(values=PAL_SENSOR3, name=NULL) + scale_fill_manual(values=PAL_SENSOR3, name=NULL) +
    labs(x=NULL, y=ylab) + theme_article(10) + theme(legend.position="bottom")
}
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_error_vs_structure_signed"), mk("err","ΔTmax error (sim − obs, °C)"), 6.9, 3.6)
ggsave_article(file.path(outdir,"Fig_error_vs_structure_abs"),    mk("aerr","|ΔTmax error| (°C)"), 6.9, 3.6)

# crossing: LAI where S2 ATBD abs-error LOESS rises above LiDAR full
g_lai <- pp[scenario %in% c("LiDAR full","S2 ATBD (dyn)")]
gx <- seq(min(g_lai$LAI_ALS), max(g_lai$LAI_ALS), length.out=200)
fit <- function(sc){ d<-g_lai[scenario==sc]; predict(loess(aerr~LAI_ALS, d, span=1), gx) }
df <- data.table(LAI=gx, lidar=fit("LiDAR full"), s2=fit("S2 ATBD (dyn)"))
cross <- df[is.finite(lidar)&is.finite(s2)][which.min(abs(s2-lidar))]
cat(sprintf("\nLAI where S2 |error| ~ LiDAR |error| (abs-error LOESS crossing): LAI ≈ %.2f\n", cross$LAI))
cat(sprintf("   below it (LAI<%.1f): S2 abs-error <= LiDAR ; above it S2 worsens\n", cross$LAI))
cat("DONE -> Table4n + 2 figs\n")
