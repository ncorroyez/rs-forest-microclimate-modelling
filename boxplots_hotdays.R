# ==============================================================================
# Chapter 3 — boxplots restricted to the 10% HOTTEST DAYS of the period (macro
# daily Tmax ≥ p90 over Apr 15 – Nov 15). Per-plot ΔTmax and slope (Tair 1 m, no
# shift) for the six scenarios + observed HOBO. Two layouts: faceted by month
# (only months containing hot days), scenarios on x. Run from z_Example root.
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt",
          DYN_S2_ATBD_NM="S2 ATBD", DYN_S2_OPT_NM="S2 opt",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt")
WIN <- as.Date(c("2021-04-15","2021-11-15")); ds_all <- seq(WIN[1], WIN[2], by="day")
mlab <- function(d) format(d, "%b")

ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, hr:=floor_date(time,"hour")]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]
thr <- quantile(mday$Tmax_macro, 0.90, na.rm=TRUE); hot <- mday[Tmax_macro>=thr]$date
cat(sprintf("p90 threshold = %.1f °C ; %d hot days\n", thr, length(hot)))
print(table(format(hot,"%b")))
macro_h <- macro_h[as.Date(time) %in% hot]; macro_h[, mo:=mlab(as.Date(time))]
mday <- mday[date %in% hot][, mo:=mlab(date)]
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_

hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% hot,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"), mo=mlab(as.Date(time)), date=as.Date(time))]
hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
o_dT <- merge(hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date,mo)], mday[,.(date,Tmax_macro)], by="date")[
  , .(value=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)), by=.(id_plot,mo)][, `:=`(source="Observed (HOBO)", metric="dTmax")]
o_sl <- hb[, .(value=slp(Tobs,Tmacro)), by=.(id_plot,mo)][, `:=`(source="Observed (HOBO)", metric="slope")]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% hot, .(Tsim=Tair_sim, hr=floor_date(time,"hour"), date=as.Date(time), mo=mlab(as.Date(time)))]; r }

per <- list(o_dT, o_sl)
for(scn in names(scen)){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
  smh <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
  sday <- merge(sim[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=.(id_plot,date,mo)], mday[,.(date,Tmax_macro)], by="date")
  d_dT <- sday[, .(value=mean(Tmax_sim-Tmax_macro,na.rm=TRUE)), by=.(id_plot,mo)][, `:=`(source=scen[[scn]], metric="dTmax")]
  d_sl <- smh[, .(value=slp(Tsim,Tmacro)), by=.(id_plot,mo)][, `:=`(source=scen[[scn]], metric="slope")]
  per <- c(per, list(d_dT, d_sl))
}
D <- rbindlist(per, fill=TRUE)
mo_present <- intersect(c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"), unique(D$mo))
D[, mo := factor(mo, levels=mo_present)]
short <- c("Observed (HOBO)"="Obs","LiDAR full"="LiDAR full","LiDAR d_opt"="LiDAR d_opt",
           "S2 ATBD"="ATBD","S2 opt"="opt","S2 ATBD ×ratio→full"="ATBD→full","S2 opt ×ratio→d_opt"="opt→d_opt")
D[, sx := factor(short[as.character(source)], levels=unname(short))]
D[, typ := fcase(source=="Observed (HOBO)","Observed", grepl("LiDAR",source),"LiDAR", default="Sentinel-2")]
fwrite(D, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4q_hotdays_boxplot_data.csv")
pal <- c(Observed="grey40", LiDAR=unname(PAL_SENSOR["LiDAR"]), "Sentinel-2"=unname(PAL_SENSOR["Sentinel-2"]))
mk <- function(met, ylab, href){
  ggplot(D[metric==met], aes(sx, value, fill=typ)) +
    { if(!is.na(href)) geom_hline(yintercept=href, linetype=3, colour="grey50") } +
    geom_boxplot(outlier.size=0.4, linewidth=0.3) +
    facet_wrap(~mo, ncol=length(mo_present)) +
    scale_fill_manual(values=pal, name=NULL) +
    labs(x=NULL, y=ylab, subtitle=sprintf("10%% hottest days of Apr–Nov (macro Tmax ≥ %.1f °C, n = %d days)", thr, length(hot))) +
    theme_article(9) + theme(legend.position="top", axis.text.x=element_text(angle=45,hjust=1,size=7))
}
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
w <- max(7, 2.2*length(mo_present))
ggsave_article(file.path(outdir,"Fig_boxplot_hotdays_dTmax"), mk("dTmax","ΔTmax (sub-canopy − macro, °C)", 0), w, 4.5)
ggsave_article(file.path(outdir,"Fig_boxplot_hotdays_slope"), mk("slope","slope-and-equilibrium slope", 1), w, 4.5)
cat("DONE -> Table4q + 2 hot-day boxplot figures\n")
