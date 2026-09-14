# ==============================================================================
# Chapter 3 — month-by-month 53-HOBO validation. Per scenario and per calendar
# month (Apr 15 -> Nov 15, 2021): ΔTmax & slope-and-equilibrium agreement across
# the 53 plots as R², RMSE, MAE and bias. Also the OBSERVED buffering regime per
# month: how many of the 53 plots buffer (obs slope<1) vs amplify (obs slope>=1),
# and the mean observed slope / ΔTmax. Metrics at 1 m, no time shift.
# Run from z_Example root:  Rscript monthly_metrics.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt",
          DYN_ALS_S2TIMING_NM="S2t·ALS/ATBD", DYN_ALSoOPT_S2TIMING_NM="S2t·ALS/opt",
          STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt")
WIN <- as.Date(c("2021-04-15","2021-11-15")); ds_all <- seq(WIN[1], WIN[2], by="day")
mlab <- function(d) format(d, "%Y-%m")

# macro hourly + daily Tmax
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, `:=`(hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]; mday[, mo:=mlab(date)]

# HOBO (obs, no shift)
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
hb[, `:=`(hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]; hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
hb_day <- merge(hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date=as.Date(time))], mday[,.(date,Tmax_macro)], by="date")
hb_day[, mo:=mlab(date)]
OBS <- merge(hb_day[, .(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)), by=.(id_plot,mo)],
             hb[, .(so=slp(Tobs,Tmacro)), by=.(id_plot,mo)], by=c("id_plot","mo"))

# observed buffering regime per month
regime <- OBS[!is.na(so), .(n=.N, n_buff=sum(so<1), n_amp=sum(so>=1),
                            mean_slope=mean(so), mean_dT=mean(do)), by=mo][order(mo)]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% ds_all, .(time, Tsim=Tair_sim)]; r }
agg <- function(s,o){ d<-data.table(s=s,o=o)[is.finite(s)&is.finite(o)]
  list(n=nrow(d), r2=if(nrow(d)>2) cor(d$s,d$o)^2 else NA_real_,
       rmse=sqrt(mean((d$s-d$o)^2)), mae=mean(abs(d$s-d$o)), bias=mean(d$s-d$o)) }

res <- list()
for(scn in names(scen)){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL)
    s[, `:=`(id_plot=id, hr=floor_date(time,"hour"), date=as.Date(time), mo=mlab(as.Date(time)))][] }), fill=TRUE)
  smh <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
  sday <- merge(sim[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=.(id_plot,date,mo)], mday[,.(date,Tmax_macro)], by="date")
  for(m in regime$mo){
    dT <- sday[mo==m, .(ds=mean(Tmax_sim-Tmax_macro,na.rm=TRUE)), by=id_plot]
    sl <- smh[mo==m, .(ss=slp(Tsim,Tmacro)), by=id_plot]
    d  <- merge(merge(dT,sl,by="id_plot"), OBS[mo==m], by="id_plot")
    aT<-agg(d$ds,d$do); aS<-agg(d$ss,d$so)
    res[[length(res)+1]] <- data.table(scenario=scen[[scn]], mo=m, n=aT$n,
      dT_r2=aT$r2, dT_rmse=aT$rmse, dT_mae=aT$mae, dT_bias=aT$bias,
      sl_r2=aS$r2, sl_rmse=aS$rmse, sl_mae=aS$mae, sl_bias=aS$bias) }
}
R <- rbindlist(res); numc <- setdiff(names(R), c("scenario","mo","n"))
R[, (numc):=lapply(.SD, round, 3), .SDcols=numc]
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4l_monthly_metrics.csv")
fwrite(regime, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4m_monthly_obs_regime.csv")

mo_levels <- regime$mo   # precompute to avoid column-name shadowing inside data.table j
cat("\n=== observed buffering regime per month (53 plots) ===\n"); print(regime)
cat("\n=== monthly metrics (ΔTmax) ===\n")
print(dcast(R, mo ~ scenario, value.var="dT_r2")[, .SD, .SDcols=patterns("mo|LiDAR full|S2 ATBD \\(dyn\\)")])

# figure: ΔTmax R²/RMSE/MAE + slope R² across months, per scenario
ml <- melt(R, id.vars=c("scenario","mo"),
           measure.vars=c("dT_r2","dT_rmse","dT_mae","sl_r2"), variable.name="metric", value.name="val")
ml[, metric:=factor(metric, levels=c("dT_r2","dT_rmse","dT_mae","sl_r2"),
                    labels=c("(a) ΔTmax R²","(b) ΔTmax RMSE (°C)","(c) ΔTmax MAE (°C)","(d) slope R²"))]
ml[, typ:=fcase(grepl("S2t·",scenario),"Fusion",
                grepl("LiDAR",scenario),"LiDAR", default="Sentinel-2")]
ml[, mo:=factor(mo, levels=mo_levels)]
g <- ggplot(ml, aes(mo, val, group=scenario, colour=typ)) +
  geom_line(aes(linetype=scenario), linewidth=0.7) + geom_point(size=1.3) +
  facet_wrap(~metric, scales="free_y", ncol=2) +
  scale_colour_manual(values=PAL_SENSOR, name=NULL) +
  scale_linetype_manual(values=c(1,2,3,4,1,2,3,4), guide=guide_legend(ncol=3)) +
  labs(x=NULL, y=NULL, linetype=NULL) + theme_article(9) +
  theme(legend.position="bottom", axis.text.x=element_text(angle=30,hjust=1))
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_monthly_metrics"), g, 6.9, 6.4)

# figure: observed buffering composition per month
rg <- melt(regime, id.vars="mo", measure.vars=c("n_buff","n_amp"), variable.name="regime", value.name="count")
rg[, regime:=factor(regime, levels=c("n_amp","n_buff"), labels=c("amplifying (slope ≥ 1)","buffering (slope < 1)"))]
rg[, mo:=factor(mo, levels=mo_levels)]
g2 <- ggplot(rg, aes(mo, count, fill=regime)) + geom_col(width=0.7) +
  scale_fill_manual(values=c("buffering (slope < 1)"=unname(PAL_GRP["buffering"]),
                             "amplifying (slope ≥ 1)"=unname(PAL_GRP["amplifying"])), name=NULL) +
  labs(x=NULL, y="number of plots (of 53)") + theme_article(10) +
  theme(legend.position="top", axis.text.x=element_text(angle=30,hjust=1))
ggsave_article(file.path(outdir,"Fig_monthly_obs_regime"), g2, 6.9, 3.4)
cat("\nDONE -> Table4l_monthly_metrics.csv, Table4m_monthly_obs_regime.csv, 2 figs\n")
