# ==============================================================================
# Chapter 3 — comprehensive 53-HOBO validation: ΔTmax & slope-and-equilibrium, all
# at 1 m, WITH vs WITHOUT the -2h sim time shift, over windows {spring, summer,
# autumn, whole=spring+summer+autumn}; the whole & summer windows also in a hottest
# -days (p90) variant. Unified from the 1 m hourly sim series (shift = parameter).
# Run from z_Example root:  Rscript extend_full.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full (dyn)", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt (dyn)",
          STATIC_ALS="LiDAR full (static)", STATIC_ALS_DOPT="LiDAR d_opt (static)")
WIN <- as.Date(c("2021-04-15","2021-11-15"))
ds_all <- seq(WIN[1], WIN[2], by="day")

# macro: hourly + daily Tmax over the whole window
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"), Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, `:=`(date=as.Date(time), hr=floor_date(time,"hour"))]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=date]

# HOBO hourly (no shift; it is the observation)
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
hb[, `:=`(date=as.Date(time), hr=floor_date(time,"hour"))]
hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
hb_day <- hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date)]
hb_day <- merge(hb_day, mday, by="date")[, dT_obs := Tmax_obs - Tmax_macro][]

# period day-sets
P <- list(
  spring   = list(days=ds_all[ds_all<as.Date("2021-06-01")], hot=FALSE),
  summer   = list(days=ds_all[ds_all>=as.Date("2021-06-01") & ds_all<=as.Date("2021-09-30")], hot=FALSE),
  autumn   = list(days=ds_all[ds_all>as.Date("2021-09-30")], hot=FALSE),
  whole    = list(days=ds_all, hot=FALSE))
hot_days <- function(days){ q<-quantile(mday[date %in% days]$Tmax_macro,0.90,na.rm=TRUE); mday[date %in% days & Tmax_macro>=q]$date }
P[["summer_hot"]] <- list(days=hot_days(P$summer$days), hot=TRUE)
P[["whole_hot"]]  <- list(days=hot_days(P$whole$days),  hot=TRUE)

slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% ds_all, .(time, Tsim=Tair_sim)]; r }

# observed per-plot metrics by period (no shift)
obs_metrics <- function(days){
  d <- hb_day[date %in% days, .(do=mean(dT_obs,na.rm=TRUE)), by=id_plot]
  s <- hb[date %in% days, .(so=slp(Tobs,Tmacro)), by=id_plot]
  merge(d,s,by="id_plot") }
OBS <- lapply(P, function(p) obs_metrics(p$days))

res <- list()
for(scn in names(scen)){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
  for(shift in c(0L,2L)){
    sm <- copy(sim); sm[, t2:=time-hours(shift)]; sm[, `:=`(date=as.Date(t2), hr=floor_date(t2,"hour"))]
    sm2 <- merge(sm, macro_h[,.(hr,Tmacro)], by="hr")                  # hourly for slope
    sday <- sm[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=.(id_plot,date)]
    sday <- merge(sday, mday, by="date")[, dT_sim:=Tmax_sim-Tmax_macro][]
    for(pn in names(P)){ days<-P[[pn]]$days
      dT <- sday[date %in% days, .(ds=mean(dT_sim,na.rm=TRUE)), by=id_plot]
      sl <- sm2[date %in% days, .(ss=slp(Tsim,Tmacro)), by=id_plot]
      d <- merge(merge(dT,sl,by="id_plot"), OBS[[pn]], by="id_plot")
      dd<-d[is.finite(ds)&is.finite(do)]; sd<-d[is.finite(ss)&is.finite(so)]
      res[[length(res)+1]] <- data.table(scenario=scen[[scn]], period=pn,
        shift=ifelse(shift==0,"no shift","-2h shift"), n=nrow(dd),
        dT_r2=cor(dd$ds,dd$do)^2, dT_bias=mean(dd$ds-dd$do),
        sl_r2=cor(sd$ss,sd$so)^2, sl_bias=mean(sd$ss-sd$so)) }
  }
}
R <- rbindlist(res); num<-c("dT_r2","dT_bias","sl_r2","sl_bias"); R[,(num):=lapply(.SD,round,3),.SDcols=num]
R[, period:=factor(period, levels=c("spring","summer","autumn","whole","summer_hot","whole_hot"))]
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4j_full_shift_periods.csv")
cat("\n=== ΔTmax & slope R² | shift vs no-shift | periods ===\n")
print(R[order(period,shift,-dT_r2), .(scenario,period,shift,dT_r2,sl_r2,dT_bias,sl_bias)], nrow=200)
cat("\nDONE -> Table4j_full_shift_periods.csv\n")
