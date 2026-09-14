# ==============================================================================
# Chapter 3 — stratified 53-HOBO validation: does any plot stratum let Sentinel-2
# close the gap with LiDAR? We split the 53 plots by (a) LiDAR LAI, (b) canopy
# height Hmax, (c) observed buffering (ΔTmax_obs), median split, and report per
# stratum the bias & RMSE (robust to range restriction) plus R² (in support, with
# the caveat that R² shrinks mechanically in a narrow stratum). Summer, 1 m, no
# shift. Run from z_Example root:  Rscript strat_metrics.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt",
          STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt")
ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")   # summer

# covariates
cov <- as.data.table(readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds"))[
  , .(id_plot=as.character(id_plot), LAI_ALS, Hmax, fCover)]

# macro hourly + daily Tmax
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, hr:=floor_date(time,"hour")]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]

# HOBO (obs, no shift)
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
hb[, hr:=floor_date(time,"hour")]; hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
hb_day <- merge(hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date=as.Date(time))], mday, by="date")
OBS <- merge(hb_day[, .(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)), by=id_plot],
             hb[, .(so=slp(Tobs,Tmacro)), by=id_plot], by="id_plot")

# sim per plot (1 m, no shift)
read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% ds_all, .(time, Tsim=Tair_sim)]; r }

per_plot <- rbindlist(lapply(names(scen), function(scn){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL)
    sm <- merge(s[, .(time,Tsim,hr=floor_date(time,"hour"),date=as.Date(time))], macro_h[,.(hr,Tmacro)], by="hr")
    sd <- merge(s[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=.(date=as.Date(time))], mday, by="date")
    data.table(scenario=scen[[scn]], id_plot=id,
               ds=mean(sd$Tmax_sim-sd$Tmax_macro,na.rm=TRUE), ss=slp(sm$Tsim,sm$Tmacro)) }), fill=TRUE)
}), fill=TRUE)
per_plot <- merge(merge(per_plot, OBS, by="id_plot"), cov, by="id_plot")

# slope-and-equilibrium class from OBSERVED slope: <1 buffering, >=1 amplifying
po <- unique(per_plot[, .(id_plot, so)])
cat(sprintf("\n=== observed slope-and-equilibrium class (summer, n=%d) ===\n", nrow(po)))
cat(sprintf("   slope < 1 (buffering): %d plots | slope >= 1 (amplifying): %d plots\n",
            sum(po$so < 1, na.rm=TRUE), sum(po$so >= 1, na.rm=TRUE)))
cat(sprintf("   observed slope range: %.2f to %.2f (median %.2f)\n",
            min(po$so,na.rm=TRUE), max(po$so,na.rm=TRUE), median(po$so,na.rm=TRUE)))
fc <- unique(per_plot[, .(id_plot, fCover, LAI_ALS)])
cat(sprintf("   fCover range: %.3f to %.3f (median %.3f) | LAI_ALS range %.2f to %.2f (median %.2f)\n",
            min(fc$fCover), max(fc$fCover), median(fc$fCover),
            min(fc$LAI_ALS), max(fc$LAI_ALS), median(fc$LAI_ALS)))

# stratify: median split on each covariate + a fixed slope=1 split on observed slope
strat_vars <- list(LAI_ALS="LiDAR LAI", fCover="fractional cover", Hmax="canopy height", do="observed ΔTmax (buffering)")
out <- list()
for(v in names(strat_vars)){
  med <- median(per_plot[[v]][!duplicated(per_plot$id_plot)], na.rm=TRUE)
  pp <- copy(per_plot); pp[, grp:=ifelse(get(v) < med, "low", "high")]
  agg <- pp[, .(n=.N,
                dT_bias=mean(ds-do), dT_rmse=sqrt(mean((ds-do)^2)),
                sl_bias=mean(ss-so,na.rm=TRUE), sl_rmse=sqrt(mean((ss-so)^2,na.rm=TRUE)),
                dT_r2=cor(ds,do)^2, sl_r2=cor(ss,so,use="complete.obs")^2),
            by=.(scenario,grp)]
  agg[, `:=`(strat=strat_vars[[v]], var=v)]
  out[[v]] <- agg
}
# slope-and-equilibrium split (observed slope vs 1)
pp <- copy(per_plot); pp[, grp:=ifelse(so < 1, "buffering (slope<1)", "amplifying (slope>=1)")]
aggS <- pp[, .(n=.N,
               dT_bias=mean(ds-do), dT_rmse=sqrt(mean((ds-do)^2)),
               sl_bias=mean(ss-so,na.rm=TRUE), sl_rmse=sqrt(mean((ss-so)^2,na.rm=TRUE)),
               dT_r2=cor(ds,do)^2, sl_r2=cor(ss,so,use="complete.obs")^2),
           by=.(scenario,grp)]
aggS[, `:=`(strat="slope-and-equilibrium class", var="so_class")]
out[["so_class"]] <- aggS
R <- rbindlist(out)
numc <- c("dT_bias","dT_rmse","sl_bias","sl_rmse","dT_r2","sl_r2")
R[, (numc):=lapply(.SD, round, 3), .SDcols=numc]
setcolorder(R, c("strat","grp","scenario","n",numc))
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4k_stratified.csv")

ord <- c("LiDAR full","S2 ATBD (dyn)","S2 opt (dyn)","S2 ATBD ×ratio→full","S2 opt ×ratio→d_opt","LiDAR d_opt")
for(v in c(names(strat_vars),"so_class")){
  lab <- if(v=="so_class") "slope-and-equilibrium class (obs slope vs 1)" else paste0(strat_vars[[v]]," (median split)")
  cat(sprintf("\n=== stratify by %s ===\n", lab))
  sub <- R[var==v]; sub[, scenario:=factor(scenario, levels=ord)]
  print(sub[order(grp,scenario), .(grp,scenario,n,dT_bias,dT_rmse,dT_r2,sl_bias,sl_rmse,sl_r2)])
}
cat("\nDONE -> Table4k_stratified.csv\n")
