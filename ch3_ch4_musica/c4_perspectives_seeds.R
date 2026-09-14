# ==============================================================================
# c4_perspectives_seeds.R — seed analyses for TWO additional General-Discussion
# perspectives (beyond GEDI), using only data already in hand:
#
#  A. BUFFERING UNDER EXTREME HEAT (2022). HOBO covers the record June/July 2022
#     heatwaves and the existing genuine-53 windcorr sims run to 2022-07-30.
#     (i) Does the observed LAI-buffering coupling hold on 2022 hot days vs
#         summer 2021?  (ii) Do the thesis forcing products (Ch3 fusion, Ch4
#         GEDI-corrected fusion) still validate in an extreme year?
#
#  B. MICROCLIMATE -> PHENOLOGY LOOP (Wu et al. 2024 at stand scale). Per-plot
#     S2 green-up date (TRS50 up-crossing of the per-plot ATBD series) vs the
#     plot's own observed spring microclimate (HOBO mean daily Tmax,
#     Mar 15 - Apr 30, 2021), with the LAI confound checked.
#
# Outputs: NC_Full/manuscripts/ch4/tables/Table23*, printout.
#   Rscript c4_perspectives_seeds.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df[,id_plot:=as.character(id_plot)];df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))

FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MAall<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
hball<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hball[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hball<-hball[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA

obs_window<-function(ds,hot_p=NA){
  MA<-MAall[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")]
  mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
  if(is.finite(hot_p)){hd<-mday[Tmx>=quantile(Tmx,hot_p)];ds<-hd$date
    MA<-MA[as.Date(time)%in%ds];mday<-mday[date%in%ds]}
  hb<-hball[as.Date(time)%in%ds,.(id_plot,Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
  hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
  merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE),n_days=uniqueN(date)),by=id_plot],
        hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
}

# ---- A(i): observed coupling, 2021 summer vs 2022 heat ------------------------
w21<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
w22<-seq(as.Date("2022-06-01"),as.Date("2022-07-29"),by="day")
res_obs<-rbindlist(lapply(list(
    list(lab="2021 JJAS (all days)",     ds=w21, hp=NA),
    list(lab="2021 JJAS (hot p90)",      ds=w21, hp=0.9),
    list(lab="2022 Jun-Jul (all days)",  ds=w22, hp=NA),
    list(lab="2022 Jun-Jul (hot p90)",   ds=w22, hp=0.9)),
  function(cfg){o<-merge(obs_window(cfg$ds,cfg$hp),df[,.(id_plot,LAI_ALS)],by="id_plot")
    data.table(window=cfg$lab,n_plots=nrow(o),
      mean_dTmax=round(mean(o$do),2), sd_dTmax=round(sd(o$do),2),
      r_dTmax_LAI=round(cor(o$do,o$LAI_ALS),3),
      r_slope_LAI=round(cor(o$so,o$LAI_ALS,use="complete.obs"),3),
      mean_slope=round(mean(o$so,na.rm=TRUE),3))}))
cat("=== A(i) observed LAI-buffering coupling, normal vs extreme ===\n")
print(res_obs)
fwrite(res_obs,file.path(TAB,"Table23_c4_heat2022_obs.csv"))

# ---- A(ii): scenario validation in the 2022 window -----------------------------
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
getpp<-function(scn,ds,OBS){fs<-list.files(file.path(d3,scn),pattern="\\.nc$",full.names=TRUE)
 MA<-MAall[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")]
 mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
 m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))]
  if(!nrow(s))return(NULL)
  sm<-merge(s,MA[,.(hr,Tm)],by="hr")
  sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date")
  data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
 merge(m,OBS,by="id_plot")}
OBS22<-obs_window(w22)
res22<-rbindlist(lapply(c("FUSION_GEDIMAX","FUSION_H","STATIC_ALS","STATIC_S2_ATBD","DYN_ALS_GEDIMAX"),
  function(s){d<-getpp(s,w22,OBS22)
   data.table(scenario=s,n=nrow(d),
     dT_R2=round(cor(d$ds,d$do)^2,3), dT_RMSE=round(sqrt(mean((d$ds-d$do)^2)),2),
     dT_bias=round(mean(d$ds-d$do),2), sl_R2=round(cor(d$ss,d$so,use="complete.obs")^2,3))}))
cat("\n=== A(ii) scenario validation, Jun 1 - Jul 29 2022 (heatwave summer) ===\n")
print(res22[order(-dT_R2)])
fwrite(res22,file.path(TAB,"Table23b_c4_heat2022_scenarios.csv"))

# ---- B: per-plot green-up vs spring microclimate -------------------------------
tb<-readRDS(file.path(CFG_C3$out_dir,"lai_prep","ts_by_plot.rds"))
greenup<-rbindlist(lapply(names(tb$atbd),function(p){a<-tb$atbd[[p]]
  if(is.null(a)||!nrow(a))return(NULL)
  half<-min(a$lai)+diff(range(a$lai))/2
  up<-which(a$lai>=half)
  data.table(pid=p,sos_s2=a$doy[up[1]])}))
spring<-seq(as.Date("2021-03-15"),as.Date("2021-04-30"),by="day")
hb<-hball[as.Date(time)%in%spring,.(id_plot,Tobs=t_hobo,date=as.Date(time))]
tspring<-hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)][,.(t_spring=mean(Tmx_o)),by=id_plot]
B<-merge(merge(greenup,df[,.(pid,id_plot,LAI_ALS,Hmax)],by="pid"),tspring,by="id_plot")
cat("\n=== B microclimate -> phenology loop (n =",nrow(B),") ===\n")
cat("cor(SOS_S2, spring Tmax_micro):",round(cor(B$sos_s2,B$t_spring),3),"\n")
cat("cor(SOS_S2, LAI_ALS):",round(cor(B$sos_s2,B$LAI_ALS),3),
    " | cor(spring T, LAI_ALS):",round(cor(B$t_spring,B$LAI_ALS),3),"\n")
pc<-residuals(lm(sos_s2~LAI_ALS,B));pt<-residuals(lm(t_spring~LAI_ALS,B))
cat("partial cor(SOS, springT | LAI_ALS):",round(cor(pc,pt),3),"\n")
fwrite(B,file.path(TAB,"Table23c_c4_pheno_loop.csv"))
cat("\nDONE\n")
