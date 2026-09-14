# ==============================================================================
# cmp_gedi_windcorr.R — Chapter 4 Part B extraction: score the GEDI-anchored
# MuSICA scenarios against HOBO (summer JJAS, genuine-53 windcorr, v3.2.3-iter)
# next to the Ch3 leaderboard heads, pooled AND stratified at the Ch2 crossover
# (LAI_ALS = 3.86). FUSION_GEDI is assembled post-hoc: dense plots (FORMS_H >=
# 17.8 m, the FigK operational gate) from STATIC_GEDI_RF nc, open plots from
# STATIC_S2_ATBD nc — no extra simulations.
# Writes NC_Full/manuscripts/ch4/tables/Table8*.csv.
#   Rscript cmp_gedi_windcorr.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))

# --- obs microclimate + macro (identical to c3_figset_data_wc.R) --------------
FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MA<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
sdo_d<-sd(OBS$do); sdo_s<-sd(OBS$so,na.rm=TRUE)

# --- per-plot sim metrics per scenario ----------------------------------------
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
getpp<-function(scn,dir=d3){fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
 m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr,Tm)],by="hr")
  sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
 merge(m,OBS,by="id_plot")}

SC<-c("STATIC_GEDI_RF","STATIC_GEDI_RATIO","NAIVE_S2_FORMSH","STATIC_S2_ATBD",
      "STATIC_ALS","STATIC_ALS_DOPT")
PP<-setNames(lapply(SC,getpp),SC)

# FUSION_GEDI: dense (FORMS_H >= 17.8) <- GEDI RF nc ; open <- raw-S2 nc
hthr<-17.8
dense_ids<-df[FORMS_H>=hthr, id_plot]
cat(sprintf("FUSION_GEDI gate FORMS_H>=%.1f m: %d dense / %d open\n",
            hthr, length(dense_ids), nrow(df)-length(dense_ids)))
PP$FUSION_GEDI<-rbind(PP$STATIC_GEDI_RF[id_plot%in%dense_ids],
                      PP$STATIC_S2_ATBD[!id_plot%in%dense_ids])

# --- pooled + stratified metrics ----------------------------------------------
met<-function(d) data.table(n=nrow(d),
  dT_R2=cor(d$ds,d$do)^2, dT_RMSE=sqrt(mean((d$ds-d$do)^2)), dT_bias=mean(d$ds-d$do),
  dT_SDrec=sd(d$ds)/sdo_d,
  sl_R2=cor(d$ss,d$so,use="complete.obs")^2, sl_bias=mean(d$ss-d$so,na.rm=TRUE))
xover<-3.86; dense_lai<-df[LAI_ALS>=xover, id_plot]
res<-rbindlist(lapply(names(PP),function(s){
  d<-PP[[s]]
  rbind(cbind(scenario=s,stratum="pooled",met(d)),
        cbind(scenario=s,stratum="dense",met(d[id_plot%in%dense_lai])),
        cbind(scenario=s,stratum="open", met(d[!id_plot%in%dense_lai])))}))
num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
cat("\n=== POOLED (summer JJAS, n=53) ===\n")
print(dcast(res[stratum=="pooled"],scenario~.,value.var=c("dT_R2","dT_RMSE","dT_bias","dT_SDrec","sl_R2"))[order(-dT_R2)])
cat("\n=== DENSE (LAI_ALS >= 3.86) ===\n")
print(dcast(res[stratum=="dense"],scenario~.,value.var=c("dT_R2","dT_RMSE","dT_bias","sl_R2"))[order(-dT_R2)])
cat("\n=== OPEN (LAI_ALS < 3.86) ===\n")
print(dcast(res[stratum=="open"],scenario~.,value.var=c("dT_R2","dT_RMSE","dT_bias","sl_R2"))[order(-dT_R2)])
fwrite(res,file.path(TAB,"Table8_c4_musica_gedi.csv"))
cat("\nwrote",file.path(TAB,"Table8_c4_musica_gedi.csv"),"\n")
