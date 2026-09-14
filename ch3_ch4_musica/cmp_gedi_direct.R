# ==============================================================================
# cmp_gedi_direct.R — Chapter 4: score STATIC_GEDI_DIRECT (magnitude = mean PAI
# of the 3 nearest power footprints, no regression) against the leaderboard
# heads, STRATIFIED BY DISTANCE to the nearest footprint (the operational-domain
# question: how close must GEDI shoot for direct use to pay off?).
# Summer JJAS, genuine-53 windcorr, v3.2.3-iter.
# Writes NC_Full/manuscripts/ch4/tables/Table11_c4_direct_musica.csv.
#   Rscript cmp_gedi_direct.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")

FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MA<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
sdo_d<-sd(OBS$do)

d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
getpp<-function(scn){fs<-list.files(file.path(d3,scn),pattern="\\.nc$",full.names=TRUE)
 m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr,Tm)],by="hr")
  sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
 merge(m,OBS,by="id_plot")}

SC<-c("STATIC_GEDI_DIRECT","STATIC_GEDI_RF","NAIVE_S2_FORMSH","STATIC_S2_ATBD","STATIC_ALS")
PP<-setNames(lapply(SC,getpp),SC)

dd<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_gedi_direct.rds")))
dd[,id_plot:=as.character(id_plot)]
strata<-list(`near_le150m`=dd[dist_fp<=150,id_plot],
             `mid_le300m` =dd[dist_fp<=300,id_plot],
             `far_gt300m` =dd[dist_fp>300,id_plot],
             `pooled`     =dd$id_plot)
met<-function(d) data.table(n=nrow(d),
  dT_R2=cor(d$ds,d$do)^2, dT_RMSE=sqrt(mean((d$ds-d$do)^2)), dT_bias=mean(d$ds-d$do),
  sl_R2=cor(d$ss,d$so,use="complete.obs")^2)
res<-rbindlist(lapply(names(PP),function(s) rbindlist(lapply(names(strata),function(st)
  cbind(scenario=s,stratum=st,met(PP[[s]][id_plot%in%strata[[st]]]))))))
num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
for(st in names(strata)){cat("\n===",st,"(n=",length(strata[[st]]),") ===\n")
  print(res[stratum==st][order(-dT_R2),.(scenario,dT_R2,dT_RMSE,dT_bias,sl_R2)])}
fwrite(res,file.path(TAB,"Table11_c4_direct_musica.csv"))
cat("\nwrote",file.path(TAB,"Table11_c4_direct_musica.csv"),"\n")
