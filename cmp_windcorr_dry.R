# ==============================================================================
# Compare wind-corrected vs uniform-wind Ch3 runs on the 3 reference scenarios.
# Per scenario: ΔTmax bias vs obs, RMSE, between-plot R², SD-recovery, slope.
#   Rscript cmp_windcorr_dry.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
UNI<-file.path(CFG_C3$out_dir,"nc_genuine53"); WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr_dry")
SCEN<-c("CONST_ALS","STATIC_ALS","DYN_S2_ATBD")

# macro (observed forcing air temp) from the pblh forcing used by the runs
ncf<-nc_open(CFG_C3$forcing_file %>% {ifelse(file.exists(.),.,"in_files/musica_in_Blois_pblh.nc")})
FORC<-if(file.exists("in_files/musica_in_Blois_pblh.nc"))"in_files/musica_in_Blois_pblh.nc" else CFG_C3$forcing_file
nc_close(ncf)
ncf<-nc_open(FORC); tu<-ncatt_get(ncf,"time","units")$value; th<-ncvar_get(ncf,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds]; mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA

# observed
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],
           hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")

getpp<-function(scn,dir){fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}

summ<-function(scn,dir,tag){P<-getpp(scn,dir); if(is.null(P)||!nrow(P))return(NULL)
  M<-merge(OBS,P,by="id_plot")
  data.table(scenario=scn,wind=tag,n=nrow(M),
    bias_dT=mean(M$ds-M$do,na.rm=TRUE), rmse_dT=sqrt(mean((M$ds-M$do)^2,na.rm=TRUE)),
    R2_dT=cor(M$ds,M$do,use="complete.obs")^2, SDrec=sd(M$ds,na.rm=TRUE)/sd(M$do,na.rm=TRUE),
    slope_bias=mean(M$ss-M$so,na.rm=TRUE), R2_slope=cor(M$ss,M$so,use="complete.obs")^2)}

out<-rbindlist(lapply(SCEN,function(s) rbind(summ(s,UNI,"uniform"),summ(s,WC,"windcorr"))),fill=TRUE)
cat(sprintf("Observed: SD(ΔTmax)=%.2f  range[%.1f,%.1f]  SD(slope)=%.3f\n\n",
    sd(OBS$do),min(OBS$do),max(OBS$do),sd(OBS$so)))
print(out[order(scenario,wind)],digits=3)
fwrite(out,"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/windcorr_dry_compare.csv")
cat("\nWROTE chapter3_S2_LAI/tables/windcorr_dry_compare.csv\n")
