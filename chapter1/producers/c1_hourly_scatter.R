# ==============================================================================
# Density scatter of ALL hourly 1 m air temps: MuSICA vs HOBO, 53 loggers pooled.
# TWO panels: (a) all summer days, (b) hottest 10% of days. 1:1, R2/RMSE/MAE.
#
# Unlike the DeltaTmax figures this compares the raw hourly series, so it does not use
# R/dtmax_convention.R; it applies OBS_CLOCK_OFFSET_H directly to the observation labels
# (forcing and simulations are on solar time = UTC+1, the loggers on UTC), which is the
# same alignment expressed at the hourly level.
#
# Reads : in_files/FR-Blo_2021_v2.nc                                  (forcing; hot-day selection)
#         in_files/Blois_data_temperature.csv                         (CFG$hobo_temp_csv)
#         out_files/musica_hobo_native20/1111/musica_out_HOBO_<id>.nc (full-model runs)
# Writes: outputs/figures_chap1/Fig_hourly_temp_scatter.{png,pdf}
#   Rscript c1_hourly_scatter.R
# ==============================================================================
# Density scatter of ALL hourly 1 m air temps: MuSICA vs HOBO, 53 loggers pooled.
# TWO panels: (a) all summer days, (b) hottest 10% of days. 1:1, R²/RMSE/MAE.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(ggplot2);library(patchwork);library(viridisLite)
  src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("pipeline/00_config.R");source("scripts/_article_style.R")})
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; SIMD<-"out_files/musica_hobo_native20/1111"; FORC<-"in_files/FR-Blo_2021_v2.nc"
# hot days = top 10% by station daily Tmax
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC");mac<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macD<-mac[as.Date(time)%in%ds,.(Tmx=max(Tm)),by=.(date=as.Date(time))];hot<-macD[Tmx>=quantile(Tmx,0.90),date]
hobo<-as.data.table(read.csv(CFG$hobo_temp_csv));hobo[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hobo<-hobo[position_sensor=="a"&as.Date(datetime)%in%ds&!id_plot%in%CFG$ids_to_remove]
OBS<-hobo[,.(Tobs=mean(t_hobo,na.rm=T)),by=.(id_plot,time=floor_date(datetime,"hour"))]
# CLOCK ALIGNMENT: the forcing (and therefore every simulation) is on solar time
# = UTC+1, the loggers are on UTC, so forcing labels run one hour ahead. Shift the
# observation labels forward to pair the same physical instant before merging.
OBS[,time:=time-lubridate::hours(OBS_CLOCK_OFFSET_H)]
#' Hourly sub-canopy temperature interpolated to the fixed height Z from one run
#' @param p path to a MuSICA output NetCDF
#' @return data.table(time, Tsim) over the summer window, or NULL if the file is missing,
#'   truncated, unreadable, or lacks the profile variables
mic<-function(p){if(!file.exists(p)||file.size(p)<1e5)return(NULL);nc<-try(nc_open(p),silent=T);if(inherits(nc,"try-error"))return(NULL);on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var)))return(NULL)
  tu<-ncatt_get(nc,"time","units")$value;t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC");th<-ncvar_get(nc,"time");Tk<-ncvar_get(nc,"Tair_z");rh<-ncvar_get(nc,"relative_height");vh<-median(ncvar_get(nc,"veget_height_top"),na.rm=T);zl<-rh*vh
  if(Z<=zl[1]){il<-1;ih<-1;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  data.table(time=floor_date(t0+dhours(th),"hour"),Tsim=((1-w)*Tk[il,]+w*Tk[ih,])-273.15)[as.Date(time)%in%ds]}
fs<-list.files(SIMD,"musica_out_HOBO_.*\\.nc$");SIM<-rbindlist(lapply(fs,function(f){d<-mic(file.path(SIMD,f));if(is.null(d))return(NULL);d[,id_plot:=sub("musica_out_HOBO_(.*)\\.nc","\\1",f)]}))
D<-merge(OBS,SIM,by=c("id_plot","time"))[is.finite(Tobs)&is.finite(Tsim)];D[,hotday:=as.Date(time)%in%hot]
#' One hex-density 1:1 panel of simulated against observed hourly temperature
#' @param dd paired hourly table
#' @param ttl panel title
#' @return a ggplot object annotated with R-squared, RMSE, MAE and n
mkp<-function(dd,ttl){r2<-cor(dd$Tsim,dd$Tobs)^2;rmse<-sqrt(mean((dd$Tsim-dd$Tobs)^2));mae<-mean(abs(dd$Tsim-dd$Tobs));lim<-range(c(dd$Tobs,dd$Tsim))
  lab<-sprintf("italic(R)^2=='%.2f'*'  RMSE='*'%.2f'*'  MAE='*'%.2f'*'  n='*'%s'",r2,rmse,mae,format(nrow(dd),big.mark=","))
  ggplot(dd,aes(Tobs,Tsim))+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey40")+geom_hex(bins=70)+
    scale_fill_viridis_c(trans="log10",name="count",option="magma",direction=-1)+
    annotate("text",x=lim[1],y=lim[2],hjust=0,vjust=1,parse=TRUE,label=lab,size=3.3)+
    coord_equal(xlim=lim,ylim=lim)+labs(x="observed T (°C, sub-canopy plots)",y="simulated T (°C, MuSICA)",title=ttl)+theme_article(11)+theme(plot.title=element_text(size=11,hjust=0.5,face="plain"))+legend_corner(0.99,0.02)}
pa<-mkp(D,"all summer days");pb<-mkp(D[hotday==TRUE],"hottest 10% of days")
fig<-(pa|pb)+plot_annotation(tag_levels="a")&theme(plot.tag=element_text(face="bold",size=14))
ggsave_article("outputs/figures_chap1/Fig_hourly_temp_scatter",fig,10,5.2)
cat(sprintf("ALL: R2=%.3f RMSE=%.2f (n=%d) | HOT: R2=%.3f RMSE=%.2f (n=%d)\nDONE\n",cor(D$Tsim,D$Tobs)^2,sqrt(mean((D$Tsim-D$Tobs)^2)),nrow(D),cor(D[hotday==T]$Tsim,D[hotday==T]$Tobs)^2,sqrt(mean((D[hotday==T]$Tsim-D[hotday==T]$Tobs)^2)),nrow(D[hotday==T])))
