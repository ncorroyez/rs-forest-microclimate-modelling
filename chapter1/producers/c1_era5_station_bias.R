# ==============================================================================
# How far is the reanalysis from the station that now forces the model? Scatter plus the
# mean diurnal bias (ERA5 minus station) by hour, with an IQR ribbon.
#
# Reads : in_files/FR-Blo_2021_v2.nc     (station-based forcing: the CHS41 Tair series)
#         in_files/musica_in_Blois.nc    (the older ERA5-forced file)
# Writes: outputs/figures_chap1/B1_era5_station_bias.{png,pdf}
#   Rscript c1_era5_station_bias.R
# ==============================================================================
# Appendix — bias of ERA5 air temperature against the CHS41 station over the study
# period (summer 2021), justifying the station-based forcing. (a) hourly 1:1 scatter;
# (b) mean diurnal bias (ERA5 - station) by hour, IQR ribbon. → FigX_era5_station_bias
suppressPackageStartupMessages({library(ncdf4);library(data.table);library(lubridate);library(ggplot2);library(patchwork);library(ggtext)
  source("scripts/_article_style.R")})
#' Hourly air-temperature series from a MuSICA forcing NetCDF, in degrees C
#' @param p path to the forcing NetCDF
#' @return data.table(time, Ta) on whole hours, NA rows dropped
rd<-function(p){nc<-nc_open(p);tu<-ncatt_get(nc,"time","units")$value
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC");th<-ncvar_get(nc,"time")
  Ta<-as.numeric(ncvar_get(nc,grep("Tair",names(nc$var),value=T)[1]))-273.15
  nc_close(nc);data.table(time=floor_date(t0+th*3600,"hour"),Ta=Ta)[!is.na(Ta)]}
s<-rd("in_files/FR-Blo_2021_v2.nc")[,.(time,sta=Ta)]; e<-rd("in_files/musica_in_Blois.nc")[,.(time,era=Ta)]
D<-merge(s,e,by="time")[as.Date(time)>=as.Date("2021-06-01")&as.Date(time)<=as.Date("2021-09-30")]
bias<-mean(D$era-D$sta); rmse<-sqrt(mean((D$era-D$sta)^2)); r<-cor(D$era,D$sta)
lim<-range(c(D$sta,D$era))
# (a) hourly 1:1 scatter
pa<-ggplot(D,aes(sta,era))+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55")+
  geom_hex(bins=45)+scale_fill_gradient(low="grey85",high="grey15",name="hours")+
  annotate("richtext",x=lim[1],y=lim[2],hjust=0,vjust=1,label.color=NA,fill=NA,size=3.5,
           label=sprintf("bias = %+.2f °C<br>RMSE = %.2f °C<br>*r* = %.2f",bias,rmse,r))+
  coord_equal(xlim=lim,ylim=lim)+
  labs(x="Station air temperature (°C, CHS41)",y="ERA5 air temperature (°C)")+
  theme_article(11)+theme(legend.position=c(0.99,0.02),legend.justification=c(1,0),
    legend.key.height=unit(0.35,"cm"),legend.key.width=unit(0.3,"cm"))
# (b) mean diurnal bias
D[,h:=hour(time)]
di<-D[,.(m=mean(era-sta),lo=quantile(era-sta,.25),hi=quantile(era-sta,.75)),by=h][order(h)]
pb<-ggplot(di,aes(h,m))+geom_hline(yintercept=0,colour="grey55",linewidth=0.4)+
  geom_ribbon(aes(ymin=lo,ymax=hi),fill="#D7191C",alpha=0.15)+
  geom_line(colour="#D7191C",linewidth=0.8)+geom_point(colour="#D7191C",size=1.6)+
  scale_x_continuous(breaks=seq(0,24,6))+
  labs(x="Hour of day (UTC)",y="ERA5 - station air temperature (°C)")+
  theme_article(11)
fig<-(pa|pb)+plot_annotation(tag_levels="a")&theme(plot.tag=element_text(face="bold",size=14))
ggsave_article("outputs/figures_chap1/B1_era5_station_bias",fig,9.5,4.4)   # → Appendix B, Fig. B1
cat(sprintf("bias=%+.2f rmse=%.2f r=%.3f | morning(6-9h) bias=%.2f evening(18h)=%.2f\nDONE\n",
  bias,rmse,r,di[h%in%6:9,mean(m)],di[h==18,m]))
