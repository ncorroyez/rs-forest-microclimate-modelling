# ==============================================================================
# c4_disc_seeds2.R — three more observation-only seeds (2021 only, per Nathan) for the General
# Discussion (no simulations):
#  B. NIGHT: does ΔTmin couple to LAI like ΔTmax does? (limit of the LAI
#     paradigm; summer 2021, macro from the forcing file)
#  C. RECIPROCAL Wu test: winter microclimate (Jan 1 - Mar 14, 2021, leafless)
#     vs per-plot S2 green-up date, LAI-partialled.
#  D. HOW MANY LOGGERS: bootstrap subsampling of the slope~LAI relation
#     (n = 10..53) -> CI width vs network size.
# Outputs: NC_Full/manuscripts/ch4/tables/Table24*.
#   Rscript c4_disc_seeds2.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df[,id_plot:=as.character(id_plot)];df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]
hb[,`:=`(id_plot=as.character(id_plot),date=as.Date(time))]

# ---- B. night vs day coupling (summer 2021, macro from forcing) --------------
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MA<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds]
mm<-MA[,.(Tmx=max(Tm,na.rm=TRUE),Tmn=min(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
hd<-hb[date%in%ds,.(Tmx_o=max(t_hobo,na.rm=TRUE),Tmn_o=min(t_hobo,na.rm=TRUE)),by=.(id_plot,date)]
hd<-merge(hd,mm,by="date")
nb<-hd[,.(dmax=mean(Tmx_o-Tmx,na.rm=TRUE),dmin=mean(Tmn_o-Tmn,na.rm=TRUE)),by=id_plot]
nb<-merge(nb,df[,.(id_plot,LAI_ALS,Hmax)],by="id_plot")
cat("\n=== B. day vs night coupling to LAI (summer 2021, n =",nrow(nb),") ===\n")
cat(sprintf("cor(dTmax, LAI) = %.3f | cor(dTmin, LAI) = %.3f | cor(dTmin, Hmax) = %.3f\n",
  nb[,cor(dmax,LAI_ALS)],nb[,cor(dmin,LAI_ALS)],nb[,cor(dmin,Hmax)]))
cat(sprintf("mean dTmax = %.2f | mean dTmin = %+.2f | sd dTmin = %.2f\n",
  nb[,mean(dmax)],nb[,mean(dmin)],nb[,sd(dmin)]))
fwrite(nb,file.path(TAB,"Table24b_c4_day_night.csv"))

# ---- C. reciprocal Wu test: winter T -> green-up ------------------------------
tb<-readRDS(file.path(CFG_C3$out_dir,"lai_prep","ts_by_plot.rds"))
greenup<-rbindlist(lapply(names(tb$atbd),function(p){a<-tb$atbd[[p]]
  if(is.null(a)||!nrow(a))return(NULL)
  half<-min(a$lai)+diff(range(a$lai))/2
  data.table(pid=p,sos_s2=a$doy[which(a$lai>=half)[1]])}))
wint<-seq(as.Date("2021-01-01"),as.Date("2021-03-14"),by="day")
tw<-hb[date%in%wint,.(Tmx_o=max(t_hobo,na.rm=TRUE)),by=.(id_plot,date)][,.(t_winter=mean(Tmx_o)),by=id_plot]
C<-merge(merge(greenup,df[,.(pid,id_plot,LAI_ALS)],by="pid"),tw,by="id_plot")
cat("\n=== C. winter microclimate -> green-up (n =",nrow(C),") ===\n")
cat(sprintf("cor(SOS, winter Tmax) = %.3f | partial | LAI = %.3f\n",
  C[,cor(sos_s2,t_winter)],
  cor(residuals(lm(sos_s2~LAI_ALS,C)),residuals(lm(t_winter~LAI_ALS,C)))))
fwrite(C,file.path(TAB,"Table24c_c4_winter_sos.csv"))

# ---- D. how many loggers? bootstrap the slope~LAI relation --------------------
slp<-function(y,x) as.numeric(coef(lm(y~x))[2])
MAh<-MA[,hr:=floor_date(time,"hour")]
hbh<-merge(hb[date%in%ds,.(id_plot,Tobs=t_hobo,hr=floor_date(time,"hour"))],MAh[,.(hr,Tm)],by="hr")
SO<-hbh[,.(so=slp(Tobs,Tm)),by=id_plot]
D<-merge(SO,df[,.(id_plot,LAI_ALS)],by="id_plot")
set.seed(42);B<-2000L
res_n<-rbindlist(lapply(c(10,15,20,30,40,53),function(nn){
  bs<-replicate(B,{i<-sample(nrow(D),nn,replace=FALSE)
    c(cor(D$so[i],D$LAI_ALS[i]),slp(D$so[i],D$LAI_ALS[i]))})
  data.table(n=nn,r_med=median(bs[1,]),r_lo=quantile(bs[1,],.025),r_hi=quantile(bs[1,],.975),
             b_med=median(bs[2,]),b_lo=quantile(bs[2,],.025),b_hi=quantile(bs[2,],.975))}))
num<-names(res_n)[sapply(res_n,is.numeric)];res_n[,(num):=lapply(.SD,round,3),.SDcols=num]
cat("\n=== D. subsampling the slope~LAI relation (95% CI across",B,"draws) ===\n")
print(res_n)
fwrite(res_n,file.path(TAB,"Table24d_c4_how_many_loggers.csv"))
cat("\nDONE\n")
