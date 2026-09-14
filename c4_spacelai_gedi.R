# ==============================================================================
# c4_spacelai_gedi.R — Chapter 4 Part B, the model-free test (FigA axis of Ch3):
# correlation of each per-plot LAI product with the OBSERVED summer microclimate
# (ΔTmax and Tmicro~Tmacro slope, 53 HOBO plots, JJAS), pooled and dense-only.
# Places the GEDI-anchored products next to the Ch3 references:
#   ALS full -0.88/-0.92 | ALS-supervised RF(S2+FORMS-H) -0.68/-0.77 | S2 ~0.
# Writes NC_Full/manuscripts/ch4/tables/Table9_c4_spacelai.csv.
#   Rscript c4_spacelai_gedi.R
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

t7<-fread(file.path(TAB,"Table7_c4_perplot_gedi_lai.csv"));t7[,id_plot:=as.character(id_plot)]
d<-merge(t7,OBS,by="id_plot");dense<-d$LAI_ALS>=3.86
prods<-c(ALS_full="LAI_ALS",S2_raw="LAI_S2_ATBD",GEDI_RF="LAI_GEDI_RF",
         GEDI_RATIO="LAI_GEDI_RATIO",FORMS_H="FORMS_H")
res<-rbindlist(lapply(names(prods),function(p){x<-d[[prods[[p]]]]
  data.table(product=p,
    r_dTmax_pooled=cor(x,d$do), r_slope_pooled=cor(x,d$so,use="complete.obs"),
    r_dTmax_dense=cor(x[dense],d$do[dense]),
    r_slope_dense=cor(x[dense],d$so[dense],use="complete.obs"),
    cor_LAI_ALS_pooled=cor(x,d$LAI_ALS), cor_LAI_ALS_dense=cor(x[dense],d$LAI_ALS[dense]))}))
num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
print(res)
fwrite(res,file.path(TAB,"Table9_c4_spacelai.csv"))
cat("wrote",file.path(TAB,"Table9_c4_spacelai.csv"),"\n")
