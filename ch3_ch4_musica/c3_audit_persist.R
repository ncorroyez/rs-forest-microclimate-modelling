# ==============================================================================
# Chapter 3 — pipeline-audit closure: persist to CSV every skill number that was
# cited in the manuscript but previously only printed by scratchpad scripts.
# Writes authoritative tables so the manuscript is fully reproducible:
#   figP_switches.csv     — LiDAR/S2/switch-LAI/switch-FORMS-H/blend/oracle (recomposed, one method)
#   figQ_krobust.csv      — LiDAR & switch at k=0.5 vs k=0.65
#   figR_optdomain.csv    — "opt" vs ATBD within the opt calibration domain (39 plots)
#   figS_moran.csv        — Moran's I on ΔTmax residual and observed field
#   figT_multiclass.csv   — k-means structural multi-class selection (k=3,4)
#   Rscript c3_audit_persist.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"; CROSS<-3.86; HGATE<-17.8
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn){fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
R2<-function(a,b)suppressWarnings(cor(a,b,use="complete.obs")^2)
df[,st:=ifelse(LAI_ALS<CROSS,"open","dense")]

# ---- figP: switches (one recomposition method) ------------------------------
A<-merge(OBS,getpp("STATIC_ALS"),by="id_plot"); S<-getpp("STATIC_S2_ATBD")
M<-merge(A,S,by="id_plot",suffixes=c(".a",".s")); M<-merge(M,df[,.(id_plot,LAI_ALS,FORMS_H,st)],by="id_plot")
mk<-function(dv,lab){D<-M$st=="dense";O<-M$st=="open";data.table(product=lab,dense=R2(dv[D],M$do[D]),open=R2(dv[O],M$do[O]),
  pooled=R2(dv,M$do),bias=mean(dv-M$do),RMSE=sqrt(mean((dv-M$do)^2)))}
orc<-ifelse(abs(M$ds.a-M$do)<=abs(M$ds.s-M$do),M$ds.a,M$ds.s)
figP<-rbind(mk(M$ds.a,"LiDAR full"),mk(M$ds.s,"S2 ATBD"),
  mk(ifelse(M$LAI_ALS>=CROSS,M$ds.a,M$ds.s),"Switch LAI-keyed"),
  mk(ifelse(M$FORMS_H>=HGATE,M$ds.a,M$ds.s),"Switch FORMS-H"),
  mk(orc,"Oracle per-plot"))
Bl<-merge(OBS,getpp("BLEND_H"),by="id_plot");Bl<-merge(Bl,df[,.(id_plot,st)],by="id_plot")
figP<-rbind(figP,{D<-Bl$st=="dense";O<-Bl$st=="open";data.table(product="Blend FORMS-H",dense=R2(Bl$ds[D],Bl$do[D]),open=R2(Bl$ds[O],Bl$do[O]),pooled=R2(Bl$ds,Bl$do),bias=mean(Bl$ds-Bl$do),RMSE=sqrt(mean((Bl$ds-Bl$do)^2)))})
fwrite(figP,file.path(TAB,"figP_switches.csv"))

# ---- figQ: k robustness -----------------------------------------------------
A5<-merge(OBS,getpp("STATIC_ALS_K05"),by="id_plot");A5<-merge(A5,df[,.(id_plot,st,LAI_ALS)],by="id_plot")
M5<-merge(A5,S,by="id_plot",suffixes=c(".a",".s"))
sw5<-ifelse(M5$LAI_ALS>=CROSS,M5$ds.a,M5$ds.s)
figQ<-data.table(quantity=c("LiDAR k=0.5 dense","LiDAR k=0.5 open","LiDAR k=0.5 pooled","LiDAR k=0.5 bias","Switch k=0.5 dense","Switch k=0.5 open","Switch k=0.5 pooled"),
  value=round(c(R2(A5[st=="dense"]$ds,A5[st=="dense"]$do),R2(A5[st=="open"]$ds,A5[st=="open"]$do),R2(A5$ds,A5$do),mean(A5$ds-A5$do),
    R2(sw5[M5$st=="dense"],M5$do[M5$st=="dense"]),R2(sw5[M5$st=="open"],M5$do[M5$st=="open"]),R2(sw5,M5$do)),3))
fwrite(figQ,file.path(TAB,"figQ_krobust.csv"))

# ---- figR: opt in-domain (LAI_ALS>2 & Hmax>10) ------------------------------
df[,inDom:=LAI_ALS>2 & Hmax>10]
ev<-function(scn){m<-merge(OBS,getpp(scn),by="id_plot");m<-merge(m,df[,.(id_plot,inDom,st)],by="id_plot");d<-m[inDom==TRUE]
  data.table(scenario=scn,dom_n=nrow(d),dom_pooled=R2(d$ds,d$do),dom_dense=R2(d[st=="dense"]$ds,d[st=="dense"]$do),dom_open=R2(d[st=="open"]$ds,d[st=="open"]$do),dom_bias=mean(d$ds-d$do))}
fwrite(rbind(ev("STATIC_S2_OPT"),ev("STATIC_S2_ATBD"),ev("STATIC_ALS"),ev("BLEND_H")),file.path(TAB,"figR_optdomain.csv"))

# ---- figS: Moran's I --------------------------------------------------------
moranI<-function(z,x,y){n<-length(z);D<-as.matrix(dist(cbind(x,y)));W<-1/D;diag(W)<-0;W<-W/rowSums(W);zc<-z-mean(z)
  I<-(n/sum(W))*(sum(W*outer(zc,zc))/sum(zc^2));set.seed(1);perm<-sapply(1:999,function(k){zp<-sample(zc);(n/sum(W))*(sum(W*outer(zp,zp))/sum(zp^2))});c(I=I,p=(1+sum(abs(perm)>=abs(I)))/1000)}
figS<-rbindlist(lapply(c("BLEND_H","FUSION_H","STATIC_ALS"),function(sc){m<-merge(OBS,getpp(sc),by="id_plot");m<-merge(m,df[,.(id_plot,x,y)],by="id_plot");r<-moranI(m$do-m$ds,m$x,m$y);data.table(field=paste0("resid_",sc),I=round(r[1],3),p=round(r[2],3))}))
mo<-merge(OBS,df[,.(id_plot,x,y)],by="id_plot");rf<-moranI(mo$do,mo$x,mo$y);figS<-rbind(figS,data.table(field="observed_dTmax",I=round(rf[1],3),p=round(rf[2],3)))
fwrite(figS,file.path(TAB,"figS_moran.csv"))

# ---- figT: multi-class structural selection (in-sample best per class) ------
X<-scale(as.matrix(df[match(M$id_plot,df$id_plot),.(LAI_ALS,Hmax,fCover,VCI,LCV)]))
figT<-rbindlist(lapply(c(3,4),function(k){set.seed(1);cl<-kmeans(X,centers=k,nstart=25)$cluster
  pick<-sapply(1:k,function(c){ix<-cl==c;if(mean((M$ds.a[ix]-M$do[ix])^2)<=mean((M$ds.s[ix]-M$do[ix])^2))"L" else "S"})
  dv<-ifelse(pick[cl]=="L",M$ds.a,M$ds.s);data.table(k=k,pooled=round(R2(dv,M$do),3),dense=round(R2(dv[M$st=="dense"],M$do[M$st=="dense"]),3),open=round(R2(dv[M$st=="open"],M$do[M$st=="open"]),3))}))
fwrite(figT,file.path(TAB,"figT_multiclass.csv"))
cat("persisted: figP_switches, figQ_krobust, figR_optdomain, figS_moran, figT_multiclass\n")
print(figP,digits=3)
