# ==============================================================================
# c4_disc_seeds3_sos_fix.R — recompute per-plot S2 SOS with the TRS50 method of
# c4_phenofit_metrics.R (weekly interpolation, phenofit::whit2 lambda = 10,
# half-range crossing) and, crucially, take the LAST upward crossing before the
# seasonal peak (walk back from the peak), so winter noise crossings (DOY 26-56 in
# Table23c) are ignored. Re-runs the vernal (spring Tmax, Mar 15 - Apr 30) and
# winter (Jan 1 - Mar 14) tests, LAI-partialled. 2021 only.
# Outputs: NC_Full/manuscripts/ch4/tables/Table23c_c4_pheno_loop_v2_",SER,".csv,
#          Table24c_c4_winter_sos_v2_",SER,".csv, Table25_c4_sos_fix_summary_",SER,".csv
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(phenofit)})
source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df[,id_plot:=as.character(id_plot)];df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]
hb[,`:=`(id_plot=as.character(id_plot),date=as.Date(time))]
tb<-readRDS(file.path(CFG_C3$out_dir,"lai_prep","ts_by_plot.rds"))

sos_trs50<-function(doy,val){
  ok<-is.finite(doy)&is.finite(val); if(sum(ok)<6) return(c(NA,NA,NA))
  wk<-seq(min(doy[ok]),max(doy[ok]),by=7); vw<-approx(doy[ok],val[ok],wk)$y
  sm<-as.numeric(whit2(vw,lambda=10))
  half<-min(sm)+(max(sm)-min(sm))/2; ipk<-which.max(sm)
  i<-ipk; while(i>1 && sm[i-1]>=half) i<-i-1          # walk back from the peak
  raw<-wk[which(sm>=half)[1]]                          # old-style first crossing
  c(sos=wk[i],sos_first=raw,peak=wk[ipk])}
SER<-Sys.getenv("SOS_SERIES","annual"); cat("series used:",SER,"\n")
gu<-rbindlist(lapply(names(tb[[SER]]),function(p){a<-tb[[SER]][[p]]
  if(is.null(a)||!nrow(a))return(NULL); s<-sos_trs50(a$doy,a$lai)
  data.table(pid=p,sos_s2=s[1],sos_first=s[2],doy_peak=s[3],n_dates=nrow(a))}))
gu<-gu[is.finite(sos_s2)]
cat("SOS v2 range:",range(gu$sos_s2),"| n =",nrow(gu),
    "| n plots where old first-crossing < DOY 60:",sum(gu$sos_first<60,na.rm=TRUE),"\n")

pcor<-function(x,y,z) cor(residuals(lm(x~z)),residuals(lm(y~z)))
mk<-function(win,lab,tag){
  tw<-hb[date%in%win,.(Tmx_o=max(t_hobo,na.rm=TRUE)),by=.(id_plot,date)][,.(t=mean(Tmx_o)),by=id_plot]
  M<-merge(merge(gu,df[,.(pid,id_plot,LAI_ALS,Hmax)],by="pid"),tw,by="id_plot")
  setnames(M,"t",paste0("t_",tag))
  out<-data.table(window=lab,n=nrow(M),
    r_sos_T=round(cor(M$sos_s2,M[[paste0("t_",tag)]]),3),
    r_partial_LAI=round(pcor(M$sos_s2,M[[paste0("t_",tag)]],M$LAI_ALS),3),
    r_partial_LAI_Hmax=round(cor(residuals(lm(sos_s2~LAI_ALS+Hmax,M)),
                                 residuals(lm(M[[paste0("t_",tag)]]~LAI_ALS+Hmax,M))),3),
    r_sos_LAI=round(cor(M$sos_s2,M$LAI_ALS),3),
    r_T_LAI=round(cor(M[[paste0("t_",tag)]],M$LAI_ALS),3),
    r_OLDsos_T=round(cor(M$sos_first,M[[paste0("t_",tag)]]),3))
  list(M=M,out=out)}
sp<-mk(seq(as.Date("2021-03-15"),as.Date("2021-04-30"),by="day"),"spring Mar15-Apr30 2021","spring")
wi<-mk(seq(as.Date("2021-01-01"),as.Date("2021-03-14"),by="day"),"winter Jan1-Mar14 2021","winter")
res<-rbind(sp$out,wi$out); print(res)
fwrite(sp$M,file.path(TAB,paste0("Table23c_c4_pheno_loop_v2_",SER,".csv")))
fwrite(wi$M,file.path(TAB,paste0("Table24c_c4_winter_sos_v2_",SER,".csv")))
fwrite(res,file.path(TAB,paste0("Table25_c4_sos_fix_summary_",SER,".csv")))
cat("\nDONE\n")
