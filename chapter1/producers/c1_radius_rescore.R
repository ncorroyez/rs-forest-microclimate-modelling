# ==============================================================================
# Footprint-radius sweep, rescored on the aligned clock
# Rescore the footprint-radius sweep against the 53 loggers under the canonical
# time-matched DeltaTmax convention WITH the aligned clock (R/dtmax_convention.R,
# OBS_CLOCK_OFFSET_H = -1). Pure re-extraction; no MuSICA run.
#
# WHY. Section 3.4 quotes r = 0.43 at 5 m ... 0.56 at 50 m from
# tab_radius_by_cluster_convB.csv, which predates the clock alignment (2026-07-28
# vs 2026-07-29) and is referenced by no script, so its clock handling could not
# be established by reading. This settles it by recomputing.
#
# Loop driven by the 53 validation loggers, never by list.files(), so n is
# identical across radii.
# Reads :
#         out_files/radius_test_corr/{5m,10m,12.5m,15m,20m,25m,50m}/musica_out_HOBO_<id>.nc
#         in_files/FR-Blo_2021_v2.nc / in_files/Blois_data_temperature.csv
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv      (stage A2: P labels)
# Out: out_files/Chapter1/tables/tab_radius_sweep_aligned.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1
FORC<-"in_files/FR-Blo_2021_v2.nc"
MREF<-macro_ref(FORC,ds); MOBS<-macro_ref_obs(FORC,ds)
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macH<-macH[as.Date(time)%in%ds][,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]

hob<-as.data.table(read.csv(CFG$hobo_temp_csv))
hob[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hob<-hob[position_sensor=="a" & as.Date(datetime)%in%ds & !id_plot%in%CFG$ids_to_remove]
HO<-hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot,time=floor_date(datetime,"hour"))]
OBS<-HO[,.(obs_dt=delta_tmax_mean(.SD[,.(time,Tmic)],MOBS,ds)),by=id_plot]
IDS<-OBS$id_plot

# LINEAGE. `musica_hobo_radius` is the WIND-UNCORRECTED set (mean canopy wind 2.89
# vs 1.25 in the reported baseline), so it is NOT what Section 3.4 describes.
# `radius_test_corr` carries the corrected baseline (u = 1.233, matching
# musica_hobo_native20/1111 at 1.255) and is the set used here.
RAD<-data.table(radius=c(5,10,12.5,15,20,25,50),
                dir=paste0("out_files/radius_test_corr/",c("5m","10m","12.5m","15m","20m","25m","50m")))
clu<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot,P)]
#' Score one footprint radius against the 53 loggers, pooled and per archetype
#' @param dir directory holding that radius's per-logger NetCDFs
#' @return list of n, pooled r with 95% CI, and one r per archetype
score<-function(dir){
  S<-rbindlist(lapply(IDS,function(id){
    m<-micro_hourly_at(file.path(dir,sprintf("musica_out_HOBO_%s.nc",id)),Z)
    data.table(id_plot=id,sim_dt=if(is.null(m)) NA_real_ else delta_tmax_mean(m,MREF,ds))}))
  D<-merge(OBS,S,by="id_plot")[is.finite(obs_dt)&is.finite(sim_dt)]
  ct<-cor.test(D$obs_dt,D$sim_dt)
  DC<-merge(D,clu,by="id_plot")
  per<-DC[,.(r=if(.N>3) cor(obs_dt,sim_dt) else NA_real_),by=P][order(P)]
  c(list(n=nrow(D),r=unname(ct$estimate),lo=ct$conf.int[1],hi=ct$conf.int[2]),
    setNames(as.list(per$r),per$P))
}
R<-cbind(RAD[,.(radius)],rbindlist(lapply(RAD$dir,function(d) as.data.table(score(d))),fill=TRUE))
fwrite(R,"out_files/Chapter1/tables/tab_radius_sweep_aligned.csv")
cat("\n=== BALAYAGE DU RAYON, horloge alignee, convention B (n identique) ===\n")
print(R[,lapply(.SD,function(x) if(is.numeric(x)) round(x,3) else x)])
cat("\nDONE\n")
