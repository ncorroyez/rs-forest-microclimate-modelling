# ==============================================================================
# c3_slope_junsep_recompute.R
#
# Window-robustness of the Ch3 coupling-slope verdicts (F2/D8 of the harsh review).
# The reframe tables computed the per-plot slope over April-October; the loggers
# are only active June-September and the LiDAR is a summer flight, so the slope is
# recomputed on June-September and compared to April-October with the SAME method.
# Pure re-extraction of existing NetCDFs; no MuSICA.
#
# Scenarios (corrected mapping, verified against reframe_perplot_scenarios.csv):
#   LiDARfixe = STATIC_ALS | S2seul = DYN_S2_ATBD |
#   S2opt = STATIC_S2_OPT (forest-tuned LUT, parametric phenology) |
#   Combinaison = DYN_S2_RESCALED
# Source: out_files/Chapter3/nc_genuine53_windcorr ; macro = musica_in_Blois_pblh.nc
#
# Out: out_files/Chapter3/tables/slope_junsep_vs_aprilct.csv (per plot, both windows)
#      + printed by-archetype summary (mean slope, ranking r, bias) per window.
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("Chapter3_config.R")

FORC<-"in_files/musica_in_Blois_pblh.nc"; Z<-1
BASE<-"out_files/Chapter3/nc_genuine53_windcorr"
MAP<-c(LiDARfixe="STATIC_ALS", S2seul="DYN_S2_ATBD", S2opt="STATIC_S2_OPT", Combinaison="DYN_S2_RESCALED")
WINS<-list(apriloct=c("2021-04-01","2021-10-31"), junsep=c("2021-06-01","2021-09-30"))

# macro (hourly, degC)
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macH<-macH[,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]

# observed micro (HOBO, 1 m proxy = position a), hourly
hob<-as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hob[,time:=floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hob<-hob[position_sensor=="a" & !(id_plot%in%CFG_C3$ids_to_remove)]
HO<-hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
IDS<-sort(unique(HO$id_plot))

slope_win<-function(dt, w){ d<-dt[as.Date(time)>=as.Date(w[1]) & as.Date(time)<=as.Date(w[2])]
  if(nrow(d)<50) return(NA_real_); as.numeric(coef(lm(Tmic~Tmac,d))[2]) }

# per-plot obs slope both windows
OBS<-rbindlist(lapply(IDS,function(id){ m<-merge(HO[id_plot==id],macH,by="time")
  data.table(id_plot=id, obs_apriloct=slope_win(m,WINS$apriloct), obs_junsep=slope_win(m,WINS$junsep)) }))

# per-plot sim slope both windows, each scenario
SIM<-rbindlist(lapply(names(MAP),function(scn){
  rbindlist(lapply(IDS,function(id){
    f<-file.path(BASE,MAP[scn],sprintf("musica_out_HOBO_%s.nc",id))
    m<-if(file.exists(f)) tryCatch(micro_hourly_at(f,Z),error=function(e)NULL) else NULL
    if(is.null(m)) return(data.table(scenario=scn,id_plot=id,sim_apriloct=NA_real_,sim_junsep=NA_real_))
    mm<-merge(m,macH,by="time")
    data.table(scenario=scn,id_plot=id, sim_apriloct=slope_win(mm,WINS$apriloct), sim_junsep=slope_win(mm,WINS$junsep)) }))
}))

# archetype labels from the reframe table
P<-unique(fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/reframe_perplot_scenarios.csv")[,.(id_plot,P)])
D<-merge(merge(SIM,OBS,by="id_plot"),P,by="id_plot",all.x=TRUE)
fwrite(D,"out_files/Chapter3/tables/slope_junsep_vs_aprilct.csv")

# by-archetype summary per window: mean sim slope, ranking r (sim vs obs), bias
summ<-function(sim,obs){ ok<-is.finite(sim)&is.finite(obs)
  list(n=sum(ok), mslope=round(mean(sim[ok]),3),
       r=if(sum(ok)>3) round(cor(sim[ok],obs[ok]),2) else NA,
       bias=round(mean(sim[ok]-obs[ok]),3)) }
for(wn in c("apriloct","junsep")){
  cat(sprintf("\n================ WINDOW: %s ================\n",wn))
  cat(sprintf("%-11s %-3s %4s %8s %6s %7s\n","scenario","P","n","mslope","r","bias"))
  for(scn in names(MAP)) for(p in paste0("P",1:4)){
    d<-D[scenario==scn & P==p]; s<-summ(d[[paste0("sim_",wn)]], d[[paste0("obs_",wn)]])
    if(s$n>0) cat(sprintf("%-11s %-3s %4d %8.3f %6s %7.3f\n",scn,p,s$n,s$mslope,ifelse(is.na(s$r),"NA",s$r),s$bias)) }
}
cat("\nwrote out_files/Chapter3/tables/slope_junsep_vs_aprilct.csv\n")
