# ==============================================================================
# RESCORE of the correction-level variants against the 53 loggers, under the
# canonical time-matched DeltaTmax convention with the aligned clock
# (R/dtmax_convention.R, OBS_CLOCK_OFFSET_H = -1).
#
# WHY. Section 4.2 and Appendix F quote a sequence over the leaf-area/wind
# correction levels, but the stale numbers were produced under the old
# convention (independent daily maxima, -2 h shift). Every variant already
# exists on disk, so this is a pure re-extraction; no MuSICA run.
#
# GATE (verified before writing this script): all six variant sets share the
# forcing clock ("hours since 2021-01-01 00:00:00", 8759 steps), so a single
# MREF/MOBS pair scores them all. The loop is driven by the 53 validation
# loggers, never by list.files(), so n is identical across variants.
#
# Reads :
#         four correction levels, each a directory of per-logger NetCDFs:
#         out_files/hobo_native20_uncorr/musica_out_HOBO_<id>.nc      (no correction)
#         out_files/hobo_native20_windonly/musica_out_HOBO_<id>.nc    (wind only)
#         out_files/musica_hobo_native20/1111/musica_out_HOBO_<id>.nc (wind + scan angle: BASELINE)
#         out_files/hobo_native20_k065/musica_out_HOBO_<id>.nc        (+ k = 0.65)
#         in_files/FR-Blo_2021_v2.nc / in_files/Blois_data_temperature.csv
# Out: out_files/Chapter1/tables/tab_correction_levels.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1
FORC<-"in_files/FR-Blo_2021_v2.nc"
MREF<-macro_ref(FORC,ds)          # forcing clock  -> SIMULATIONS
MOBS<-macro_ref_obs(FORC,ds)      # logger clock   -> OBSERVATIONS
macH<-{nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
  d<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
  d[as.Date(time)%in%ds][,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]}

# ---- OBSERVATIONS (identical recipe to the headline validation) ---------------
hob<-as.data.table(read.csv(CFG$hobo_temp_csv))
hob[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hob<-hob[position_sensor=="a" & as.Date(datetime)%in%ds & !id_plot%in%CFG$ids_to_remove]
hob[,time:=floor_date(datetime,"hour")]
HO<-hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot,time)]
OBS<-HO[,{m<-.SD[,.(time,Tmic)]
          mm<-merge(m,macH,by="time")
          .(obs_dt=delta_tmax_mean(m,MOBS,ds),
            obs_sl=if(nrow(mm)>50) as.numeric(coef(lm(Tmic~Tmac,mm))[2]) else NA_real_)},by=id_plot]
IDS<-OBS$id_plot                                   # the 53, fixed across variants

# ---- variant ledger -----------------------------------------------------------
# Two axes are deliberately kept apart. Wind: uncorr -> windonly -> full.
# Leaf-area retrieval: full (k = 0.5 + scan angle) -> k065.
VAR<-data.table(
  variant=c("uncorr","windonly","full","k065"),
  dir  =c("out_files/hobo_native20_uncorr","out_files/hobo_native20_windonly",
          "out_files/musica_hobo_native20/1111","out_files/hobo_native20_k065"),
  axis =c("wind","wind","wind+leaf area","leaf area"),
  label=c("no correction","wind only","wind + scan angle (baseline)","+ k = 0.65"))

PER<-list()
#' Score one correction variant against the 53 loggers (DeltaTmax and slope)
#' @param dir directory of that variant's per-logger NetCDFs
#' @return one-row data.table of n, r + CI, bias, RMSE, amplitude and slope agreement.
#'   Side effect: appends the per-plot values to PER, which feeds Figure F2 (stage B16).
score<-function(dir){
  S<-rbindlist(lapply(IDS,function(id){
    m<-micro_hourly_at(file.path(dir,sprintf("musica_out_HOBO_%s.nc",id)),Z)
    if(is.null(m)) return(data.table(id_plot=id,sim_dt=NA_real_,sim_sl=NA_real_))
    mm<-merge(m,macH,by="time")
    data.table(id_plot=id, sim_dt=delta_tmax_mean(m,MREF,ds),
               sim_sl=if(nrow(mm)>50) as.numeric(coef(lm(Tmic~Tmac,mm))[2]) else NA_real_)}))
  D<-merge(OBS,S,by="id_plot")
  PER[[length(PER)+1L]]<<-copy(D)[,dir:=dir]      # keep per-plot values for Fig. F2
  ok<-is.finite(D$obs_dt)&is.finite(D$sim_dt); o<-D$obs_dt[ok]; s<-D$sim_dt[ok]
  ok2<-is.finite(D$obs_sl)&is.finite(D$sim_sl)
  ct<-cor.test(o,s)
  data.table(n=sum(ok), r=ct$estimate, r_lo=ct$conf.int[1], r_hi=ct$conf.int[2],
             bias=mean(s-o), rmse=sqrt(mean((s-o)^2)),
             amp=100*coef(lm(s~o))[2], sim_sd=sd(s), obs_sd=sd(o),
             slope_r=cor(D$obs_sl[ok2],D$sim_sl[ok2]), slope_sim=mean(D$sim_sl[ok2]))
}
R<-cbind(VAR,rbindlist(lapply(VAR$dir,score)))
PP<-merge(rbindlist(PER),VAR[,.(dir,variant,label)],by="dir")
fwrite(PP,"out_files/Chapter1/tables/tab_correction_levels_perplot.csv")
fwrite(R,"out_files/Chapter1/tables/tab_correction_levels.csv")

cat("\n=== NIVEAUX DE CORRECTION, horloge alignee, n identique (convention B) ===\n")
cat(sprintf("%-30s %-14s %6s %16s %7s %6s %8s\n","variante","axe","r","IC95","biais","ampl","pente r"))
for(i in seq_len(nrow(R))) cat(sprintf("%-30s %-14s %6.3f  [%5.2f,%5.2f] %+6.2f %5.0f%% %8.3f\n",
  R$label[i],R$axis[i],R$r[i],R$r_lo[i],R$r_hi[i],R$bias[i],R$amp[i],R$slope_r[i]))
cat("DONE\n")

# ---- density split: is the correction effect concentrated in open stands? -----
clu<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot,LAI)]
md<-median(merge(OBS,clu,by="id_plot")$LAI,na.rm=TRUE)
#' Same variant, split at the median LAI: is the correction concentrated in open stands?
#' @param dir directory of that variant's per-logger NetCDFs
#' @return one-row data.table of bias and r in the dense and in the sparse half
split_score<-function(dir){
  S<-rbindlist(lapply(IDS,function(id){
    m<-micro_hourly_at(file.path(dir,sprintf("musica_out_HOBO_%s.nc",id)),Z)
    if(is.null(m)) return(data.table(id_plot=id,sim_dt=NA_real_))
    data.table(id_plot=id,sim_dt=delta_tmax_mean(m,MREF,ds))}))
  D<-merge(merge(OBS,S,by="id_plot"),clu,by="id_plot")[is.finite(obs_dt)&is.finite(sim_dt)&is.finite(LAI)]
  d<-D[LAI>=md]; s<-D[LAI<md]
  data.table(bias_dense=mean(d$sim_dt-d$obs_dt), r_dense=cor(d$obs_dt,d$sim_dt),
             bias_sparse=mean(s$sim_dt-s$obs_dt), r_sparse=cor(s$obs_dt,s$sim_dt))
}
SP<-cbind(VAR[,.(label)],rbindlist(lapply(VAR$dir,split_score)))
fwrite(SP,"out_files/Chapter1/tables/tab_correction_levels_split.csv")
cat("\n=== split densite (mediane LAI) ===\n")
cat(sprintf("%-30s %12s %8s %13s %9s\n","variante","biais dense","r dense","biais sparse","r sparse"))
for(i in seq_len(nrow(SP))) cat(sprintf("%-30s %+11.2f %8.3f %+12.2f %9.3f\n",
  SP$label[i],SP$bias_dense[i],SP$r_dense[i],SP$bias_sparse[i],SP$r_sparse[i]))
cat("DONE2\n")
