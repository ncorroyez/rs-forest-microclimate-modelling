# ==============================================================================
# RE-EXTRACTION of the native-unit sensitivity design on TWO metrics x TWO periods.
#
# WHY. The 3200 NetCDFs of the units design are still on disk
# (out_files/Chapter1/nc_sensitivity_perplot_units, 8 tags x 400 plots), but the
# original extraction (c1_sensitivity_perplot_units.R) kept only ONE number per
# run: DeltaTmax over the full summer. The buffering slope and the hot-day subset
# were therefore missing from the attribution, which forced the Methods to be
# narrowed and the hot-day companion figure to be dropped. This recovers both by
# re-reading the existing outputs. NO MuSICA run.
#
# CONVENTIONS, deliberately identical to the rest of the chapter:
#   DeltaTmax : time-matched, R/dtmax_convention.R, simulations on the forcing clock.
#   slope     : coef(lm(Tmic ~ Tmac))[2] on the hourly series, the SAME definition
#               used by the forward inclusion (scripts/c1_fig5_j1_convB.R:53) and by
#               metrics6 (c1_sensitivity_percluster.R:90), so the attribution and the
#               forward inclusion remain commensurable.
#   hot days  : macro daily maxima >= the 90th percentile, computed exactly as in
#               c1_hot_extract.R, so S2 shares its 13 days with S3/S4/S5.
#
# Writes a NEW table; the original is left untouched until the regression check
# below passes (the recomputed DeltaTmax columns must reproduce it exactly).
# Reads :
#         out_files/Chapter1/nc_sensitivity_perplot_units/S####_<tag>.nc  (400 x 8 = 3200)
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds (design traits)
#         in_files/FR-Blo_2021_v2.nc                                      (macro reference)
# Out: out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(parallel)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"in_files/FR-Blo_2021_v2.nc"; NCDIR<-"out_files/Chapter1/nc_sensitivity_perplot_units"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z_FIX<-1.0
MREF<-macro_ref(FORC,ds)
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
mac<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macH<-mac[as.Date(time)%in%ds][,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]
macD<-macH[,.(Tmx=max(Tmac)),by=.(date=as.Date(time))]
HOT<-macD[Tmx>=quantile(Tmx,0.90),date]
cat(sprintf("hot days: %d\n",length(HOT)))

# --- four numbers per simulation ------------------------------------------------
#' Four scores for one simulation: DeltaTmax and buffering slope, all summer and hot days
#' @param ncp path to a design NetCDF
#' @return named numeric c(dt_all, dt_hot, sl_all, sl_hot); all NA if the run is unreadable
metrics_of<-function(ncp){
  m<-micro_hourly_at(ncp,Z_FIX)
  if(is.null(m)) return(c(dt_all=NA,dt_hot=NA,sl_all=NA,sl_hot=NA))
  mm<-merge(m,macH,by="time"); mh<-mm[as.Date(time)%in%HOT]
  c(dt_all=delta_tmax_mean(m,MREF,ds),
    dt_hot=delta_tmax_mean(m,MREF,HOT,min_days=5),
    sl_all=if(nrow(mm)>50) as.numeric(coef(lm(Tmic~Tmac,mm))[2]) else NA_real_,
    sl_hot=if(nrow(mh)>30) as.numeric(coef(lm(Tmic~Tmac,mh))[2]) else NA_real_)
}
TAGS<-c("base","unifLAD","LAIp","LAIm","Hmaxp","Hmaxm","fCovp","fCovm")
pids<-sort(unique(sub("_[A-Za-z]+\\.nc$","",list.files(NCDIR,"\\.nc$"))))
cat(sprintf("plots: %d x %d tags = %d files\n",length(pids),length(TAGS),length(pids)*length(TAGS)))

jobs<-CJ(pid=pids,tag=TAGS,sorted=FALSE)
res<-mclapply(seq_len(nrow(jobs)),function(i)
  metrics_of(file.path(NCDIR,sprintf("%s_%s.nc",jobs$pid[i],jobs$tag[i]))),mc.cores=8)
M<-cbind(jobs,rbindlist(lapply(res,function(x) as.data.table(as.list(x)))))
W<-dcast(melt(M,id.vars=c("pid","tag")),pid~tag+variable,value.var="value")

# --- steps, from the same sample and the same rules as the original design -------
samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub<-samp[,.SD[sample(.N,min(.N,100))],by=Cluster]; sub[,pid:=sprintf("S%04d",.I)]
D_LAI<-0.5; D_HMAX<-1.0; D_FCOV<-0.10; LAI_FLOOR<-0.1; HMAX_FLOOR<-3; FCOV_MIN<-0.5; FCOV_MAX<-1
S<-sub[,.(pid,Cluster,LAI,Hmax,fCover)]
S[,`:=`(s_LAIp=(LAI+D_LAI)-LAI, s_LAIm=pmax(LAI-D_LAI,LAI_FLOOR)-LAI,
        s_Hmaxp=(Hmax+D_HMAX)-Hmax, s_Hmaxm=pmax(Hmax-D_HMAX,HMAX_FLOOR)-Hmax,
        s_fCovp=pmin(fCover+D_FCOV,FCOV_MAX)-fCover, s_fCovm=pmax(fCover-D_FCOV,FCOV_MIN)-fCover)]
X<-merge(S,W,by="pid")
#' Per-unit sensitivity: (perturbed - base) divided by the step actually taken
#' @param sim metric of the perturbed run
#' @param base metric of the base run
#' @param step the realised trait step (may be truncated by a floor or ceiling)
#' @return numeric per-unit slope, NA where any input is non-finite or the step is ~0
pu<-function(sim,base,step) fifelse(is.finite(sim)&is.finite(base)&is.finite(step)&abs(step)>1e-9,(sim-base)/step,NA_real_)
for(mt in c("dt_all","dt_hot","sl_all","sl_hot")){
  b<-X[[paste0("base_",mt)]]
  for(tr in c("LAI","Hmax","fCov")) for(sd in c("p","m")){
    X[[sprintf("%s_%s_%s",tr,ifelse(sd=="p","up","dn"),mt)]]<-pu(X[[sprintf("%s%s_%s",tr,sd,mt)]],b,X[[paste0("s_",tr,sd)]])
  }
  X[[paste0("LAD_",mt)]]<-b-X[[paste0("unifLAD_",mt)]]     # real minus uniform swap
}
fwrite(X,"out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")

# --- REGRESSION CHECK against the original single-metric table --------------------
O<-rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot_units","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
CMP<-merge(O,X,by="pid",suffixes=c(".old",".new"))
pairs<-list(c("base","base_dt_all"),c("LAI_up","LAI_up_dt_all"),c("LAI_dn","LAI_dn_dt_all"),
            c("Hmax_up","Hmax_up_dt_all"),c("Hmax_dn","Hmax_dn_dt_all"),
            c("fCov_up","fCov_up_dt_all"),c("fCov_dn","fCov_dn_dt_all"),c("dT_LAD","LAD_dt_all"))
cat("\n=== CONTROLE DE NON-REGRESSION (DeltaTmax doit reproduire l'existant) ===\n")
worst<-0
for(p in pairs){
  a<-CMP[[p[1]]]; b<-CMP[[p[2]]]; d<-max(abs(a-b),na.rm=TRUE); worst<-max(worst,d)
  cat(sprintf("  %-10s vs %-18s  ecart max = %.3e   (n=%d)\n",p[1],p[2],d,sum(is.finite(a)&is.finite(b))))
}
cat(sprintf("\n  ECART MAXIMAL GLOBAL = %.3e  -> %s\n",worst,if(worst<1e-9)"IDENTIQUE" else "DIVERGENCE, NE PAS UTILISER"))
cat("DONE\n")
