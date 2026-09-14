# ==============================================================================
# Chapter 3 — extract the k = 0.5 vs k = 0.65 counterfactual (Appendix).
# Same extraction as c3_junsep_full_recompute_CHS41.R: CHS41 macro, June-Sept,
# z = 1 m. Pure re-extraction, no simulation. Answers the question Chapter 1
# refers here: does a uniform 23 % leaf-area reduction change the simulated
# coupling, and does the between-plot ranking survive it?
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source)); source("Chapter3_config_CHS41.R")
Z <- 1; WIN <- c("2021-06-01","2021-09-30"); MINSZ <- 1e6
TREES <- c(k050 = file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
           k065 = file.path(CFG_C3$out_dir,"nc_k065","STATIC_ALS_K065"))
OUTDIR <- file.path(CFG_C3$out_dir,"tables")

read_macro <- function(p){ nc<-nc_open(p); tu<-ncatt_get(nc,"time","units")$value; th<-ncvar_get(nc,"time")
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
  d<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15); nc_close(nc)
  d[,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time] }
mets <- function(dt){ d<-dt[as.Date(time)>=as.Date(WIN[1]) & as.Date(time)<=as.Date(WIN[2])]
  if(nrow(d)<50) return(list(slope=NA_real_,dtmax=NA_real_))
  sl<-as.numeric(coef(lm(Tmic~Tmac,d))[2])
  dd<-d[,.(mxi=max(Tmic,na.rm=TRUE),mxa=max(Tmac,na.rm=TRUE)),by=.(date=as.Date(time))]
  list(slope=sl,dtmax=mean(dd$mxi-dd$mxa,na.rm=TRUE)) }

hob <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hob[, time := floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hob <- hob[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
HO  <- hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
IDS <- sort(unique(HO$id_plot)); stopifnot(length(IDS)==53)
macH <- read_macro(CFG_C3$forcing_file)

OBS <- rbindlist(lapply(IDS,function(id){ r<-mets(merge(HO[id_plot==id],macH,by="time"))
  data.table(id_plot=id,obs_slope=r$slope,obs_dtmax=r$dtmax) }))
SIM <- rbindlist(lapply(names(TREES),function(k){ rbindlist(lapply(IDS,function(id){
  f<-file.path(TREES[[k]],sprintf("musica_out_HOBO_%s.nc",id)); stopifnot(file.exists(f),file.size(f)>=MINSZ)
  r<-mets(merge(micro_hourly_at(f,Z),macH,by="time"))
  data.table(k=k,id_plot=id,sim_slope=r$slope,sim_dtmax=r$dtmax) })) }))
P <- unique(fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/reframe_perplot_scenarios.csv")[,.(id_plot,P)])
D <- merge(merge(SIM,OBS,by="id_plot"),P,by="id_plot",all.x=TRUE)
stopifnot(nrow(D)==53*2, !any(is.na(D$sim_slope)))
fwrite(D, file.path(OUTDIR,"k065_perplot_macro-chs41.csv"))

W <- dcast(D, id_plot+P+obs_slope+obs_dtmax ~ k, value.var=c("sim_slope","sim_dtmax"))
W[, d_slope := sim_slope_k065 - sim_slope_k050][, d_dtmax := sim_dtmax_k065 - sim_dtmax_k050]
fwrite(W, file.path(OUTDIR,"k065_paired_macro-chs41.csv"))

cat("\n=== Does the 23 % reduction change the simulated coupling? (n = 53, June-Sept, CHS41) ===\n")
cat(sprintf("slope  : mean shift %+.4f  [P1 %+.4f | P2 %+.4f | P3 %+.4f | P4 %+.4f]\n",
  mean(W$d_slope), mean(W[P=="P1"]$d_slope), mean(W[P=="P2"]$d_slope),
  mean(W[P=="P3"]$d_slope), mean(W[P=="P4"]$d_slope)))
cat(sprintf("dTmax  : mean shift %+.3f C [P1 %+.3f | P2 %+.3f | P3 %+.3f | P4 %+.3f]\n",
  mean(W$d_dtmax), mean(W[P=="P1"]$d_dtmax), mean(W[P=="P2"]$d_dtmax),
  mean(W[P=="P3"]$d_dtmax), mean(W[P=="P4"]$d_dtmax)))
set.seed(1); B <- 2000
bo <- function(x){ q<-quantile(replicate(B,mean(sample(x,length(x),TRUE))),c(.025,.975)); sprintf("[%+.4f, %+.4f]",q[1],q[2]) }
cat(sprintf("slope shift 95%% CI %s   dTmax shift 95%% CI %s\n", bo(W$d_slope), bo(W$d_dtmax)))

cat("\n=== Does the between-plot RANKING survive? r(sim, obs) ===\n")
for (m in c("slope","dtmax")) for (k in c("k050","k065"))
  cat(sprintf("  %-5s %s : r = %.3f\n", m, k,
      cor(W[[paste0("sim_",m,"_",k)]], W[[paste0("obs_",m)]], use="complete.obs")))
cat(sprintf("\n  r(sim_slope_k050, sim_slope_k065) = %.4f   (rank rho = %.4f)\n",
    cor(W$sim_slope_k050,W$sim_slope_k065), cor(W$sim_slope_k050,W$sim_slope_k065,method="spearman")))
cat(sprintf("  r(sim_dtmax_k050, sim_dtmax_k065) = %.4f   (rank rho = %.4f)\n",
    cor(W$sim_dtmax_k050,W$sim_dtmax_k065), cor(W$sim_dtmax_k050,W$sim_dtmax_k065,method="spearman")))
cat("\n=== Warm bias (sim - obs), by k ===\n")
for (k in c("k050","k065")) cat(sprintf("  %s : dTmax bias %+.3f C\n", k, mean(W[[paste0("sim_dtmax_",k)]]-W$obs_dtmax)))
cat("\nDONE\n")
