# ==============================================================================
# c3_junsep_full_recompute_CHS41.R
#
# Same extraction as c3_junsep_full_recompute.R, on the CHS41-Rmerge / NO-wind
# tree (out_files/Chapter3_CHS41/nc, completed by gen_fix_CHS41_dyn.R whose
# outputs live in .../nc_fix and are read as a fallback). This is the convention
# Chapter 1 uses. Two macro references are produced side by side:
#   macro = chs41  -> the forcing that drove the runs, aligned with Chapter 1
#   macro = pblh   -> musica_in_Blois_pblh.nc, what the published extraction used
# June-September only. Pure re-extraction of existing NetCDFs, no simulation.
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source)); source("Chapter3_config_CHS41.R")
Z <- 1; WIN <- c("2021-06-01","2021-09-30")
BASE <- file.path(CFG_C3$out_dir,"nc"); FIX <- file.path(CFG_C3$out_dir,"nc_fix")
MAP  <- c(LiDARfixe="STATIC_ALS", S2seul="DYN_S2_ATBD", S2opt="STATIC_S2_OPT", Combinaison="DYN_S2_RESCALED")
MACROS <- list(chs41=CFG_C3$forcing_file, pblh="in_files/musica_in_Blois_pblh.nc")
OUTDIR <- file.path(CFG_C3$out_dir,"tables"); dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
MINSZ <- 1e6

read_macro <- function(p){ nc<-nc_open(p); tu<-ncatt_get(nc,"time","units")$value; th<-ncvar_get(nc,"time")
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
  d<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15); nc_close(nc)
  d[,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time] }

hob <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hob[, time := floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hob <- hob[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
HO  <- hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
IDS <- sort(unique(HO$id_plot)); stopifnot(length(IDS)==53)

mets <- function(dt){ d<-dt[as.Date(time)>=as.Date(WIN[1]) & as.Date(time)<=as.Date(WIN[2])]
  if(nrow(d)<50) return(list(slope=NA_real_,dtmax=NA_real_,dtmin=NA_real_))
  sl<-as.numeric(coef(lm(Tmic~Tmac,d))[2])
  dd<-d[,.(mxi=max(Tmic,na.rm=TRUE),mxa=max(Tmac,na.rm=TRUE),
           mni=min(Tmic,na.rm=TRUE),mna=min(Tmac,na.rm=TRUE)),by=.(date=as.Date(time))]
  list(slope=sl,dtmax=mean(dd$mxi-dd$mxa,na.rm=TRUE),dtmin=mean(dd$mni-dd$mna,na.rm=TRUE)) }

pick <- function(scn,id){ a<-file.path(BASE,MAP[scn],sprintf("musica_out_HOBO_%s.nc",id))
  b<-file.path(FIX ,MAP[scn],sprintf("musica_out_HOBO_%s.nc",id))
  if(file.exists(a)&&file.size(a)>=MINSZ) a else if(file.exists(b)&&file.size(b)>=MINSZ) b else NA_character_ }

P <- unique(fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/reframe_perplot_scenarios.csv")[,.(id_plot,P)])

for (mk in names(MACROS)) {
  macH <- read_macro(MACROS[[mk]])
  OBS <- rbindlist(lapply(IDS,function(id){ r<-mets(merge(HO[id_plot==id],macH,by="time"))
    data.table(id_plot=id,obs_slope_junsep=r$slope,obs_dtmax_junsep=r$dtmax,obs_dtmin_junsep=r$dtmin) }))
  SIM <- rbindlist(lapply(names(MAP),function(scn){ rbindlist(lapply(IDS,function(id){
    f<-pick(scn,id)
    if(is.na(f)) return(data.table(scenario=scn,id_plot=id,sim_slope_junsep=NA_real_,sim_dtmax_junsep=NA_real_,sim_dtmin_junsep=NA_real_))
    m<-tryCatch(micro_hourly_at(f,Z),error=function(e)NULL)
    if(is.null(m)) return(data.table(scenario=scn,id_plot=id,sim_slope_junsep=NA_real_,sim_dtmax_junsep=NA_real_,sim_dtmin_junsep=NA_real_))
    r<-mets(merge(m,macH,by="time"))
    data.table(scenario=scn,id_plot=id,sim_slope_junsep=r$slope,sim_dtmax_junsep=r$dtmax,sim_dtmin_junsep=r$dtmin) })) }))
  D <- merge(merge(SIM,OBS,by="id_plot"),P,by="id_plot",all.x=TRUE)
  stopifnot(nrow(D)==53*4, sum(is.na(D$sim_dtmax_junsep))==0)
  fwrite(D, file.path(OUTDIR,sprintf("junsep_chs41_perplot_macro-%s.csv",mk)))

  LONG <- rbindlist(lapply(c("slope","dtmax","dtmin"), function(met)
    D[,.(scenario,id_plot,P,metric=met,sim=get(sprintf("sim_%s_junsep",met)),obs=get(sprintf("obs_%s_junsep",met)))]))
  SUMM <- LONG[is.finite(sim)&is.finite(obs),
    .(n=.N, mean_sim=round(mean(sim),3), mean_obs=round(mean(obs),3),
      r=if(.N>3) round(cor(sim,obs),2) else NA_real_,
      bias=round(mean(sim-obs),3), sd_res=round(sd(sim-obs),3),
      disp=round(sd(sim)/sd(obs),3)), by=.(metric,scenario,P)][order(metric,scenario,P)]
  ALLP <- LONG[is.finite(sim)&is.finite(obs),
    .(P="ALL", n=.N, mean_sim=round(mean(sim),3), mean_obs=round(mean(obs),3),
      r=round(cor(sim,obs),2), bias=round(mean(sim-obs),3), sd_res=round(sd(sim-obs),3),
      disp=round(sd(sim)/sd(obs),3)), by=.(metric,scenario)]
  SUMM <- rbind(SUMM, ALLP, use.names=TRUE)[order(metric,scenario,P)]
  fwrite(SUMM, file.path(OUTDIR,sprintf("junsep_chs41_by_archetype_macro-%s.csv",mk)))

  set.seed(1); B <- 2000
  paired <- function(x){ x<-x[is.finite(x)]; if(length(x)<3) return(c(NA,NA,NA))
    bs<-replicate(B,mean(sample(x,length(x),TRUE))); c(mean(x),unname(quantile(bs,c(.025,.975)))) }
  CON <- rbindlist(lapply(c("slope","dtmax"), function(met){
    W <- dcast(LONG[metric==met], id_plot+P~scenario, value.var="sim")
    rbindlist(lapply(setdiff(names(MAP),"LiDARfixe"), function(scn){
      W[,{ s<-paired(get(scn)-LiDARfixe)
        .(metric=met, contrast=paste0(scn,"-LiDARfixe"), n=sum(is.finite(get(scn)-LiDARfixe)),
          mean=round(s[1],3), lo=round(s[2],3), hi=round(s[3],3)) }, by=P] })) }))
  fwrite(CON, file.path(OUTDIR,sprintf("junsep_chs41_paired_contrasts_macro-%s.csv",mk)))
  cat(sprintf("macro=%s written (n=%d rows)\n", mk, nrow(D)))
}
cat("done\n")
