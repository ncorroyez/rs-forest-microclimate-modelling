# ==============================================================================
# DIAGNOSTIC ONLY — re-extract the Ch3 June-September validation from the
# CHS41-Rmerge / no-wind NetCDF tree (out_files/Chapter3_CHS41/nc), and compare it
# with the published extraction, which came from the wind-corrected July tree
# (out_files/Chapter3/nc_genuine53_windcorr). Also extracts the night-time offset
# dTmin, which exists nowhere in the Ch3 outputs. Writes nothing into any
# manuscript. Pure re-extraction of existing NetCDFs, no simulation.
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
setwd(path.expand("~/Documents/z_Example_rmusica_31012025"))
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source)); source("Chapter3_config.R")
Z <- 1
WIN <- c("2021-06-01","2021-09-30")
TREES <- list(windcorr = "out_files/Chapter3/nc_genuine53_windcorr",
              chs41    = "out_files/Chapter3_CHS41/nc")
MAP <- c(LiDARfixe="STATIC_ALS", S2seul="DYN_S2_ATBD", S2opt="STATIC_S2_OPT",
         Combinaison="DYN_S2_RESCALED")
MACROS <- list(pblh  = "in_files/musica_in_Blois_pblh.nc",
               chs41 = "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc")

read_macro <- function(p){
  nc <- nc_open(p); tu <- ncatt_get(nc,"time","units")$value; th <- ncvar_get(nc,"time")
  t0 <- as.POSIXct(sub(".*since ","",tu),tz="UTC")
  d <- data.table(time=floor_date(t0+th*3600,"hour"),
                  Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15); nc_close(nc)
  d[,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]
}

hob <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hob[, time := floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hob <- hob[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
HO  <- hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
IDS <- sort(unique(HO$id_plot)); cat("n plots:",length(IDS),"\n")

mets <- function(dt){
  d <- dt[as.Date(time)>=as.Date(WIN[1]) & as.Date(time)<=as.Date(WIN[2])]
  if (nrow(d) < 50) return(list(slope=NA_real_,dtmax=NA_real_,dtmin=NA_real_))
  sl <- as.numeric(coef(lm(Tmic~Tmac,d))[2])
  dd <- d[,.(mxi=max(Tmic,na.rm=TRUE),mxa=max(Tmac,na.rm=TRUE),
             mni=min(Tmic,na.rm=TRUE),mna=min(Tmac,na.rm=TRUE)),by=.(date=as.Date(time))]
  list(slope=sl, dtmax=mean(dd$mxi-dd$mxa,na.rm=TRUE), dtmin=mean(dd$mni-dd$mna,na.rm=TRUE))
}

OUT <- rbindlist(lapply(names(MACROS), function(mk){
  macH <- read_macro(MACROS[[mk]])
  obs <- rbindlist(lapply(IDS,function(id){
    m <- merge(HO[id_plot==id],macH,by="time"); r <- mets(m)
    data.table(id_plot=id, obs_slope=r$slope, obs_dtmax=r$dtmax, obs_dtmin=r$dtmin)}))
  rbindlist(lapply(names(TREES), function(tk){
    rbindlist(lapply(names(MAP), function(scn){
      sim <- rbindlist(lapply(IDS,function(id){
        f <- file.path(TREES[[tk]],MAP[scn],sprintf("musica_out_HOBO_%s.nc",id))
        if(!file.exists(f)) return(data.table(id_plot=id,sim_slope=NA_real_,sim_dtmax=NA_real_,sim_dtmin=NA_real_))
        m <- tryCatch(micro_hourly_at(f,Z),error=function(e)NULL)
        if(is.null(m)) return(data.table(id_plot=id,sim_slope=NA_real_,sim_dtmax=NA_real_,sim_dtmin=NA_real_))
        r <- mets(merge(m,macH,by="time"))
        data.table(id_plot=id,sim_slope=r$slope,sim_dtmax=r$dtmax,sim_dtmin=r$dtmin)}))
      cbind(data.table(macro=mk,tree=tk,scenario=scn), merge(obs,sim,by="id_plot"))
    }))}))}))
fwrite(OUT, "/tmp/claude-1001/-home-corroyez-Documents-NC-Full/784edf00-1161-462e-8f45-19b8ca8dbce1/scratchpad/c3_chs41_junsep_diag_perplot.csv")
cat("rows:",nrow(OUT),"\n")
