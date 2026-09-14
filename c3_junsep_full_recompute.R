# ==============================================================================
# c3_junsep_full_recompute.R
#
# Full re-extraction of the Ch3 validation on BOTH windows (April-October, the old
# reframe window, and June-September, the logger/LiDAR summer window), same method,
# so the manuscript can be moved to June-September consistently. Metrics per plot,
# per scenario: coupling slope and ΔTmax, simulated and observed. Aggregated by
# archetype (ranking r, mean, bias, residual dispersion) and as paired contrasts
# against LiDAR fixe (bootstrap CI). Pure re-extraction of existing NetCDFs.
#
# Scenarios (verified mapping): LiDARfixe=STATIC_ALS, S2seul=DYN_S2_ATBD,
# S2opt=STATIC_S2_OPT, Combinaison=DYN_S2_RESCALED.
# Source: nc_genuine53_windcorr; macro = musica_in_Blois_pblh.nc; T at 1 m, no shift.
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)})
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]; invisible(lapply(src,source))
source("Chapter3_config.R")
FORC <- "in_files/musica_in_Blois_pblh.nc"; Z <- 1
BASE <- "out_files/Chapter3/nc_genuine53_windcorr"
MAP  <- c(LiDARfixe="STATIC_ALS", S2seul="DYN_S2_ATBD", S2opt="STATIC_S2_OPT", Combinaison="DYN_S2_RESCALED")
WINS <- list(apriloct=c("2021-04-01","2021-10-31"), junsep=c("2021-06-01","2021-09-30"))

nc <- nc_open(FORC); tu <- ncatt_get(nc,"time","units")$value; th <- ncvar_get(nc,"time")
t0 <- as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH <- data.table(time=floor_date(t0+th*3600,"hour"), Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15); nc_close(nc)
macH <- macH[,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]

# ---- observed micro (HOBO 1 m) ----
hob <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hob[, time := floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hob <- hob[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
HO  <- hob[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
IDS <- sort(unique(HO$id_plot))

# slope and dTmax on a window from an hourly (time,Tmic)+(Tmac) merged table
metrics_win <- function(dt, w){
  d <- dt[as.Date(time)>=as.Date(w[1]) & as.Date(time)<=as.Date(w[2])]
  if (nrow(d) < 50) return(list(slope=NA_real_, dtmax=NA_real_))
  sl <- as.numeric(coef(lm(Tmic~Tmac,d))[2])
  dd <- d[,.(mx_mic=max(Tmic,na.rm=TRUE), mx_mac=max(Tmac,na.rm=TRUE)), by=.(date=as.Date(time))]
  list(slope=sl, dtmax=mean(dd$mx_mic-dd$mx_mac, na.rm=TRUE))
}

# per-plot observed
OBS <- rbindlist(lapply(IDS,function(id){ m <- merge(HO[id_plot==id],macH,by="time")
  a<-metrics_win(m,WINS$apriloct); j<-metrics_win(m,WINS$junsep)
  data.table(id_plot=id, obs_slope_apriloct=a$slope, obs_dtmax_apriloct=a$dtmax,
             obs_slope_junsep=j$slope, obs_dtmax_junsep=j$dtmax) }))

# per-plot simulated, each scenario
SIM <- rbindlist(lapply(names(MAP),function(scn){
  rbindlist(lapply(IDS,function(id){
    f <- file.path(BASE,MAP[scn],sprintf("musica_out_HOBO_%s.nc",id))
    m <- if (file.exists(f)) tryCatch(micro_hourly_at(f,Z),error=function(e)NULL) else NULL
    if (is.null(m)) return(data.table(scenario=scn,id_plot=id,sim_slope_apriloct=NA_real_,sim_dtmax_apriloct=NA_real_,sim_slope_junsep=NA_real_,sim_dtmax_junsep=NA_real_))
    mm <- merge(m,macH,by="time"); a<-metrics_win(mm,WINS$apriloct); j<-metrics_win(mm,WINS$junsep)
    data.table(scenario=scn,id_plot=id, sim_slope_apriloct=a$slope, sim_dtmax_apriloct=a$dtmax,
               sim_slope_junsep=j$slope, sim_dtmax_junsep=j$dtmax) })) }))

P <- unique(fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/reframe_perplot_scenarios.csv")[,.(id_plot,P)])
D <- merge(merge(SIM,OBS,by="id_plot"),P,by="id_plot",all.x=TRUE)
fwrite(D,"out_files/Chapter3/tables/junsep_full_perplot.csv")

# ---- long form over the 4 (metric, window) combinations, then group ----
combos <- CJ(met=c("slope","dtmax"), win=c("apriloct","junsep"))
LONG <- rbindlist(lapply(seq_len(nrow(combos)), function(k){
  met <- combos$met[k]; win <- combos$win[k]
  D[, .(scenario, id_plot, P, metric=met, window=win,
        sim=get(sprintf("sim_%s_%s",met,win)), obs=get(sprintf("obs_%s_%s",met,win)))]
}))

# ---- by-archetype summary ----
SUMM <- LONG[is.finite(sim) & is.finite(obs),
  .(n=.N, mean_sim=round(mean(sim),3),
    r=if(.N>3) round(cor(sim,obs),2) else NA_real_,
    bias=round(mean(sim-obs),3), sd_res=round(sd(sim-obs),3)),
  by=.(window,metric,scenario,P)][order(window,metric,scenario,P)]
fwrite(SUMM,"out_files/Chapter3/tables/junsep_full_by_archetype.csv")

# ---- paired contrasts vs LiDARfixe (bootstrap CI) ----
set.seed(1); B <- 2000
paired <- function(x){ x<-x[is.finite(x)]; if(length(x)<3) return(c(NA,NA,NA))
  bs<-replicate(B, mean(sample(x,length(x),TRUE))); c(mean(x), unname(quantile(bs,c(.025,.975)))) }
CON <- rbindlist(lapply(seq_len(nrow(combos)), function(k){
  met <- combos$met[k]; win <- combos$win[k]
  W <- dcast(LONG[metric==met & window==win], id_plot+P~scenario, value.var="sim")
  rbindlist(lapply(setdiff(names(MAP),"LiDARfixe"), function(scn){
    W[, { s <- paired(get(scn)-LiDARfixe)
      .(window=win, metric=met, contrast=paste0(scn,"-LiDARfixe"),
        n=sum(is.finite(get(scn)-LiDARfixe)), mean=round(s[1],3), lo=round(s[2],3), hi=round(s[3],3)) },
      by=P] }))
}))
fwrite(CON,"out_files/Chapter3/tables/junsep_full_paired_contrasts.csv")

cat("\n=== SLOPE ranking r, by scenario x archetype (apriloct -> junsep) ===\n")
S<-SUMM[metric=="slope"]; W<-dcast(S,scenario+P~window,value.var="r")
print(W[order(scenario,P)])
cat("\nwrote junsep_full_perplot.csv / _by_archetype.csv / _paired_contrasts.csv\n")
