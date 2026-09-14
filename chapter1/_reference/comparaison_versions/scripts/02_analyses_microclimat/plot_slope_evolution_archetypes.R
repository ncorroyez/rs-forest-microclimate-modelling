# ==============================================================================
# Per-scenario evolution of the slope-and-equilibrium coupling slope
# beta = coef(lm(T_micro ~ T_macro)) on ALL hourly temperatures (Gril 2023; no
# time shift; beta<1 = buffering), computed PER MONTH, for the 4 archetype HOBO
# plots (P1 open -> P4 dense), both binaries (v3.2.0 + v3.2.3 iter), the 8
# monthly scenarios. Sim micro = Tair at 1 m (get_tair_at_z, identical to
# cmp_v320_v323iter_ch3.R); macro = free-air forcing Tair. Observed slope from the
# HOBO loggers is the black reference. Window Apr-Oct so the seasonal evolution
# (leaf-out tightens coupling toward summer, relaxes at senescence) is visible.
#   Rscript plot_slope_evolution_archetypes.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr); library(data.table)
  library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

WIN <- as.Date(c("2021-04-01","2021-10-31"))
mo  <- function(d) factor(month.abb[month(d)], levels=month.abb)
beta <- function(y,x){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<24) return(NA_real_)
  as.numeric(coef(lm(y[ok]~x[ok]))[2]) }
SCN <- c("STATIC_ALS","STATIC_ALS_DOPT","DYN_S2_ATBD_NM","DYN_S2_OPT_NM",
         "DYN_ATBD_rfull_NM","DYN_OPT_rdopt_NM","DYN_ALS_S2TIMING_NM","DYN_ALSoOPT_S2TIMING_NM")

# ---- macro hourly (free air) -------------------------------------------------
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2]][, hr:=floor_date(time,"hour")]

# ---- 4 archetype HOBO plots (same selection as the time-series figure) -------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time:=as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
cand <- df[id_plot %in% unique(as.character(hb$id_plot)) & !(id_plot %in% CFG_C3$ids_to_remove)][order(LAI_ALS)]
qi  <- round(quantile(seq_len(nrow(cand)), c(0.10,0.40,0.65,0.90)))
sel <- cand[qi][, plab := sprintf("P%d  (LAI=%.1f)", .I, LAI_ALS)]
cat("Archetype plots:\n"); print(sel[,.(id_plot,LAI_ALS=round(LAI_ALS,2),plab)])

# ---- observed monthly slope --------------------------------------------------
ho <- hb[id_plot %in% sel$id_plot & as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2],
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"))]
ho <- merge(ho, macro_h[,.(hr,Tmacro)], by="hr")
ho[, mon:=mo(as.Date(hr))]
OBS <- merge(ho[, .(slope=beta(Tobs,Tmacro)), by=.(id_plot,mon)], sel[,.(id_plot,plab)], by="id_plot")

# ---- simulated monthly slope, both binaries ----------------------------------
read_micro <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc)
  if(is.null(r)||!nrow(r))return(NULL); r<-as.data.table(r)[as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2]]
  r[, hr:=floor_date(time,"hour")][, .(hr, Tsim=Tair_sim)] }
roots <- c("v3.2.0"=file.path(CFG_C3$out_dir,"nc"), "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter"))
SIM <- rbindlist(lapply(names(roots), function(v) rbindlist(lapply(seq_len(nrow(sel)), function(k){
  id<-sel$id_plot[k]; rbindlist(lapply(SCN, function(s){
    f<-file.path(roots[v],s,sprintf("musica_out_HOBO_%s.nc",id)); if(!file.exists(f))return(NULL)
    m<-read_micro(f); if(is.null(m))return(NULL); m<-merge(m, macro_h[,.(hr,Tmacro)], by="hr")
    m[, .(slope=beta(Tsim,Tmacro)), by=.(mon=mo(as.Date(hr)))][, `:=`(version=v, plab=sel$plab[k], scenario=s)]
  })) }))))

# ---- family / variant encoding (4 colours x 2 linetypes) ---------------------
FAM <- c(STATIC_ALS="LiDAR", STATIC_ALS_DOPT="LiDAR", DYN_S2_ATBD_NM="Sentinel-2",
         DYN_S2_OPT_NM="Sentinel-2", DYN_ATBD_rfull_NM="S2 rescaled->LiDAR",
         DYN_OPT_rdopt_NM="S2 rescaled->LiDAR", DYN_ALS_S2TIMING_NM="Fusion (ALS x S2-timing)",
         DYN_ALSoOPT_S2TIMING_NM="Fusion (ALS x S2-timing)")
VAR <- c(STATIC_ALS="full", STATIC_ALS_DOPT="d_opt", DYN_S2_ATBD_NM="ATBD", DYN_S2_OPT_NM="opt",
         DYN_ATBD_rfull_NM="ATBD", DYN_OPT_rdopt_NM="opt", DYN_ALS_S2TIMING_NM="ALS",
         DYN_ALSoOPT_S2TIMING_NM="ALSxopt")
PAL_FAM <- c("LiDAR"="#0072B2","Sentinel-2"="#D55E00","S2 rescaled->LiDAR"="#E69F00",
             "Fusion (ALS x S2-timing)"="#CC79A7")
LT <- c("full"="solid","d_opt"="22","ATBD"="solid","opt"="22","ALS"="solid","ALSxopt"="22")
SIM[, family:=factor(FAM[scenario], levels=names(PAL_FAM))][, variant:=VAR[scenario]]
SIM[, version:=factor(version, levels=c("v3.2.0","v3.2.3 iter"))]

OUT_TAB <- "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables"
OUT_FIG <- "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/figures"
fwrite(SIM, file.path(OUT_TAB,"slope_evolution_archetypes_sim.csv"))
fwrite(OBS, file.path(OUT_TAB,"slope_evolution_archetypes_obs.csv"))

g <- ggplot(SIM[!is.na(slope)], aes(mon, slope, colour=family, linetype=variant, group=scenario)) +
  geom_hline(yintercept=1, colour="grey55", linetype="dashed", linewidth=0.4) +
  geom_line(data=OBS[!is.na(slope)], aes(mon,slope,group=1), colour="black", linewidth=0.7, inherit.aes=FALSE) +
  geom_point(data=OBS[!is.na(slope)], aes(mon,slope), colour="black", size=1.4, inherit.aes=FALSE) +
  geom_line(linewidth=0.6) + geom_point(size=1.3) +
  facet_grid(version ~ plab) +
  scale_colour_manual(values=PAL_FAM, name="Scenario family") +
  scale_linetype_manual(values=LT, guide="none") +
  labs(x="2021 (month)", y="Micro-macro coupling slope  beta = dT_micro / dT_macro",
       subtitle=paste0("Monthly slope-and-equilibrium beta (all hourly temperatures, no shift). ",
                       "Dashed grey = beta=1 (no buffering); below = buffering.",
                       "\nBlack = observed HOBO. Linetype: solid = full/ATBD/ALS, dashed = d_opt/opt. ",
                       "4 HOBO plots binned by LAI_ALS.")) +
  theme_article(10) + theme(legend.position="top")
ggsave_article(file.path(OUT_FIG,"Fig_slope_evolution_archetypes"), g, 10, 5.6)
cat("\nDONE -> Fig_slope_evolution_archetypes (.png/.pdf) + slope_evolution_archetypes_{sim,obs}.csv\n")
