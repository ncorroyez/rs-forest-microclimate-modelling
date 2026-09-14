# ==============================================================================
# Diagnostic time series: LAI(t) input + daily Tmax(t) per scenario, for 4 HOBO
# plots binned by LAI_ALS (P1 open -> P4 dense; these are HOBO loggers, NOT the
# cLHS Ch1 archetypes), under both binaries (v3.2.0 + v3.2.3 iter), the 8 monthly
# scenarios. LAI is recomputed from the EXACT mk_stat/mk_dyn used in
# run_v323iter_ch3_full.R (one-sided Leaf_area, seed/×2 dropped); temperature is
# read straight from nc Tair_z at 1 m. Raw S2 points are overlaid on the
# DYN_S2_ATBD curve as a keying check; macro free-air Tmax is the buffering ref.
#   Rscript plot_ts_lai_temp_archetypes.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

WIN  <- as.Date(c("2021-06-01","2021-09-30"))
LY   <- 2020:2022
OUT_FIG <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
OUT_TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"

# ---- rebuild df + S2 series + 8 scenarios (identical to run_v323iter_ch3_full) -
inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("seed row"); d$Julian_day<-366L; rbind(ph,d) }
mk_dyn  <- function(ts) function(pr){ pid<-plot_id_from_row(pr)
  if(!pid%in%names(ts)) return(NULL); inject_seed(make_phenology_from_s2(ts[[pid]],LY)) }
mk_stat <- function(col) function(pr) inject_seed(calc_phenology(list.year=LY,nleafage=1,
  budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,
  LAI_max_per_cohort=as.numeric(pr[[col]])))

prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE)
df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt]
df[, LAI_A_OPTCAL := LAI_ALS * LAI_S2_ATBD / LAI_S2_DOPT]
df[is.na(LAI_A_OPTCAL), LAI_A_OPTCAL := LAI_ALS]
rA  <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD); bod <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT)
ropt_p <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
fp <- setNames(df$LAI_ALS/df$LAI_S2_ATBD, df$pid); fA <- setNames(df$LAI_ALS/df$LAI_S2_DOPT, df$pid)

nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk <- rast(files); dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[,-1,drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
mkser <- function(scale_fun, fb_col) setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]; a<-atbd[[p]]
  if(is.null(a)) return(data.frame(doy=1:365,lai=rep(df[[fb_col]][i],365)))
  data.frame(doy=a$doy, lai=a$lai*scale_fun(i,p)) }), df$pid)
ts_atbd<-atbd; ts_opt<-mkser(function(i,p)ropt_p[[p]],"LAI_S2_DOPT"); ts_rfull<-mkser(function(i,p)rA,"LAI_ALS")
ts_ordopt<-mkser(function(i,p)ropt_p[[p]]*bod,"LAI_ALS_DOPT")
ts_fp<-mkser(function(i,p)fp[[p]],"LAI_ALS"); ts_fA<-mkser(function(i,p)fA[[p]],"LAI_A_OPTCAL")

.lf <- function(col) function(pr) as.numeric(pr[[col]])
SC <- list(
  STATIC_ALS              = mk_stat("LAI_ALS"),
  STATIC_ALS_DOPT         = mk_stat("LAI_ALS_DOPT"),
  DYN_S2_ATBD_NM          = mk_dyn(ts_atbd),
  DYN_S2_OPT_NM           = mk_dyn(ts_opt),
  DYN_ATBD_rfull_NM       = mk_dyn(ts_rfull),
  DYN_OPT_rdopt_NM        = mk_dyn(ts_ordopt),
  DYN_ALS_S2TIMING_NM     = mk_dyn(ts_fp),
  DYN_ALSoOPT_S2TIMING_NM = mk_dyn(ts_fA))
SCN <- names(SC)

# ---- choose 4 archetype HOBO plots by ascending LAI_ALS -----------------------
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove)]
hobo_ids <- unique(as.character(hb$id_plot))
cand <- df[id_plot %in% hobo_ids & !(id_plot %in% CFG_C3$ids_to_remove)][order(LAI_ALS)]
qi <- round(quantile(seq_len(nrow(cand)), c(0.10,0.40,0.65,0.90)))
sel <- cand[qi]
sel[, plab := sprintf("P%d  (LAI=%.1f)", .I, LAI_ALS)]
cat("Selected archetype HOBO plots:\n"); print(sel[,.(id_plot,pid,LAI_ALS=round(LAI_ALS,2),Hmax=round(Hmax,1),fCover=round(fCover,2),plab)])

# ---- LAI(t) per plot x scenario (one-sided Leaf_area, year 2021) --------------
lai_dt <- rbindlist(lapply(seq_len(nrow(sel)), function(k){
  prow <- as.data.frame(sel[k]); rbindlist(lapply(SCN, function(s){
    ph <- SC[[s]](prow); if(is.null(ph)) return(NULL); ph <- as.data.table(ph)
    lac <- grep("Leaf_area", names(ph), value=TRUE)[1]
    ph[, date := as.Date(Julian_day-1, origin=sprintf("%d-01-01", year))]
    ph[year==2021 & date>=WIN[1] & date<=WIN[2], .(plab=sel$plab[k], scenario=s, date, lai=get(lac))]
  })) }))

# raw S2 extraction points for the DYN_S2_ATBD keying check
raw_pts <- merge(long[date>=WIN[1] & date<=WIN[2]], sel[,.(pid,plab)], by.x="plot_id", by.y="pid")
raw_pts[, scenario := "DYN_S2_ATBD_NM"]

# ---- daily Tmax per plot x scenario x version + obs + macro -------------------
read_tmax <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc)
  if(is.null(r)||!nrow(r))return(NULL); r<-as.data.table(r)[, d:=as.Date(time)]
  r[d>=WIN[1] & d<=WIN[2], .(tmax=max(Tair_sim,na.rm=TRUE)), by=d] }
roots <- c("v3.2.0"=file.path(CFG_C3$out_dir,"nc"), "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter"))
tmp_dt <- rbindlist(lapply(names(roots), function(v) rbindlist(lapply(seq_len(nrow(sel)), function(k){
  id<-sel$id_plot[k]; rbindlist(lapply(SCN, function(s){
    f<-file.path(roots[v],s,sprintf("musica_out_HOBO_%s.nc",id)); if(!file.exists(f))return(NULL)
    t<-read_tmax(f); if(is.null(t))return(NULL); t[, .(version=v, plab=sel$plab[k], scenario=s, date=d, tmax)]
  })) }))))

# observed HOBO daily Tmax
obs <- hb[id_plot %in% sel$id_plot, .(tmax=max(t_hobo,na.rm=TRUE)), by=.(id_plot, d=as.Date(time))]
obs <- merge(obs, sel[,.(id_plot,plab)], by="id_plot")[d>=WIN[1] & d<=WIN[2], .(plab, date=d, tmax)]

# macro free-air daily Tmax (forcing)
ncf<-nc_open(CFG_C3$forcing_file); macro<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
  Tair=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro <- macro[, .(tmax=max(Tair,na.rm=TRUE)), by=.(date=as.Date(time))][date>=WIN[1] & date<=WIN[2]]

# ---- save underlying data -----------------------------------------------------
dir.create(OUT_TAB, recursive=TRUE, showWarnings=FALSE)
fwrite(lai_dt, file.path(OUT_TAB,"ts_archetypes_lai.csv"))
fwrite(tmp_dt, file.path(OUT_TAB,"ts_archetypes_tmax.csv"))
fwrite(obs,    file.path(OUT_TAB,"ts_archetypes_obs.csv"))

# ---- family / variant encoding (4 colours x 2 linetypes = 8 scenarios) --------
FAM <- c(STATIC_ALS="LiDAR", STATIC_ALS_DOPT="LiDAR",
         DYN_S2_ATBD_NM="Sentinel-2", DYN_S2_OPT_NM="Sentinel-2",
         DYN_ATBD_rfull_NM="S2 rescaled->LiDAR", DYN_OPT_rdopt_NM="S2 rescaled->LiDAR",
         DYN_ALS_S2TIMING_NM="Fusion (ALS x S2-timing)", DYN_ALSoOPT_S2TIMING_NM="Fusion (ALS x S2-timing)")
VAR <- c(STATIC_ALS="full", STATIC_ALS_DOPT="d_opt",
         DYN_S2_ATBD_NM="ATBD", DYN_S2_OPT_NM="opt",
         DYN_ATBD_rfull_NM="ATBD", DYN_OPT_rdopt_NM="opt",
         DYN_ALS_S2TIMING_NM="ALS", DYN_ALSoOPT_S2TIMING_NM="ALSxopt")
PAL_FAM <- c("LiDAR"="#0072B2","Sentinel-2"="#D55E00",
             "S2 rescaled->LiDAR"="#E69F00","Fusion (ALS x S2-timing)"="#CC79A7")
LT <- c("full"="solid","d_opt"="22","ATBD"="solid","opt"="22","ALS"="solid","ALSxopt"="22")
add_enc <- function(d){ d[, family:=factor(FAM[scenario],levels=names(PAL_FAM))]
  d[, variant:=VAR[scenario]]; d[] }
lai_dt<-add_enc(lai_dt); tmp_dt<-add_enc(tmp_dt)

# ---- assemble faceted long frame: 3 rows x 4 cols -----------------------------
ROWS <- c("LAI one-sided (m2/m2)", "Tmax @1m  v3.2.0 (degC)", "Tmax @1m  v3.2.3 iter (degC)")
L <- lai_dt[, .(plab, date, family, variant, scenario, row=ROWS[1], value=lai)]
T0 <- tmp_dt[version=="v3.2.0",      .(plab, date, family, variant, scenario, row=ROWS[2], value=tmax)]
T3 <- tmp_dt[version=="v3.2.3 iter", .(plab, date, family, variant, scenario, row=ROWS[3], value=tmax)]
main <- rbind(L,T0,T3); main[, row:=factor(row, levels=ROWS)]
# obs + macro replicated onto both temperature rows
obs2   <- rbind(obs[,.(plab,date,value=tmax,row=ROWS[2])], obs[,.(plab,date,value=tmax,row=ROWS[3])])
macro2 <- rbind(cbind(macro[,.(date,value=tmax)], row=ROWS[2]), cbind(macro[,.(date,value=tmax)], row=ROWS[3]))
macro2 <- macro2[rep(1:.N, each=nrow(sel))][, plab:=rep(sel$plab, .N/nrow(sel))]  # broadcast to all cols
obs2[, row:=factor(row,levels=ROWS)]; macro2[, row:=factor(row,levels=ROWS)]
raw_pts2 <- raw_pts[, .(plab, date, value=lai, row=factor(ROWS[1],levels=ROWS))]

g <- ggplot(main, aes(date, value)) +
  geom_line(data=macro2, aes(date,value), colour="grey60", linetype="dotted", linewidth=0.4, inherit.aes=FALSE) +
  geom_line(data=obs2,   aes(date,value), colour="black", linewidth=0.7, inherit.aes=FALSE) +
  geom_point(data=raw_pts2, aes(date,value), colour="#D55E00", size=0.7, alpha=0.5, inherit.aes=FALSE) +
  geom_line(aes(colour=family, linetype=variant, group=scenario), linewidth=0.55) +
  facet_grid(row ~ plab, scales="free_y", switch="y") +
  scale_colour_manual(values=PAL_FAM, name="Scenario family") +
  scale_linetype_manual(values=LT, guide="none") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  labs(x="2021", y=NULL,
       subtitle=paste0("LAI input (top) and daily Tmax at 1 m per scenario, 4 HOBO plots binned by LAI_ALS (P1 open -> P4 dense).",
                       "\nBlack = observed HOBO; grey dotted = free-air (macro) Tmax; orange points = raw S2 LAI (DYN_S2_ATBD keying check).",
                       "\nLinetype: solid = full/ATBD/ALS variant, dashed = d_opt/opt variant within each family.")) +
  theme_article(10) + theme(legend.position="top", strip.placement="outside",
                            strip.text.y.left=element_text(angle=90))
dir.create(OUT_FIG, recursive=TRUE, showWarnings=FALSE)
ggsave_article(file.path(OUT_FIG,"Fig_ts_lai_temp_archetypes"), g, 11, 7.2)

# ==============================================================================
# Variant 2 — buffering anomaly ΔTmax = Tmax_sim − Tmax_macro (7-day rolling),
# which removes the 25 degC weather swing and reveals the ~0.4-1.5 degC
# per-scenario / per-plot separation buried in the raw Tmax band.
# ==============================================================================
roll <- function(v) data.table::frollmean(v, 7, align="center", na.rm=TRUE)
dM <- merge(tmp_dt, macro[,.(date, tmac=tmax)], by="date")
dM[, danom := tmax - tmac]
setorder(dM, version, scenario, plab, date)
dM[, danom_s := roll(danom), by=.(version, scenario, plab)]
oM <- merge(obs, macro[,.(date, tmac=tmax)], by="date")[, danom:=tmax-tmac][order(plab,date)]
oM[, danom_s := roll(danom), by=plab]

A0 <- dM[version=="v3.2.0",      .(plab, date, family, variant, scenario, row=ROWS[2], value=danom_s)]
A3 <- dM[version=="v3.2.3 iter", .(plab, date, family, variant, scenario, row=ROWS[3], value=danom_s)]
amain <- rbind(L, A0, A3); amain[, row:=factor(row, levels=ROWS)]
aobs <- rbind(oM[,.(plab,date,value=danom_s,row=ROWS[2])], oM[,.(plab,date,value=danom_s,row=ROWS[3])])
aobs[, row:=factor(row, levels=ROWS)]
zline <- data.table(row=factor(c(ROWS[2],ROWS[3]), levels=ROWS))

ga <- ggplot(amain, aes(date, value)) +
  geom_hline(data=zline, aes(yintercept=0), colour="grey60", linetype="dotted", linewidth=0.4) +
  geom_line(data=aobs, aes(date,value), colour="black", linewidth=0.7, inherit.aes=FALSE) +
  geom_point(data=raw_pts2, aes(date,value), colour="#D55E00", size=0.7, alpha=0.5, inherit.aes=FALSE) +
  geom_line(aes(colour=family, linetype=variant, group=scenario), linewidth=0.6) +
  facet_grid(row ~ plab, scales="free_y", switch="y") +
  scale_colour_manual(values=PAL_FAM, name="Scenario family") +
  scale_linetype_manual(values=LT, guide="none") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  labs(x="2021", y=NULL,
       subtitle=paste0("LAI input (top) and buffering anomaly dTmax = Tmax_sim - Tmax_macro at 1 m (7-day rolling), per scenario.",
                       "\nBelow grey-dotted 0 = cooler than free air (buffering). Black = observed HOBO dTmax; orange = raw S2 LAI.",
                       "\nLinetype: solid = full/ATBD/ALS, dashed = d_opt/opt within each family. 4 HOBO plots binned by LAI_ALS.")) +
  theme_article(10) + theme(legend.position="top", strip.placement="outside",
                            strip.text.y.left=element_text(angle=90))
ggsave_article(file.path(OUT_FIG,"Fig_ts_lai_dTmax_archetypes"), ga, 11, 7.2)

# ==============================================================================
# Variant 3 — MONTHLY means (x = month). Collapses the daily weather noise so the
# per-scenario separation is fully legible; LAI monthly mean (top) + monthly-mean
# buffering anomaly dTmax per binary (rows). One point per month per scenario.
# ==============================================================================
mo <- function(d) factor(month.abb[month(d)], levels=month.abb)
lai_m <- lai_dt[, .(value=mean(lai)), by=.(plab, family, variant, scenario, mon=mo(date))][, row:=ROWS[1]]
d0_m  <- dM[version=="v3.2.0",      .(value=mean(danom,na.rm=TRUE)), by=.(plab,family,variant,scenario,mon=mo(date))][, row:=ROWS[2]]
d3_m  <- dM[version=="v3.2.3 iter", .(value=mean(danom,na.rm=TRUE)), by=.(plab,family,variant,scenario,mon=mo(date))][, row:=ROWS[3]]
mmain <- rbind(lai_m, d0_m, d3_m); mmain[, row:=factor(row, levels=ROWS)]
obs_m <- oM[, .(value=mean(danom,na.rm=TRUE)), by=.(plab, mon=mo(date))]
obs_m2 <- rbind(copy(obs_m)[, row:=ROWS[2]], copy(obs_m)[, row:=ROWS[3]])[, row:=factor(row,levels=ROWS)]
fwrite(mmain, file.path(OUT_TAB,"ts_archetypes_monthly.csv"))

gm <- ggplot(mmain, aes(mon, value, colour=family, linetype=variant, group=scenario)) +
  geom_hline(data=zline, aes(yintercept=0), colour="grey60", linetype="dotted", linewidth=0.4, inherit.aes=FALSE) +
  geom_line(data=obs_m2, aes(mon,value,group=1), colour="black", linewidth=0.7, inherit.aes=FALSE) +
  geom_point(data=obs_m2, aes(mon,value), colour="black", size=1.4, inherit.aes=FALSE) +
  geom_line(linewidth=0.6) + geom_point(size=1.3) +
  facet_grid(row ~ plab, scales="free_y", switch="y") +
  scale_colour_manual(values=PAL_FAM, name="Scenario family") +
  scale_linetype_manual(values=LT, guide="none") +
  labs(x="2021 (month)", y=NULL,
       subtitle=paste0("Monthly means: LAI input (top) and buffering anomaly dTmax = Tmax_sim - Tmax_macro at 1 m, per scenario.",
                       "\nBelow grey-dotted 0 = cooler than free air. Black = observed HOBO dTmax. ",
                       "Linetype: solid = full/ATBD/ALS, dashed = d_opt/opt. 4 HOBO plots binned by LAI_ALS.")) +
  theme_article(10) + theme(legend.position="top", strip.placement="outside",
                            strip.text.y.left=element_text(angle=90))
ggsave_article(file.path(OUT_FIG,"Fig_monthly_lai_dTmax_archetypes"), gm, 10, 7.2)

cat("\nDONE -> Fig_ts_lai_temp_archetypes + Fig_ts_lai_dTmax_archetypes + Fig_monthly_lai_dTmax_archetypes (.png/.pdf)\n")
cat(sprintf("LAI rows=%d  Tmax rows=%d  obs rows=%d | per-scenario daily-Tmax SD=%.2f degC\n",
            nrow(lai_dt), nrow(tmp_dt), nrow(obs), mean(tmp_dt[,sd(tmax),by=.(version,plab,date)]$V1,na.rm=TRUE)))
