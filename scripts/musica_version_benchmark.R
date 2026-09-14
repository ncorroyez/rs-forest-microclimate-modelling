# ==============================================================================
# MuSICA version benchmark — STANDALONE, reproducible, shareable.
# Compares two MuSICA runs (e.g. v3.2.0 vs v3.2.3) against HOBO loggers and
# against each other, on several temperature metrics, day subsets, and skill
# scores. No external project files — everything is configured below.
#
# ------------------------------ REQUIRED INPUTS -------------------------------
# Prepare and ship these three things (paths set in the CONFIG block):
#
# 1) HOBO_CSV — CSV of understory air-temperature loggers. Needed columns:
#       id_plot          plot id, matching the MuSICA file names (e.g. "41_01")
#       position_sensor  sensor position; we keep "a" (the ~1 m air sensor)
#       datetime         timestamp "YYYY-MM-DD HH:MM:SS" (UTC)
#       t_hobo           air temperature (°C)
#
# 2) SIM — one FOLDER per MuSICA version. Each folder holds one NetCDF per plot,
#    named  musica_out_HOBO_<id_plot>.nc , with variables:
#       Tair_z(time, nair)      air temperature profile (K)
#       relative_height(nair)   layer height / canopy height (-)
#       veget_height_top(time)  canopy height (m)
#       time                    "hours since YYYY-MM-DD ..."
#    Use the SAME coalition (e.g. full-canopy REF) in both folders, same plots.
#
# 3) FORCING — NetCDF macroclimate forcing (above-canopy reference). Needs:
#       Tair (K)  and  time ("hours since ..." or "seconds since ...")
#
# All inputs must cover the same plots and the DATES period.
# Packages: data.table, ncdf4, ggplot2, lubridate
# Run:  Rscript musica_version_benchmark.R
# ==============================================================================

suppressMessages({ library(data.table); library(ncdf4); library(ggplot2); library(lubridate) })

## ------------------------------ CONFIG (edit) --------------------------------
HOBO_CSV <- "in_files/Blois_data_temperature.csv"
SIM      <- c("v3.2.0" = "out_files/musica_hobo_z05/1111",        # MuSICA run A (folder of NetCDFs)
              "v3.2.3" = "out_files/musica_hobo_z05_v323/1111")   # MuSICA run B
FORCING  <- "in_files/musica_in_Blois.nc"
IDS_REMOVE <- c("41_13","41_14","41_20","41_41","41_50","41_51","41_53")  # plots to drop ("" for none)
DATES    <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day") # study period
Z        <- 1.0     # extraction height (m), interpolated — match the HOBO sensor height
HOT_Q    <- 0.90    # hottest-day quantile (on macro Tmax)
OUTDIR   <- "outputs/figures_pipeline/annex"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---------------- helper: interpolate hourly Tair to height Z ----------------
tair_at_Z <- function(path, z = Z) {
  nc <- tryCatch(nc_open(path), error = function(e) NULL); if (is.null(nc)) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec <- floor_date(t0 + dhours(ncvar_get(nc, "time")), "hour")
  Tk <- ncvar_get(nc, "Tair_z")                       # [nair, time]
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh
  if (z <= zl[1]) { i<-1L; j<-1L; w<-0 } else if (z >= zl[length(zl)]) { i<-length(zl); j<-i; w<-0 } else {
    i <- max(which(zl <= z)); j <- i+1L; w <- (z - zl[i])/(zl[j]-zl[i]) }
  data.table(time = tvec, t = ((1-w)*Tk[i,] + w*Tk[j,]) - 273.15)
}
id_from <- function(f) sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
sim_daily <- function(dir) rbindlist(lapply(list.files(dir,"\\.nc$",full.names=TRUE), function(f){
  d <- tair_at_Z(f); if (is.null(d)) return(NULL)
  dd <- d[as.Date(time) %in% DATES, .(Tmin=min(t),Tmax=max(t),Tmean=mean(t)), by=.(date=as.Date(time))]
  dd$id_plot <- id_from(f); dd }))
sim_hourly <- function(dir) rbindlist(lapply(list.files(dir,"\\.nc$",full.names=TRUE), function(f){
  d <- tair_at_Z(f); if (is.null(d)) return(NULL); d <- d[as.Date(time) %in% DATES]
  d$id_plot <- id_from(f); d }))

## ---------------------------- read inputs ------------------------------------
# HOBO -> hourly + daily (1 m sensor)
H <- fread(HOBO_CSV)
H[, datetime := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
H <- H[position_sensor=="a" & as.Date(datetime) %in% DATES & !id_plot %in% IDS_REMOVE]
H[, time := floor_date(datetime, "hour")]
hobo_h <- H[, .(t = mean(t_hobo, na.rm=TRUE)), by=.(id_plot, time)]
hobo_d <- hobo_h[, .(Tmin=min(t),Tmax=max(t),Tmean=mean(t)), by=.(id_plot, date=as.Date(time))]

# macro daily from forcing
ncf <- nc_open(FORCING); tu <- ncatt_get(ncf,"time","units")$value; tt <- ncvar_get(ncf,"time")
t0  <- if (max(tt) > 1e5) as.POSIXct(sub("seconds since ","",tu),tz="UTC")+tt else as.POSIXct(sub("hours since ","",tu),tz="UTC")+dhours(tt)
macro_d <- data.table(date=as.Date(floor_date(t0,"hour")), t=ncvar_get(ncf,"Tair")-273.15)[date %in% DATES,
            .(Tmin_M=min(t),Tmax_M=max(t),Tmean_M=mean(t),Amp_M=max(t)-min(t)), by=date]; nc_close(ncf)
hot_days <- macro_d[Tmax_M >= quantile(Tmax_M, HOT_Q, na.rm=TRUE), date]
message(sprintf("Plots HOBO=%d | hottest 10%%: %d/%d days (>= %.1f°C)",
                uniqueN(hobo_d$id_plot), length(hot_days), nrow(macro_d), quantile(macro_d$Tmax_M,HOT_Q,na.rm=TRUE)))

## --------- per (plot,date) metric values for one version vs HOBO -------------
metrics_daily <- function(simd){
  s <- copy(simd); setnames(s, c("Tmin","Tmax","Tmean"), c("Tmin_s","Tmax_s","Tmean_s"))
  o <- copy(hobo_d); setnames(o, c("Tmin","Tmax","Tmean"), c("Tmin_o","Tmax_o","Tmean_o"))
  m <- merge(merge(s, o, by=c("id_plot","date")), macro_d, by="date")
  data.table(id_plot=m$id_plot, date=m$date,
    Tmax_offset_sim=m$Tmax_s-m$Tmax_M,   Tmax_offset_obs=m$Tmax_o-m$Tmax_M,
    Tmin_offset_sim=m$Tmin_s-m$Tmin_M,   Tmin_offset_obs=m$Tmin_o-m$Tmin_M,
    Tmean_offset_sim=m$Tmean_s-m$Tmean_M,Tmean_offset_obs=m$Tmean_o-m$Tmean_M,
    amplitude_sim=m$Tmax_s-m$Tmin_s,     amplitude_obs=m$Tmax_o-m$Tmin_o,
    buffering_sim=m$Amp_M-(m$Tmax_s-m$Tmin_s), buffering_obs=m$Amp_M-(m$Tmax_o-m$Tmin_o)) }

DMETS <- c("Tmax_offset","Tmin_offset","Tmean_offset","amplitude","buffering")
MET   <- lapply(SIM, function(dir) metrics_daily(sim_daily(dir)))
names(MET) <- names(SIM)

## ---------------------- assemble long obs/sim tables -------------------------
# (a) daily, vs HOBO, all + hot
long_hobo <- rbindlist(lapply(names(SIM), function(v) rbindlist(lapply(DMETS, function(mt)
  rbindlist(lapply(c("all","hot10"), function(per){
    d <- MET[[v]]; if (per=="hot10") d <- d[date %in% hot_days]
    data.table(version=v, metric=mt, period=per, obs=d[[paste0(mt,"_obs")]], sim=d[[paste0(mt,"_sim")]]) }))))))
# (b) hourly absolute Tair, vs HOBO (all period)
SH <- lapply(SIM, sim_hourly); names(SH) <- names(SIM)
long_hourly <- rbindlist(lapply(names(SIM), function(v){
  m <- merge(SH[[v]], hobo_h, by=c("id_plot","time"))
  data.table(version=v, metric="Tair_hourly", period="all", obs=m$t.y, sim=m$t.x) }))
LONG <- rbind(long_hobo, long_hourly)

## ----------------------------- skill scores ---------------------------------
perf <- function(x,y){ ok<-is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]
  if (length(x)<3) return(data.table(r=NA,RMSE=NA,bias=NA,MAE=NA,n=length(x)))
  data.table(r=cor(x,y), RMSE=sqrt(mean((y-x)^2)), bias=mean(y-x), MAE=mean(abs(y-x)), n=length(x)) }

# vs-HOBO scores
SC_hobo <- LONG[, cbind(axis=paste(version,"vs HOBO"), perf(obs,sim)), by=.(metric,period,version)][, version:=NULL]
# model-model scores (same plot×date / plot×time)
mm_daily <- rbindlist(lapply(DMETS, function(mt) rbindlist(lapply(c("all","hot10"), function(per){
  a<-MET[[1]]; b<-MET[[2]]; if (per=="hot10"){a<-a[date%in%hot_days]; b<-b[date%in%hot_days]}
  m<-merge(a[,.(id_plot,date,A=get(paste0(mt,"_sim")))], b[,.(id_plot,date,B=get(paste0(mt,"_sim")))], by=c("id_plot","date"))
  cbind(metric=mt, period=per, axis=paste(names(SIM)[1],"vs",names(SIM)[2]), perf(m$A,m$B)) }))))
mmh <- merge(SH[[1]][,.(id_plot,time,A=t)], SH[[2]][,.(id_plot,time,B=t)], by=c("id_plot","time"))
mm_hourly <- cbind(metric="Tair_hourly", period="all", axis=paste(names(SIM)[1],"vs",names(SIM)[2]), perf(mmh$A,mmh$B))
SCORES <- rbind(SC_hobo, mm_daily, mm_hourly)
SCORES[, (c("r","RMSE","bias","MAE")) := lapply(.SD, round, 3), .SDcols=c("r","RMSE","bias","MAE")]
fwrite(SCORES, file.path(OUTDIR,"tab_musica_version_benchmark.csv"))
cat("\n=== Benchmark scores ===\n"); print(SCORES[order(metric,period,axis)])

## ----------------------------- FIG 1: heatmap -------------------------------
S <- copy(SCORES); S[, mp := factor(paste(metric,period))]
S[, axis := factor(axis, levels=c(paste(names(SIM)[1],"vs HOBO"), paste(names(SIM)[2],"vs HOBO"), paste(names(SIM)[1],"vs",names(SIM)[2])))]
ph <- ggplot(S, aes(axis, mp, fill=r)) +
  geom_tile(colour="white", linewidth=1) +
  geom_text(aes(label=sprintf("r=%.2f\nRMSE=%.2f\nbias=%+.2f", r, RMSE, bias)), size=3) +
  scale_fill_gradient2(low="#B2182B", mid="#F7F7F7", high="#2166AC", midpoint=0.6, limits=c(0,1), na.value="grey90") +
  labs(x=NULL,y=NULL,title="MuSICA version benchmark — skill heatmap",
       caption=sprintf("Offsets = micro−macro daily; amplitude=Tmax−Tmin; buffering=Amp_macro−Amp_micro; Tair_hourly=absolute. Height %.1f m.",Z)) +
  theme_minimal(base_size=11) + theme(axis.text=element_text(face="bold"), panel.grid=element_blank(),
       plot.caption=element_text(size=8,colour="grey40"))
ggsave(file.path(OUTDIR,"fig_musica_version_benchmark_heatmap.png"), ph, width=9, height=11, dpi=300, bg="white")

## ----------------------------- FIG 2: scatterplots --------------------------
# sim vs HOBO, all-period, metric (rows) × version (cols)
SP <- LONG[period=="all"]
SP[, metric := factor(metric, levels=c(DMETS,"Tair_hourly"))]
ann <- SP[, { ok<-is.finite(obs)&is.finite(sim)
  .(lab=sprintf("italic(r)==%.2f~~RMSE==%.2f", cor(obs[ok],sim[ok]), sqrt(mean((sim[ok]-obs[ok])^2)))) }, by=.(metric,version)]
ps <- ggplot(SP, aes(obs,sim)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50") +
  geom_hex(bins=55) + scale_fill_viridis_c(trans="log10", guide="none") +
  geom_text(data=ann, aes(x=-Inf,y=Inf,label=lab), parse=TRUE, hjust=-0.06, vjust=1.4, size=3.2, inherit.aes=FALSE) +
  facet_grid(metric ~ version, scales="free") +
  labs(x="Observed (HOBO)", y="Simulated (MuSICA)", title="Sim vs HOBO — all metrics × versions (all period)") +
  theme_bw(base_size=12) + theme(strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
       panel.grid.minor=element_blank())
ggsave(file.path(OUTDIR,"fig_musica_version_benchmark_scatter.png"), ps, width=8.5, height=15, dpi=200, bg="white")

cat(sprintf("\nSaved: %s/{tab_musica_version_benchmark.csv, fig_..._heatmap.png, fig_..._scatter.png}\n", OUTDIR))
