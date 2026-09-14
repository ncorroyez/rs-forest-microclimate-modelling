# ==============================================================================
# Sim vs HOBO daily validation for ΔTmin, diurnal amplitude, and buffering
# (amplitude attenuation), REF on z05+×2, BOTH MuSICA versions (v3.2.0 / v3.2.3).
# All micro temps at FIXED 1 m (interp). Metrics per plot×day:
#   ΔTmin   = Tmin_micro − Tmin_macro
#   Amp     = Tmax_micro − Tmin_micro            (diurnal range)
#   Buffer  = Amp_macro − Amp_micro              (>0 = forest damps the range)
# Output: outputs/figures_pipeline/annex/fig_sim_vs_obs_tmin_amp_buff.{png,pdf}
# ==============================================================================

suppressMessages({ library(data.table); library(here); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools); library(patchwork) })
source(here::here("R/config.R")); source(here::here("R/io.R")); source(here::here("R/musica.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
VERS <- list("v3.2.0 (legacy)"="out_files/musica_hobo_z05/1111", "v3.2.3"="out_files/musica_hobo_z05_v323/1111")

# ---- macro daily Tmin/Tmax from forcing ------------------------------------
nc <- nc_open(CFG$forcing_file)
tu <- ncatt_get(nc,"time","units")$value; tt<-ncvar_get(nc,"time")
t0 <- if (max(tt)>1e5) as.POSIXct(sub("seconds since ","",tu),tz="UTC")+tt else as.POSIXct(sub("hours since ","",tu),tz="UTC")+dhours(tt)
tair_macro <- ncvar_get(nc,"Tair")-273.15; nc_close(nc)
macro <- data.table(date=as.Date(floor_date(t0,"hour")), t=tair_macro)[date %in% CFG$date_seq,
          .(Tmin_macro=min(t), Tmax_macro=max(t)), by=date][, Amp_macro:=Tmax_macro-Tmin_macro]

# ---- HOBO daily Tmin/Tmax @1 m ---------------------------------------------
hobo <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")) %>%
  filter(position_sensor=="a", as.Date(datetime) %in% CFG$date_seq, !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time=floor_date(datetime,"hour")) %>% group_by(id_plot,time) %>%
  summarise(t=mean(t_hobo,na.rm=TRUE),.groups="drop") %>% as.data.table()
obs <- hobo[, .(Tmin_obs=min(t), Tmax_obs=max(t)), by=.(id_plot, date=as.Date(time))]

# ---- sim daily Tmin/Tmax @1 m (interp) per plot, per version ----------------
sim_daily <- function(dir){ rows<-list()
  for (f in list.files(here::here(dir),"\\.nc$",full.names=TRUE)){
    id<-sub("musica_out_HOBO_(.+)\\.nc$","\\1",basename(f))
    nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc)) next
    r<-tryCatch(get_tair_at_z(nc,z_target=1.0),error=function(e)NULL); nc_close(nc)
    if(is.null(r)||nrow(r)==0) next
    d<-as.data.table(r)[, .(date=as.Date(floor_date(time,"hour")), Tair_sim)][date %in% CFG$date_seq,
        .(Tmin_sim=min(Tair_sim), Tmax_sim=max(Tair_sim)), by=date]; d$id_plot<-id; rows[[id]]<-d }
  rbindlist(rows) }

fitlab <- function(x,y){ ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  sprintf("italic(r)==%.2f~~RMSE==%.2f~~bias==%+.2f", cor(x,y), sqrt(mean((y-x)^2)), mean(y-x)) }

LONG<-list()
for (v in names(VERS)){
  s <- sim_daily(VERS[[v]])
  m <- Reduce(function(a,b) merge(a,b,by=intersect(names(a),names(b))),
              list(s, obs, macro))   # by id_plot+date, then date
  m <- merge(merge(s, obs, by=c("id_plot","date")), macro, by="date")
  m[, `:=`(dTmin_sim=Tmin_sim-Tmin_macro, dTmin_obs=Tmin_obs-Tmin_macro,
           Amp_sim=Tmax_sim-Tmin_sim,     Amp_obs=Tmax_obs-Tmin_obs)]
  m[, `:=`(Buf_sim=Amp_macro-Amp_sim,     Buf_obs=Amp_macro-Amp_obs)]
  LONG[[v]] <- rbindlist(list(
    data.table(version=v, metric="Delta*T[min]~(degree*C)",        obs=m$dTmin_obs, sim=m$dTmin_sim),
    data.table(version=v, metric="Diurnal~amplitude~(degree*C)",   obs=m$Amp_obs,   sim=m$Amp_sim),
    data.table(version=v, metric="Buffering~(Amp[macro]-Amp[micro])", obs=m$Buf_obs, sim=m$Buf_sim)))
}
D <- rbindlist(LONG)
D[, version := factor(version, levels=names(VERS))]
D[, metric  := factor(metric, levels=c("Delta*T[min]~(degree*C)","Diurnal~amplitude~(degree*C)","Buffering~(Amp[macro]-Amp[micro])"))]
ann <- D[, .(lab=fitlab(obs,sim)), by=.(version,metric)]

# per-metric square limits via blank
rng <- D[, .(lo=min(c(obs,sim),na.rm=TRUE), hi=max(c(obs,sim),na.rm=TRUE)), by=metric]
p <- ggplot(D, aes(obs,sim)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50") +
  geom_hex(bins=60) + scale_fill_viridis_c(trans="log10", guide="none") +
  geom_text(data=ann, aes(x=-Inf,y=Inf,label=lab), parse=TRUE, hjust=-0.05, vjust=1.4, size=3.6, inherit.aes=FALSE) +
  facet_grid(metric ~ version, scales="free", labeller=labeller(metric=label_parsed)) +
  labs(x="Observed (HOBO)", y="Simulated (REF @1 m)",
       title="Daily ΔTmin / amplitude / buffering — sim vs HOBO, z05+×2",
       caption="ΔTmin: micro−macro min. Amplitude: daily Tmax−Tmin micro. Buffering: macro−micro diurnal range (>0 = forest damps).") +
  theme_bw(base_size=14) +
  theme(strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
        panel.grid.minor=element_blank(), plot.caption=element_text(size=10,colour="grey40"))
ggsave(file.path(OUT,"fig_sim_vs_obs_tmin_amp_buff.png"), p, width=11, height=13, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_sim_vs_obs_tmin_amp_buff.pdf"), p, width=11, height=13, device=cairo_pdf)
cli::cli_alert_success("Saved fig_sim_vs_obs_tmin_amp_buff")
print(ann)
