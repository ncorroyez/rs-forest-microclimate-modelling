# ==============================================================================
# PIPELINE STAGE 04 — Extract per-plot × per-coalition metrics from the HOBO set.
# Two DISTINCT conventions (constraint #4), never unified:
#   • absolute ΔTmax / ΔVPDmax : fixed 1 m interp + legacy -2 h shift (HOBO-comparable)
#   • relative micro/macro slope : ALSO fixed 1 m, no shift (NOT nair==1; legacy 45-8 superseded)
# Guard (constraint #3): NEVER read DT_daily_HOBO_floor05_* — that cache carries a
# THIRD (validation) baseline, not mean.
# Outputs (OUT_DATA): coal_metrics.rds, ref_validation.rds, clusters.rds
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
suppressMessages({ library(terra); library(musica.tools) })
cli_h1("STAGE 04 — extract metrics (all at fixed 1 m ; slope no shift)")

era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)
df_macro    <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq))
hot_days    <- df_macro[Tmax_macro >= quantile(Tmax_macro, PIPE$HOT_QUANTILE, na.rm=TRUE), date]
cli_alert("Hottest 10% : {length(hot_days)} days (macro Tmax >= {round(quantile(df_macro$Tmax_macro,PIPE$HOT_QUANTILE),1)}°C)")

hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
ids <- hobo_pts$id_plot

# ---- daily Tmax & VPDmax at fixed 1 m (interp) + -2 h shift ------------------
daily_1m <- function(path) {
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec <- t0 + dhours(ncvar_get(nc, "time")) - lubridate::hours(PIPE$TMAX_SHIFT_HR)
  Tk <- ncvar_get(nc, "Tair_z"); wmr <- ncvar_get(nc, "wair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh * vh; Z <- PIPE$Z_FIX
  if (Z <= zl[1]) { ilo<-1L; ihi<-1L; w<-0 } else if (Z >= zl[length(zl)]) { ilo<-length(zl); ihi<-ilo; w<-0 } else {
    ilo <- max(which(zl <= Z)); ihi <- ilo+1L; w <- (Z-zl[ilo])/(zl[ihi]-zl[ilo]) }
  Tc <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  wv <- (1-w)*wmr[ilo,] + w*wmr[ihi,]
  vpd <- pmax(esat_hpa(Tc) - (wv/(1+wv))*PIPE$P_HPA, 0)/10
  data.table(date = as.Date(floor_date(tvec,"hour")), Tc = Tc, vpd = vpd)[
    date %in% CFG$date_seq, .(Tmax = max(Tc), VPDmax = max(vpd)), by = date]
}
nc_of <- function(bit, id) file.path(PIPE$MUSICA_DIR_HOBO, bit, sprintf("musica_out_HOBO_%s.nc", id))

# ---- per-coalition metrics (for Shapley) ------------------------------------
rows <- list()
for (bit in PIPE$BITS) for (id in ids) {
  d <- daily_1m(nc_of(bit, id)); if (is.null(d) || nrow(d)==0) next
  rows[[length(rows)+1L]] <- data.table(id_plot=id, bit=bit,
    Tmax_all=mean(d$Tmax), Tmax_hot=mean(d[date %in% hot_days, Tmax]),
    VPD_all=mean(d$VPDmax),  VPD_hot=mean(d[date %in% hot_days, VPDmax]))
}
coal <- rbindlist(rows)
pipe_assert(coal[, uniqueN(bit), by=id_plot][, all(V1==16)], "All plots have 16 coalitions")
saveRDS(coal, file.path(PIPE$OUT_DATA, "coal_metrics.rds"))

# ---- REF (1111) validation metrics ------------------------------------------
# (a) per-plot micro/macro log slope at fixed 1 m, no shift (vs HOBO)
ref_slope <- rbindlist(lapply(ids, function(id) {
  r <- tryCatch(extract_hourly_slope_one(nc_of(PIPE$REF_BIT, id), era5_hourly, CFG$date_seq,
                 use_nair1 = PIPE$SLOPE_USE_NAIR1, time_shift_hr = PIPE$SLOPE_SHIFT_HR), error=function(e) NULL)
  if (is.null(r) || nrow(r)==0) return(NULL)
  data.table(id_plot=id, slope_sim=r$slope, log_slope_sim=log(r$slope)) }))
# (b) per-plot daily ΔTmax (1 m, -2 h) for sim-vs-obs r
ref_daily <- rbindlist(lapply(ids, function(id) {
  d <- daily_1m(nc_of(PIPE$REF_BIT, id)); if (is.null(d)) return(NULL)
  m <- merge(d, df_macro, by="date"); m[, .(id_plot=id, date, Delta_sim = Tmax - Tmax_macro)] }))

# ---- HOBO observations (daily Δobs + per-plot slope) ------------------------
hobo_raw <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) %>%
  filter(position_sensor=="a", as.Date(datetime) %in% CFG$date_seq, !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime,"hour")) %>%
  group_by(id_plot, time) %>% summarise(t_hobo=mean(t_hobo, na.rm=TRUE), .groups="drop") %>% as.data.table()
obs_slope <- merge(hobo_raw, era5_hourly, by="time")[, {
  f<-lm(t_hobo~Tair_era5); .(slope_obs=as.numeric(coef(f)[2]), log_slope_obs=log(as.numeric(coef(f)[2]))) }, by=id_plot]
obs_daily <- hobo_raw[, .(id_plot, date=as.Date(time), t_hobo)][, .(Tmax_obs=max(t_hobo)), by=.(id_plot,date)]
obs_daily <- merge(obs_daily, df_macro, by="date")[, .(id_plot, date, Delta_obs = Tmax_obs - Tmax_macro)]

# ---- clusters (FROZEN assignment, constraint #1) ----------------------------
df_forest <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))[, .(x,y,Cluster=as.character(Cluster))]
hxy <- sf::st_coordinates(hobo_pts)
clu <- data.table(id_plot=ids, x=hxy[,"X"], y=hxy[,"Y"])
clu[, Cluster := relabel_cluster(sapply(seq_len(.N), function(i)
  df_forest$Cluster[which.min((df_forest$x-x[i])^2+(df_forest$y-y[i])^2)]))]
clu[, Cluster := factor(Cluster, levels=paste0("P",1:4))]

saveRDS(list(ref_slope=ref_slope, ref_daily=ref_daily, obs_slope=obs_slope,
             obs_daily=obs_daily, hot_days=hot_days),
        file.path(PIPE$OUT_DATA, "ref_validation.rds"))
saveRDS(clu[, .(id_plot, Cluster)], file.path(PIPE$OUT_DATA, "clusters.rds"))
cli_alert_success("Stage 04 done — coal_metrics, ref_validation, clusters cached.")
