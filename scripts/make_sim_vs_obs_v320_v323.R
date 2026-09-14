# ==============================================================================
# Sim vs HOBO obs scatters, REF (1111) on z05+×2 inputs, for BOTH MuSICA versions
# (v3.2.0 legacy vs v3.2.3). Two scales:
#   (A) daily ΔTmax (micro−macro offset)  — from the 04-extract caches
#   (B) all hourly air temperature at 1 m (fixed 1 m interp) — sim Tmicro vs HOBO Tmicro
# Output: outputs/figures_pipeline/annex/fig_sim_vs_obs_v320_v323_1m.{png,pdf}
# ==============================================================================

suppressMessages({ library(data.table); library(here); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools); library(patchwork) })
source(here::here("R/config.R")); source(here::here("R/io.R")); source(here::here("R/musica.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

VERS <- list("v3.2.0 (legacy)" = "out_files/musica_hobo_z05/1111",
             "v3.2.3"          = "out_files/musica_hobo_z05_v323/1111")

# ---- HOBO hourly obs (1 m sensor) ------------------------------------------
hobo <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) %>%
  filter(position_sensor=="a", as.Date(datetime) %in% CFG$date_seq, !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime,"hour")) %>%
  group_by(id_plot, time) %>% summarise(t_obs = mean(t_hobo, na.rm=TRUE), .groups="drop") %>% as.data.table()

# ---- sim hourly Tair at FIXED 1 m (interpolated, HOBO height) per plot ------
sim_hourly <- function(dir) {
  rows <- list()
  for (f in list.files(here::here(dir), "\\.nc$", full.names=TRUE)) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$","\\1",basename(f))
    nc <- tryCatch(nc_open(f), error=function(e) NULL); if (is.null(nc)) next
    res <- tryCatch(get_tair_at_z(nc, z_target = 1.0), error=function(e) NULL); nc_close(nc)
    if (is.null(res)||nrow(res)==0) next
    s <- as.data.table(res)[, .(time=floor_date(time,"hour"), t_sim=Tair_sim)]   # Tair_sim already °C @1 m
    s$id_plot <- id; rows[[id]] <- s
  }
  rbindlist(rows)
}

# ---- assemble daily (from caches) + hourly (extract) ------------------------
fit_lab <- function(x,y){ ok<-is.finite(x)&is.finite(y); x<-x[ok];y<-y[ok]
  sprintf("italic(r)==%.2f~~RMSE==%.2f~~bias==%+.2f~~n==%d", cor(x,y), sqrt(mean((y-x)^2)), mean(y-x), length(x)) }

daily <- list(); hourly <- list()
for (v in names(VERS)) {
  br <- if (v=="v3.2.0 (legacy)") "z05" else "z05_v323"
  V <- readRDS(here::here(sprintf("outputs/figures_pipeline_%s/data/ref_validation.rds", br)))
  d <- merge(as.data.table(V$ref_daily), as.data.table(V$obs_daily), by=c("id_plot","date"))
  daily[[v]] <- d[is.finite(Delta_sim)&is.finite(Delta_obs), .(version=v, obs=Delta_obs, sim=Delta_sim)]
  h <- merge(sim_hourly(VERS[[v]]), hobo, by=c("id_plot","time"))
  hourly[[v]] <- h[is.finite(t_sim)&is.finite(t_obs), .(version=v, obs=t_obs, sim=t_sim)]
}
DD <- rbindlist(daily); HH <- rbindlist(hourly)
DD[, version := factor(version, levels=names(VERS))]; HH[, version := factor(version, levels=names(VERS))]
annD <- DD[, .(lab=fit_lab(obs,sim)), by=version]; annH <- HH[, .(lab=fit_lab(obs,sim)), by=version]

mk <- function(D, ann, xlab, ylab, ttl, hexfill){
  lim <- range(c(D$obs,D$sim));
  ggplot(D, aes(obs,sim)) + geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50") +
    geom_hex(bins=70) + scale_fill_viridis_c(trans="log10", option=hexfill, guide="none") +
    geom_text(data=ann, aes(x=-Inf,y=Inf,label=lab), parse=TRUE, hjust=-0.05, vjust=1.4, size=4, inherit.aes=FALSE) +
    facet_wrap(~version, nrow=1) + coord_equal(xlim=lim,ylim=lim) +
    labs(x=xlab, y=ylab, title=ttl) + theme_bw(base_size=15) +
    theme(strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"), panel.grid.minor=element_blank())
}
pA <- mk(DD, annD, expression("Observed daily "*Delta*T[max]~"("*degree*"C)"),
         expression("Simulated daily "*Delta*T[max]~"("*degree*"C)"), "Daily ΔTmax — REF vs 53 HOBO (z05+×2)", "C")
pB <- mk(HH, annH, expression("Observed hourly "*T[air]^{1*m}~"("*degree*"C)"),
         expression("Simulated hourly "*T[air]^{1*m}~"("*degree*"C)"), "All hourly air temperature @1 m (interp) — REF vs HOBO", "D")
p <- pA / pB
ggsave(file.path(OUT,"fig_sim_vs_obs_v320_v323_1m.png"), p, width=12, height=11, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_sim_vs_obs_v320_v323_1m.pdf"), p, width=12, height=11, device=cairo_pdf)
cli::cli_alert_success("Saved fig_sim_vs_obs_v320_v323_1m")
print(rbind(cbind(scale="daily ΔTmax", annD), cbind(scale="hourly Tair", annH)))
