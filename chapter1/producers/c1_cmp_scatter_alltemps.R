# ==============================================================================
# Scatter MuSICA vs HOBO sur TOUTES les températures horaires à 1 m (été 2021,
# 53 loggers × toutes les heures), pour v3.2.0 (legacy) vs v3.2.3 iter (yoyo, real BLH).
# Complément de la validation par-plot (ΔTmax) : montre l'ajustement horaire complet.
#   Rscript c1_cmp_scatter_alltemps.R
# Out: FigCmp_scatter_hourly_v320_v323iter.png (+ copie dans comparaison_versions/figures)
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("pipeline/00_config.R")

# --- HOBO horaire (capteur air, 1 m) ---
hobo <- as.data.table(read.csv(CFG$hobo_temp_csv))
hobo <- hobo[position_sensor=="a" & !(id_plot %in% CFG$ids_to_remove)]
hobo[, datetime := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hobo <- hobo[as.Date(datetime) %in% CFG$date_seq]
hobo[, time := floor_date(datetime, "hour")]
hobo <- hobo[, .(t_obs = mean(t_hobo, na.rm=TRUE)), by=.(id_plot, time)]

sim_hourly <- function(dir) {
  rbindlist(lapply(list.files(dir, "\\.nc$", full.names=TRUE), function(f) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$","\\1", basename(f))
    nc <- tryCatch(nc_open(f), error=function(e) NULL); if (is.null(nc)) return(NULL)
    res <- tryCatch(get_tair_at_z(nc, z_target=1.0), error=function(e) NULL); nc_close(nc)
    if (is.null(res) || nrow(res)==0) return(NULL)
    s <- as.data.table(res)[, .(time=floor_date(time,"hour"), t_sim=Tair_sim)]; s$id_plot <- id; s
  }), fill=TRUE)
}
VERS <- c("v3.2.0"="out_files/musica_hobo_z05/1111", "v3.2.3 iter"="out_files/musica_hobo_z05_iter/1111")
D <- rbindlist(lapply(names(VERS), function(v) {
  m <- merge(sim_hourly(VERS[[v]]), hobo, by=c("id_plot","time"))
  m[is.finite(t_sim) & is.finite(t_obs), .(version=v, t_obs, t_sim)]
}))
lab <- D[, .(r2=cor(t_sim,t_obs)^2, rmse=sqrt(mean((t_sim-t_obs)^2)), bias=mean(t_sim-t_obs), n=.N), by=version]
lab[, txt := sprintf("R² = %.2f\nRMSE = %.2f °C\nbias = %+.2f °C\nn = %s h", r2, rmse, bias, format(n, big.mark=" "))]
cat("=== validation horaire (toutes températures 1 m) ===\n"); print(lab[,.(version,r2=round(r2,3),rmse=round(rmse,2),bias=round(bias,2),n)])
D[, version := factor(version, levels=names(VERS))]; lab[, version := factor(version, levels=names(VERS))]
rng <- range(c(D$t_obs, D$t_sim), na.rm=TRUE)
p <- ggplot(D, aes(t_obs, t_sim)) +
  geom_hex(bins=70) + scale_fill_viridis_c(trans="log10", name="n (heures)", option="magma") +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey40") +
  geom_smooth(method="lm", se=FALSE, colour="#1A9850", linewidth=0.7) +
  geom_text(data=lab, aes(x=rng[1], y=rng[2], label=txt), hjust=0, vjust=1, size=3, colour="grey15") +
  facet_wrap(~version, nrow=1) + coord_equal(xlim=rng, ylim=rng) +
  labs(x="T air observée (HOBO, 1 m, °C)", y="T air simulée (MuSICA, 1 m, °C)",
       title="MuSICA vs HOBO — toutes les températures horaires (été 2021, 53 loggers)",
       subtitle="v3.2.0 (legacy validé) vs v3.2.3 iter (yoyo, vraie BLH). 1:1 tireté, ajustement vert.") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11,face="bold"), plot.subtitle=element_text(size=8.3,colour="grey35"))
out <- "out_files/Chapter1/figures/FigCmp_scatter_hourly_v320_v323iter.png"
ggsave(out, p, width=9.5, height=5.2, dpi=200, bg="white")
file.copy(out, "Chapitre1/comparaison_versions/figures/", overwrite=TRUE)
cat("DONE ->", out, "(+ copie comparaison_versions/figures)\n")
