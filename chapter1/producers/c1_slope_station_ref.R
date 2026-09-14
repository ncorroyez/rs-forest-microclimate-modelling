# ==============================================================================
# Réu 2026-06-26 #6/#7 — Pentes de tamponnement MESURÉES référencées à la station
# CHS 41 (1.5 m), boxplot par cluster. Méthode canonique du chapitre (04_extract.R:73):
#   β = coef(lm(t_hobo_horaire ~ Tair_macro_horaire))  par plot.
# Ici macro = CHS 41 (1.5 m) au lieu d'ERA5. Contrairement à ΔTmax (offset pur),
# β dépend de la covariance horaire => la station change réellement β.
# β<1 = tamponnement (le sous-bois suit moins l'open-field).
#   Rscript c1_slope_station_ref.R
# Out: tab_slope_station_ref.csv + FigStation_slope_boxplot.png
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(data.table); library(lubridate); library(dplyr); library(sf)
  library(ggplot2); library(patchwork)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("pipeline/00_config.R")
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)
DS  <- CFG$date_seq

# ---- macro horaire : ERA5 (forçage) et CHS 41 (1.5 m) -------------------------
era5_h <- as.data.table(build_era5_hourly(CFG$forcing_file, DS))            # time, Tair_era5
raw <- fread("MetHor2021.txt", sep=";", header=TRUE, encoding="Latin-1"); setnames(raw,1,"code")
tcol <- grep("instantan", names(raw), value=TRUE)[1]
chs <- raw[code=="CHS 41", .(Date, H=`Heure (TU)`, Tair_chs=as.numeric(get(tcol)))]
chs[, time := as.POSIXct(Date, format="%d/%m/%Y", tz="UTC") + lubridate::hours(H %/% 100)]
chs_h <- chs[as.Date(time) %in% DS & is.finite(Tair_chs), .(time, Tair_chs)]

# ---- HOBO horaire (mirror 04_extract.R:68-72) --------------------------------
hobo <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime=as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")) %>%
  filter(position_sensor=="a", as.Date(datetime) %in% DS, !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time=floor_date(datetime,"hour")) %>%
  group_by(id_plot, time) %>% summarise(t_hobo=mean(t_hobo,na.rm=TRUE), .groups="drop") %>% as.data.table()

# ---- common-hour merge, β per plot vs each macro -----------------------------
M <- Reduce(function(a,b) merge(a,b,by="time"), list(hobo, era5_h, chs_h))
slp <- M[, .(b_era5 = coef(lm(t_hobo ~ Tair_era5))[2],
             b_chs  = coef(lm(t_hobo ~ Tair_chs))[2], n_hr=.N), by=id_plot]
cl <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, P)]
P  <- merge(slp, cl, by="id_plot"); P[, P := factor(P, levels=names(PAL))]
fwrite(P, "out_files/Chapter1/tables/tab_slope_station_ref.csv")

# ---- sanity : β vs ERA5 doit retomber sur V$obs_slope ------------------------
V <- readRDS("outputs/figures_pipeline_z05/data/ref_validation.rds")
chk <- merge(P[, .(id_plot, b_era5)], as.data.table(V$obs_slope)[, .(id_plot, slope_obs)], by="id_plot")
cat(sprintf("sanity β(ERA5) vs V$obs_slope : r=%.4f  max|Δ|=%.4f  (n_hr médian=%d)\n",
            cor(chk$b_era5, chk$slope_obs), max(abs(chk$b_era5-chk$slope_obs)), median(P$n_hr)))
cat("\n== β médian par cluster ==\n")
print(P[, .(med_era5=round(median(b_era5),3), med_chs=round(median(b_chs),3), n=.N), by=P][order(P)])
cat(sprintf("\nΔβ obs (CHS41 − ERA5) : médiane %+.3f  (range %.3f..%.3f)\n",
            median(P$b_chs-P$b_era5), min(P$b_chs-P$b_era5), max(P$b_chs-P$b_era5)))

# ---- Figure : (a) boxplot β mesurées par cluster ; (b) effet réf -------------
pA <- ggplot(P, aes(P, b_chs, colour=P, fill=P)) +
  geom_hline(yintercept=1, linetype=2, colour="grey45") +
  geom_boxplot(width=0.55, alpha=0.18, outlier.shape=NA) +
  geom_jitter(width=0.12, height=0, size=2, alpha=0.85) +
  scale_colour_manual(values=PAL, guide="none") + scale_fill_manual(values=PAL, guide="none") +
  labs(x=NULL, y=expression(beta~"= pente "*T[micro]*" ~ "*T[CHS41]),
       title="(a) Pentes de tamponnement mesurées par cluster (réf. station 1,5 m)",
       subtitle="β<1 = tamponnement ; le gradient ouvert→dense (P1→P4) se déploie sur les obs") +
  theme_bw(base_size=11) + theme(plot.title=element_text(size=10,face="bold"),
       plot.subtitle=element_text(size=8,colour="grey35"))

rng <- range(c(P$b_era5, P$b_chs))
pB <- ggplot(P, aes(b_era5, b_chs, colour=P)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey50") +
  geom_point(size=2.2, alpha=0.85) +
  scale_colour_manual(values=PAL, name=NULL) + coord_equal(xlim=rng, ylim=rng) +
  labs(x="β (réf. ERA5 forçage)", y="β (réf. station CHS 41)",
       title="(b) Pente quasi-inchangée par la référence",
       subtitle="points ~ sur la 1:1 (Δβ médian −0,003) : ERA5 et CHS 41 très corrélés à l'heure") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))

ggsave("out_files/Chapter1/figures/FigStation_slope_boxplot.png",
       pA + pB + plot_layout(widths=c(1.05, 1)), width=11, height=5.2, dpi=200, bg="white")
cat("DONE -> tab_slope_station_ref.csv + FigStation_slope_boxplot.png\n")
