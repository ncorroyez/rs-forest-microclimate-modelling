# ==============================================================================
# Réu 2026-06-26 #1/#4 — ΔTmax référencé à la station open-field CHS 41 (1.5 m)
# au lieu du forçage free-air ERA5. HOBO à 1 m, station à 1.5 m (écart 0.5 m).
#
# Référence spatialement uniforme + agrégation par moyenne => le changement de
# référence est un OFFSET uniforme exact :
#   ΔTmax_CHS41(plot) = ΔTmax_ERA5(plot) + mean(ERA5 - CHS41)
# => r, ρ, biais modèle, SD, classement INVARIANTS ; seuls bougent le NIVEAU
#    absolu (+offset) et le passage à zéro (compte amplificateurs/tampons).
# Résultat : l'amplification ne disparaît pas, elle grossit (station plus froide
# qu'ERA5). #4 "si pas d'amp" ne vient donc PAS du choix de référence.
#   Rscript c1_deltatmax_station_ref.R
# Out: tab_deltatmax_station_ref.csv + FigStation_deltatmax_ref.png
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(data.table); library(lubridate); library(dplyr)
  library(ggplot2); library(patchwork)
})
source("R/config.R"); source("R/io.R")
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)
WIN <- CFG$date_seq

# ---- δ(date) = ERA5_max − CHS41_max (1.5 m, ≥20 h/jour) -----------------------
era5 <- as.data.table(extract_macro_daily(CFG$forcing_file, WIN))[, .(date, ERA5 = Tmax_macro)]
raw  <- fread("MetHor2021.txt", sep=";", header=TRUE, encoding="Latin-1")
setnames(raw, 1, "code"); tcol <- grep("instantan", names(raw), value=TRUE)[1]
chs  <- raw[code=="CHS 41", .(Date, Tair=as.numeric(get(tcol)))]
chs[, date := as.Date(Date, format="%d/%m/%Y")]
chs_d <- chs[date %in% WIN & is.finite(Tair), .(CHS41=max(Tair), n_hr=.N), by=date][n_hr>=20]
delta <- merge(era5, chs_d, by="date")[, .(date, d = ERA5 - CHS41)]
OFF <- mean(delta$d)
cat(sprintf("δ=ERA5−CHS41 (max j): %+.3f °C (sd %.3f, n=%d j)\n", OFF, sd(delta$d), nrow(delta)))

# ---- re-référencement par plot (obs + sim) -----------------------------------
V <- readRDS("outputs/figures_pipeline_z05/data/ref_validation.rds")
obs <- merge(as.data.table(V$obs_daily), delta, by="date")[
  , .(obs_era5=mean(Delta_obs), obs_chs=mean(Delta_obs + d)), by=id_plot]
sim <- merge(as.data.table(V$ref_daily), delta, by="date")[
  , .(sim_era5=mean(Delta_sim), sim_chs=mean(Delta_sim + d)), by=id_plot]
cl <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, P, LAI, VCI)]
P  <- Reduce(function(a,b) merge(a,b,by="id_plot"), list(obs, sim, cl))
P[, P := factor(P, levels=names(PAL))]
P[, amp_obs := obs_chs > 0][, amp_era5 := obs_era5 > 0]

fwrite(P, "out_files/Chapter1/tables/tab_deltatmax_station_ref.csv")

# ---- métriques (invariance) ---------------------------------------------------
met <- function(s,o) sprintf("r=%.3f  ρ=%.3f\nbiais=%+.2f °C\nRMSE=%.2f\nn=%d",
  cor(s,o), cor(s,o,method="spearman"), mean(s-o), sqrt(mean((s-o)^2)), length(o))
cat("\n== validation CHS41 ref ==\n", met(P$sim_chs,P$obs_chs), "\n")
cat(sprintf("\namplificateurs OBS : ERA5 %d -> CHS41 %d  (sur %d)\n",
            sum(P$amp_era5), sum(P$amp_obs), nrow(P)))
cat("amplificateurs CHS41 ref par cluster :\n"); print(P[amp_obs==TRUE, .N, by=P][order(P)])

# ---- Figure ------------------------------------------------------------------
rng <- range(c(P$obs_chs, P$sim_chs))
pA <- ggplot(P, aes(obs_chs, sim_chs, colour=P)) +
  geom_hline(yintercept=0, linetype=3, colour="grey60") +
  geom_vline(xintercept=0, linetype=3, colour="grey60") +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey45") +
  geom_point(size=2.4, alpha=0.85) +
  annotate("text", x=rng[1], y=rng[2], hjust=0, vjust=1, size=3, colour="grey25",
           label=met(P$sim_chs,P$obs_chs)) +
  scale_colour_manual(values=PAL, name=NULL) + coord_equal(xlim=rng, ylim=rng) +
  labs(x="ΔTmax observé (HOBO 1 m − CHS 41 1,5 m, °C)",
       y="ΔTmax simulé (MuSICA 1 m − CHS 41, °C)",
       title="(a) Validation référencée à la station open-field") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(size=10,face="bold"))

pB <- ggplot(P, aes(P, obs_chs, colour=P, fill=P)) +
  geom_hline(yintercept=0, linetype=2, colour="grey45") +
  geom_boxplot(width=0.5, alpha=0.18, outlier.shape=NA) +
  geom_jitter(width=0.12, height=0, size=2, alpha=0.85) +
  scale_colour_manual(values=PAL, guide="none") + scale_fill_manual(values=PAL, guide="none") +
  labs(x=NULL, y="ΔTmax observé (HOBO − CHS 41, °C)",
       title="(b) Amplification (ΔTmax>0) concentrée dans les plots ouverts (P1–P2)",
       subtitle=sprintf("réf. station = +%.2f °C vs ERA5 ; amplificateurs 9→%d (rang & SD inchangés)", OFF, sum(P$amp_obs))) +
  theme_bw(base_size=11) + theme(plot.title=element_text(size=10,face="bold"),
       plot.subtitle=element_text(size=8,colour="grey35"))

ggsave("out_files/Chapter1/figures/FigStation_deltatmax_ref.png",
       pA + pB + plot_layout(widths=c(1, 0.95)), width=11, height=5.2, dpi=200, bg="white")
cat("DONE -> tab_deltatmax_station_ref.csv + FigStation_deltatmax_ref.png\n")
