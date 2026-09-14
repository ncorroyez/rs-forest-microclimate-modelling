# ==============================================================================
# R2ter — Validation EMPIRIQUE de H2 sur les 53 placettes (données réelles)
#
# H2 : à LAI comparable, les profils dont la biomasse foliaire est plus HAUTE
# (centre de masse élevé) tamponnent davantage le ΔTmax. Ici, AUCUN profil
# synthétique : on utilise les VRAIS profils LAD LiDAR (ALS z0.5) par placette,
# le buffering OBSERVÉ (HOBO) et les clusters P1–P4.
#
#   buffering observé  = ΔTmax_obs moyen (micro − macro), plus négatif = + tamponné
#   top-heaviness      = centre de masse du profil LAD réel (fraction de Hmax)
#   test               = effet partiel de COM sur le buffering, À LAI CONTRÔLÉ
#
# Sortie : outputs/figures_pipeline/annex/fig_h2_field_validation.{png,pdf}
#          + tab_h2_field_validation.csv
# ==============================================================================

suppressMessages({ library(here); library(data.table); library(tidyverse); library(patchwork) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 1. centre de masse du profil LAD RÉEL par placette ----------------------
lad <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))
ladc <- grep("^LAD_Layer_", names(lad), value = TRUE)
hts  <- as.numeric(gsub("LAD_Layer_", "", ladc))            # hauteurs des couches (m)
com_of <- function(row) {                                   # centre de masse / Hmax
  d <- as.numeric(row[ladc]); d[is.na(d)] <- 0
  if (sum(d) <= 0) return(NA_real_)
  (sum(d * hts) / sum(d)) / as.numeric(row[["Hmax"]])
}
lad[, COM := apply(.SD, 1, com_of), .SDcols = c(ladc, "Hmax")]
plots <- lad[, .(id_plot, LAI, Hmax, COM)]

# ---- 2. buffering OBSERVÉ (HOBO) + cluster par placette ----------------------
V   <- readRDS(here::here("outputs/figures_pipeline_z05/data/ref_validation.rds"))
obs <- as.data.table(V$obs_daily)[, .(dTmax_obs = mean(Delta_obs, na.rm = TRUE)), by = id_plot]
slp <- as.data.table(V$obs_slope)[, .(id_plot, slope_obs)]   # pente micro–macro
clu <- unique(as.data.table(readRDS(here::here("outputs/figures_pipeline_z05/data/shapley_perplot.rds")))[, .(id_plot, Cluster)])
fc  <- fread(here::here("transfer/musica_version_benchmark/run_simulations/inputs/fcover_per_plot.csv"))

D <- Reduce(function(a, b) merge(a, b, by = "id_plot"), list(plots, obs, slp, clu, fc))
D <- D[is.finite(COM) & is.finite(dTmax_obs) & is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
fwrite(D, file.path(OUT, "tab_h2_field_validation.csv"))

# ---- 3. effet partiel de COM, à TOUS les autres traits contrôlés -------------
# Comparaison juste avec φ_LAD (Shapley contrôle LAI+Hmax+fCover). Deux niveaux :
#   (a) contrôle LAI seul ;  (b) contrôle LAI+Hmax+fCover (= conditions du Shapley).
part_r <- function(ctrl) {                       # corr. partielle COM~buffering | ctrl
  fx <- as.formula(paste("COM ~", paste(ctrl, collapse = "+")))
  fy <- as.formula(paste("dTmax_obs ~", paste(ctrl, collapse = "+")))
  rx <- resid(lm(fx, data = D)); ry <- resid(lm(fy, data = D))
  c(list(rx = rx, ry = ry), cor.test(rx, ry)[c("estimate","p.value")])
}
pa <- part_r("LAI")
pb <- part_r(c("LAI","Hmax","fCover"))
fitFull <- summary(lm(dTmax_obs ~ LAI + Hmax + fCover + COM, data = D))$coefficients

cat(sprintf("\nN=%d | colinéarité COM~LAI r=%.2f | COM~Hmax r=%.2f | COM~fCover r=%.2f\n",
            nrow(D), cor(D$COM,D$LAI), cor(D$COM,D$Hmax), cor(D$COM,D$fCover)))
cat(sprintf("(a) COM | LAi seul          : r_part=%+.2f  p=%.3f\n", pa$estimate, pa$p.value))
cat(sprintf("(b) COM | LAI+Hmax+fCover   : r_part=%+.2f  p=%.3f  <-- comparable au Shapley φ_LAD\n", pb$estimate, pb$p.value))
cat(sprintf("    coef(COM) régression complète = %+.3f (p=%.3f)\n", fitFull["COM","Estimate"], fitFull["COM","Pr(>|t|)"]))
rx <- pb$rx; ry <- pb$ry; ct <- list(estimate = pb$estimate, p.value = pb$p.value)  # figure = contrôle complet

# ---- 4. figures : (A) brut par cluster  (B) added-variable plot (effet partiel) --
labA <- sprintf("raw: r=%.2f", cor(D$COM, D$dTmax_obs))

pA <- ggplot(D, aes(COM, dTmax_obs, colour = Cluster)) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
  labs(x = "Real LAD profile centre of mass  (fraction of Hmax)",
       y = expression(bar(Delta*T[max])~obs~"HOBO  ("*degree*"C)"),
       subtitle = sprintf("Raw (clusters P1–P4)  |  %s", labA)) +
  theme_article() + legend_corner()

avp <- data.table(rx = rx, ry = ry, Cluster = D$Cluster)
labB <- sprintf("partial | LAI+Hmax+fCover: r=%.2f  (p=%.3f)", ct$estimate, ct$p.value)
pB <- ggplot(avp, aes(rx, ry, colour = Cluster)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
  labs(x = "COM  |  LAI+Hmax+fCover  (residuals)", y = expression(bar(Delta*T[max])~obs~"| traits (residuals)"),
       subtitle = sprintf("Partial effect of vertical position, ALL traits controlled  |  %s", labB)) +
  theme_article() + theme(legend.position = "none")

p <- patchwork::wrap_plots(pA, pB, ncol = 2)
ggsave_article(file.path(OUT, "fig_h2_field_validation"), p, 12, 5.2)
cli::cli_alert_success("Saved fig_h2_field_validation (png + pdf) + table")
