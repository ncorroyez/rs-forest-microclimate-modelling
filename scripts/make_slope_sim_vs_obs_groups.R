# ==============================================================================
# Pente micro–macro : MuSICA (y) vs HOBO (x) par placette, séparée ET regroupée.
# Facettes {Buffering, Amplifying, All} — groupes définis sur le TERRAIN
# (slope_obs < 1 = buffering ; > 1 = amplifying ; le split 45/8).
# Montre que le modèle suit le classement des tamponnantes mais plafonne (~1.1)
# sur les amplificatrices que le terrain pousse jusqu'à ~1.4.
# Branche z05 (v3.2.0, article).
#
# ⚠ À LIRE AVEC make_residual_vs_topo.R (A10) : les 8 amplificatrices = les
# placettes à trouées/pente/exposition, dont l'amplification est topo/gap-driven
# (PAS un input MuSICA). Le "sous-amplifie" est donc une limite de PÉRIMÈTRE, pas
# un défaut structurel. Le seul écart imputable au modèle = la compression du
# buffering sur les placettes fermées (gap≈0). Ne pas présenter A11 seule comme
# une preuve indépendante d'échec du modèle.
#
# Sortie : outputs/figures_pipeline/annex/fig_slope_sim_vs_obs_groups.{png,pdf}
# ==============================================================================

suppressMessages({ library(here); library(data.table); library(tidyverse) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

V <- readRDS(here::here("outputs/figures_pipeline_z05/data/ref_validation.rds"))
D <- merge(as.data.table(V$ref_slope), as.data.table(V$obs_slope), by = "id_plot")
D[, grp := ifelse(slope_obs > 1, "Amplifying (obs>1)", "Buffering (obs<1)")]

# "séparé ET regroupé" : on duplique tout dans une facette All
DD <- rbind(D[, .(slope_obs, slope_sim, facet = grp)],
            D[, .(slope_obs, slope_sim, facet = "All (n=53)")])
DD[, facet := factor(facet, levels = c("Buffering (obs<1)", "Amplifying (obs>1)", "All (n=53)"))]

ann <- DD[, .(lab = sprintf("r=%.2f  n=%d", cor(slope_obs, slope_sim), .N)), by = facet]
lim <- range(c(DD$slope_obs, DD$slope_sim), na.rm = TRUE)

p <- ggplot(DD, aes(slope_obs, slope_sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +     # 1:1
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey70") +
  geom_point(aes(colour = slope_obs > 1), size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = c(`FALSE` = unname(PAL_GRP["buffering"]),
                                 `TRUE`  = unname(PAL_GRP["amplifying"])), guide = "none") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.4,
            size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~ facet, nrow = 1) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "HOBO micro–macro slope (per plot)", y = "MuSICA micro–macro slope (per plot)") +
  theme_article() + theme(aspect.ratio = 1)
ggsave_article(file.path(OUT, "fig_slope_sim_vs_obs_groups"), p, 12, 4.6)
saveRDS(p, file.path(OUT, "_rds_A8b_slope.rds"))      # pour le composite A8

cat("\n=== pente modèle vs HOBO par groupe ===\n")
print(D[, .(n = .N, slope_obs_med = round(median(slope_obs),3), slope_sim_med = round(median(slope_sim),3),
            r = round(cor(slope_obs, slope_sim),2)), by = grp])
cat(sprintf("All : r = %.2f (n=%d)\n", cor(D$slope_obs, D$slope_sim), nrow(D)))
cli::cli_alert_success("Saved fig_slope_sim_vs_obs_groups")
