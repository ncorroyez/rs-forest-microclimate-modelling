# ==============================================================================
# A8 — Composite "glass ceiling" : (a) résidu ~ topo/hétéro + (b) pente buff/amp.
# Assemble les deux panneaux (objets ggplot sauvés par make_residual_vs_topo.R
# et make_slope_sim_vs_obs_groups.R) en une seule figure multi-panneaux.
#   Prérequis : avoir lancé les deux scripts (génèrent les .rds).
# Sortie : outputs/figures_pipeline/annex/fig_A8_glass_ceiling.{png,pdf}
# ==============================================================================
suppressMessages({ library(here); library(ggplot2); library(patchwork) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex")

pA <- readRDS(file.path(OUT, "_rds_A8a_topo.rds"))     # résidu vs topo (1 panneau)
pB <- readRDS(file.path(OUT, "_rds_A8b_slope.rds"))    # pente buff/amp/All (3 facettes)

comp <- (pA / pB) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  title = "Glass ceiling: (a) residual explained by topography/heterogeneity; (b) slope compression",
                  theme = theme(plot.title = element_text(size = 12, face = "bold")))
ggsave_article(file.path(OUT, "fig_A8_glass_ceiling"), comp, 12, 10)
cli::cli_alert_success("Saved fig_A8_glass_ceiling (composite a+b)")
