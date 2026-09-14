# ==============================================================================
# Pictogrammes des 4 archétypes (P1..P4) — petites silhouettes d'arbre
# encodant LAI (densité du feuillage) + Hmax (hauteur).
# Sortie : outputs/figs_MEB2026_final/pictogram_P{1..4}.png (transparent)
#         + outputs/figs_MEB2026_final/pictograms_strip.png (les 4 côte à côte)
# ==============================================================================

suppressMessages({
  library(here); library(ggplot2); library(grid); library(data.table)
})

OUT <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(here::here("R/cluster_relabel.R"))

# Archetype values from V9 pipeline reordered by ASCENDING LAI
# (sparsest P1 → densest P4). Old cluster IDs in parentheses.
arch <- data.table(
  Cluster = c("P1", "P2", "P3", "P4"),
  # P1 = old C3 sparse ; P2 = old C2 ; P3 = old C4 ; P4 = old C1 dense
  LAI     = c(2.35,  6.05,  7.72, 10.80),
  Hmax    = c(20.67, 14.80, 28.46, 33.70),
  fCover  = c(0.55,  0.96,  0.96,  0.99),
  colour  = unname(PAL_CLUSTER[c("P1","P2","P3","P4")])
)

# Normalisation pour le dessin (height/canopy width)
arch[, h_norm     := Hmax / max(Hmax)]                # 0..1
arch[, crown_w    := pmin(1.0, 0.35 + 0.6 * (LAI / max(LAI)))]
arch[, crown_alpha := 0.25 + 0.7 * (LAI / max(LAI))]
arch[, trunk_h    := 0.20 + 0.10 * h_norm]
arch[, crown_centre := trunk_h + (h_norm - trunk_h) / 2]
arch[, crown_h    := (h_norm - trunk_h)]

draw_one <- function(i) {
  a <- arch[i]
  # Trunk : rectangle gris foncé
  trunk <- data.frame(
    x = c(-0.05, 0.05, 0.05, -0.05),
    y = c(0, 0, a$trunk_h, a$trunk_h))
  # Crown : ellipse approximée par polygon
  theta <- seq(0, 2*pi, length.out = 60)
  cx <- a$crown_w/2 * cos(theta)
  cy <- a$crown_centre + a$crown_h/2 * sin(theta)
  crown <- data.frame(x = cx, y = cy)
  ggplot() +
    # base ground line
    geom_segment(aes(x = -0.6, xend = 0.6, y = 0, yend = 0),
                 colour = "grey50", linewidth = 0.6) +
    # trunk
    geom_polygon(data = trunk, aes(x = x, y = y),
                 fill = "#5a3a2a", colour = "#3e2818", linewidth = 0.3) +
    # crown
    geom_polygon(data = crown, aes(x = x, y = y),
                 fill = a$colour, colour = "grey25",
                 alpha = a$crown_alpha, linewidth = 0.5) +
    # archetype label only (no metrics)
    annotate("text", x = 0, y = 1.10, label = a$Cluster,
              fontface = "bold", size = 9, colour = "grey15") +
    coord_fixed(xlim = c(-0.6, 0.6), ylim = c(-0.10, 1.18)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA))
}

# Single PNGs (transparent background, slide-friendly, 600 dpi)
for (i in 1:4) {
  p <- draw_one(i)
  ggsave(file.path(OUT, sprintf("pictogram_%s.png", arch$Cluster[i])),
          p, width = 2.0, height = 3.0, dpi = 600, bg = "transparent")
  ggsave(file.path(OUT, sprintf("pictogram_%s.pdf", arch$Cluster[i])),
          p, width = 2.0, height = 3.0, device = cairo_pdf)
}

# Strip of 4 (for inline use)
suppressMessages(library(patchwork))
strip <- draw_one(1) + draw_one(2) + draw_one(3) + draw_one(4) +
  patchwork::plot_layout(nrow = 1)
ggsave(file.path(OUT, "pictograms_strip.png"),
        strip, width = 9, height = 3.2, dpi = 600, bg = "transparent")
ggsave(file.path(OUT, "pictograms_strip.pdf"),
        strip, width = 9, height = 3.2, device = cairo_pdf)

cat("Saved pictograms (transparent backgrounds):\n")
for (f in c(sprintf("pictogram_P%d.png", 1:4), "pictograms_strip.png")) {
  cat("  -", file.path(OUT, f), "\n")
}
