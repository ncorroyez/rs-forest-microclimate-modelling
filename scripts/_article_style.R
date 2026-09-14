# ==============================================================================
# Style partagé des figures annexes (Chap 1) — sourcé par les scripts make_*.R
# pour garantir une présentation homogène (thème, palettes sémantiques, tailles).
# ==============================================================================
suppressMessages(library(ggplot2))

# --- palettes sémantiques (mêmes couleurs pour le même concept partout) -------
PAL_GRP     <- c(buffering = "#1A9850", amplifying = "#D7191C")          # buff = vert, amp = rouge
PAL_CLUSTER <- c(P1 = "#D7191C", P2 = "#FDAE61", P3 = "#74C476", P4 = "#1A9850")  # clair→dense
PAL_SEQ2    <- c("#1b9e77", "#d95f02")                                   # 2 niveaux ordonnés (ex. LAI)

# --- thème maison -------------------------------------------------------------
theme_article <- function(base = 12) {
  theme_bw(base_size = base) +
    theme(panel.grid       = element_blank(),
          strip.background = element_blank(),   # house style: no grey strip fill (clean facets)
          strip.text       = element_text(face = "bold", margin = margin(b = 4)),
          plot.subtitle    = element_blank(),   # house style: no in-figure subtitles (captions live in the manuscript)
          plot.caption     = element_text(size = 9,  colour = "grey45"),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key       = element_blank())
}

# --- légende intégrée dans un coin (placement homogène) -----------------------
legend_corner <- function(x = 0.99, y = 0.99) {
  theme(legend.position = c(x, y),
        legend.justification = c(if (x > 0.5) 1 else 0, if (y > 0.5) 1 else 0))
}

# --- sauvegarde homogène (png + pdf, 300 dpi, fond blanc) ---------------------
ggsave_article <- function(path_noext, plot, width, height) {
  ggsave(paste0(path_noext, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
}
