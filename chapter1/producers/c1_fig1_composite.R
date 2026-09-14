# ==============================================================================
# Figure 1 (COMPOSITE) — site + typology in one panel set:
#   (a) massif map: 400 cLHS pixels + 53 HOBO loggers, coloured by archetype P1-P4,
#       shape distinguishes cLHS (dot) vs HOBO (diamond); scale bar + north arrow.
#   (b) mean LAD profile per archetype (line + IQR ribbon).
#   (c) VCI per archetype: cLHS boxplots + HOBO points overlaid (marked).
# Sources _article_style.R (PAL_CLUSTER, theme_article, ggsave_article).
#   Rscript c1_fig1_composite.R
# Out: out_files/Chapter1/figures/FigComposite_typology_site.{png,pdf}
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(sf)
  source("R/cluster_relabel.R"); source("scripts/_article_style.R")
})
LV <- c("P1","P2","P3","P4")

# ---- data --------------------------------------------------------------------
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, P := factor(relabel_cluster(Cluster), levels = LV)]
g <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
gc <- data.table(id_plot = g$id_plot, x = g$coord_x_utm31n, y = g$coord_y_utm31n)
hp <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, P = factor(P, levels = LV), VCI)]
hobo <- merge(hp, gc, by = "id_plot")

# ---- (a) map -----------------------------------------------------------------
rng <- function(v) c(min(v), max(v))
xr <- rng(c(samp$x, hobo$x)); yr <- rng(c(samp$y, hobo$y))
sb0 <- xr[1] + 0.04*diff(xr); sby <- yr[1] + 0.04*diff(yr)         # 1 km scale bar
nx  <- xr[2] - 0.04*diff(xr); ny0 <- yr[2] - 0.14*diff(yr)         # north arrow
pA <- ggplot() +
  geom_point(data = samp, aes(x, y, colour = P), shape = 16, size = 0.7, alpha = 0.55) +
  geom_point(data = hobo, aes(x, y, fill = P), shape = 23, size = 2.3, colour = "grey15", stroke = 0.4) +
  annotate("segment", x = sb0, xend = sb0 + 1000, y = sby, yend = sby, linewidth = 0.8) +
  annotate("text", x = sb0 + 500, y = sby, label = "1 km", vjust = -0.6, size = 3) +
  annotate("segment", x = nx, xend = nx, y = ny0, yend = ny0 + 0.10*diff(yr),
           arrow = arrow(length = unit(2.2, "mm")), linewidth = 0.7) +
  annotate("text", x = nx, y = ny0 + 0.115*diff(yr), label = "N", fontface = "bold", size = 3.3) +
  scale_colour_manual(values = PAL_CLUSTER, name = "Archetype") +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  coord_equal(xlim = xr, ylim = yr, expand = FALSE) +
  labs(x = "Easting (m, UTM 31N)", y = "Northing (m)",
       subtitle = "Blois oak forest (47.6° N, 1.3° E). Dots: 400 cLHS pixels. Diamonds: 53 HOBO loggers.") +
  theme_article(11) + theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
                            legend.background = element_rect(fill = "white", colour = "grey80"))

# France locator inset (top-left of the map), if rnaturalearth is available
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  fr <- sf::st_crop(rnaturalearth::ne_countries(country = "France", scale = "medium", returnclass = "sf"),
                    xmin = -5, xmax = 10, ymin = 41, ymax = 52)
  inset <- ggplot(fr) + geom_sf(fill = "grey92", colour = "grey55", linewidth = 0.2) +
    annotate("point", x = 1.33, y = 47.6, colour = "#D7191C", size = 1.7) +
    coord_sf(expand = FALSE) + theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3))
  pA <- pA + patchwork::inset_element(inset, left = 0.0, bottom = 0.66, right = 0.24, top = 1.0, align_to = "panel", ignore_tag = TRUE)
}

# ---- (b) mean LAD profiles ---------------------------------------------------
ladc <- grep("^LAD_Layer_", names(samp), value = TRUE); h <- as.numeric(sub("LAD_Layer_", "", ladc))
L <- melt(samp[, c("P", ladc), with = FALSE], id.vars = "P", variable.name = "lay", value.name = "lad")
L[, z := h[match(lay, ladc)]][is.na(lad), lad := 0]
prof <- L[, .(m = mean(lad), lo = quantile(lad, .25), hi = quantile(lad, .75)), by = .(P, z)]
pB <- ggplot(prof, aes(m, z, colour = P, fill = P)) +
  geom_ribbon(aes(xmin = lo, xmax = hi), alpha = 0.18, colour = NA) +
  geom_path(linewidth = 0.9) +
  scale_colour_manual(values = PAL_CLUSTER, guide = "none") + scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = expression(LAD~(m^2~m^-3)), y = "Height (m)",
       subtitle = "Mean leaf-area-density profile per archetype (ribbon: IQR)") +
  theme_article(11)

# ---- (c) VCI per archetype: cLHS box + HOBO points ---------------------------
pC <- ggplot(samp[is.finite(VCI)], aes(P, VCI)) +
  geom_boxplot(aes(fill = P), width = 0.62, outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(data = hobo[is.finite(VCI)], aes(P, VCI), shape = 23, fill = "white",
              colour = "grey15", size = 1.7, width = 0.12, stroke = 0.35) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = NULL, y = "VCI (vertical complexity)",
       subtitle = "VCI per archetype. Boxes: 400 cLHS pixels. Diamonds: 53 HOBO loggers.") +
  theme_article(11)

# ---- (d) typical ALS point clouds per archetype (rendered externally) --------
suppressPackageStartupMessages({ library(png); library(grid) })
pc_file <- "outputs/figures_chap1_pointclouds/panel_typical.png"
pD <- if (file.exists(pc_file)) {
  patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(pc_file), interpolate = TRUE)) +
    labs(title = "Typical ALS point cloud per archetype (closest to the cluster-mean profile)")
} else patchwork::plot_spacer()

# ---- assemble ----------------------------------------------------------------
fig <- pD / pA / (pB | pC) + plot_layout(heights = c(0.62, 1.15, 1)) + plot_annotation(tag_levels = "a")
ggsave_article("out_files/Chapter1/figures/FigComposite_typology_site", fig, 9.2, 12.4)
cat("DONE -> FigComposite_typology_site.{png,pdf}  (cLHS", nrow(samp), "+ HOBO", nrow(hobo), ")\n")
