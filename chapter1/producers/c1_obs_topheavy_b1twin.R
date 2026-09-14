# ==============================================================================
# Observational twin of Figure B1: does OBSERVED buffering respond to how high the
# leaf area sits, over the range real stands actually span?
#
# Committee 2026-08-09 (JO): "une figure séparée mais identique [a B1], avec en x la
# hauteur ponderee du LAI, en y le dTmax observe et en couleur le LAI observe".
#
# B1 sweeps the LAD centroid over the whole 0-1 range in silico at fixed leaf area.
# This figure puts the 53 loggers on the SAME x-axis, so the narrowness of the real
# range is read directly against the swept one. No simulation is involved.
#
# The statistical companion (partial correlation of buffering on top-heaviness given
# LAI and fCover) is scripts/c1_obs_profileshape_test.R; this is the raw picture.
#
# Reads : in_files/lad_z05/Blois_lad_z05_r25.csv                   (LAD profiles)
#         in_files_native20/{lai_z1,fCover,max}_res_10_m.tif       (logger traits)
#         in_files/data_Blois_utm31n.geojson                       (logger locations)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv   (stage A2, READ ONLY)
# Writes: outputs/figures_chap1/fig_obs_topheavy_b1twin.{png,pdf}
#         outputs/figures_chap1/tab_obs_topheavy_b1twin.csv
#   Rscript scripts/c1_obs_topheavy_b1twin.R
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(terra); library(sf)
})
setwd(here::here())
source("scripts/_article_style.R")
OUT <- "outputs/figures_chap1"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- logger traits, native-20 m lineage (same extraction as Fig. 7) ----------
g  <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
r0 <- rast("in_files_native20/lai_z1_res_10_m.tif")
g  <- st_transform(g, crs(r0)); v <- vect(g)
TR <- data.table(id_plot = g$id_plot,
  LAI    = terra::extract(rast("in_files_native20/lai_z1_res_10_m.tif"), v)[, 2],
  fCover = terra::extract(rast("in_files_native20/fCover_res_10_m.tif"), v)[, 2],
  Hmax   = terra::extract(rast("in_files_native20/max_res_10_m.tif"),    v)[, 2])

# ---- leaf-area-weighted height, as a fraction of Hmax (the B1 x-axis) --------
lad <- fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
lyr <- grep("^LAD_Layer_", names(lad), value = TRUE)
hh  <- as.numeric(sub("LAD_Layer_", "", lyr))
Lm  <- as.matrix(lad[, ..lyr]); Lm[!is.finite(Lm)] <- 0
lad[, cen := (Lm %*% hh)[, 1] / pmax(rowSums(Lm), 1e-9)]        # centroid height, m
# Normalize by the LAD file's OWN Hmax, not the native-20 m raster one: the centroid
# and the canopy top must come from the same footprint, and this is the convention
# behind the manuscript's quoted median of 0.48 (Section 4.1, Fig. 7).
lad[, topheavy := cen / pmax(Hmax, 1e-9)]                       # 0 = bottom, 1 = top
TR <- merge(TR, lad[, .(id_plot, cen, topheavy)], by = "id_plot")

obs <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, dTmax_obs)]
D <- merge(TR, obs, by = "id_plot")
D <- D[is.finite(dTmax_obs) & is.finite(topheavy) & is.finite(LAI) & is.finite(fCover)]

qs <- quantile(D$topheavy, c(.25, .5, .75))
cat(sprintf("n = %d | topheaviness median %.2f, IQR %.2f-%.2f, full range %.2f-%.2f\n",
            nrow(D), qs[2], qs[1], qs[3], min(D$topheavy), max(D$topheavy)))
cat(sprintf("occupied fraction of the B1 swept range (0-1): %.0f%% by IQR, %.0f%% by full range\n",
            100 * (qs[3] - qs[1]), 100 * diff(range(D$topheavy))))

# simple and partial association, both reported on the figure
r_s  <- cor(D$dTmax_obs, D$topheavy)
p_s  <- cor.test(D$dTmax_obs, D$topheavy)$p.value
ry   <- residuals(lm(dTmax_obs ~ LAI + fCover, D))
rx   <- residuals(lm(topheavy  ~ LAI + fCover, D))
r_p  <- cor(ry, rx); p_p <- cor.test(ry, rx)$p.value
r_lt <- cor(D$LAI, D$topheavy)
cat(sprintf("simple r(dTmax, topheaviness) = %+.2f (p = %.3f)\n", r_s, p_s))
cat(sprintf("partial r given LAI + fCover  = %+.2f (p = %.3f)\n", r_p, p_p))
cat(sprintf("collinearity r(LAI, topheaviness) = %+.2f\n", r_lt))

# ---- figure ------------------------------------------------------------------
# x spans the FULL B1 sweep (0-1) on purpose: the empty margins are the message.
#' Format a p value in house style
#' @param p numeric p value
#' @return character, "< 0.001" below that threshold and two decimals above
fp  <- function(p) if (p < 0.001) "*p* < 0.001" else sprintf("*p* = %.2f", p)
ann <- sprintf("simple *r* = %+.2f (%s); partial *r* = %+.2f given LAI and fCover (%s)",
               r_s, fp(p_s), r_p, fp(p_p))
p <- ggplot(D, aes(topheavy, dTmax_obs)) +
  annotate("rect", xmin = qs[1], xmax = qs[3], ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.7) +
  geom_point(aes(fill = LAI), shape = 21, size = 3.1, stroke = 0.3, colour = "grey20") +
  scale_fill_viridis_c(name = "observed LAI\n(one-sided)", option = "viridis", direction = -1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "leaf-area-weighted height  (LAD centroid / *H*<sub>max</sub>: 0 = bottom-heavy, 1 = top-heavy)",
       y = expression(observed~Delta*T[max]~"(°C)"),
       subtitle = ann) +
  theme_article(12) +
  theme(axis.title.x = ggtext::element_markdown(),
        plot.subtitle = ggtext::element_markdown(),
        legend.position = c(0.90, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))
ggsave_article(file.path(OUT, "fig_obs_topheavy_b1twin"), p, 7.6, 5.4)
fwrite(D[, .(id_plot, dTmax_obs, LAI, fCover, Hmax, cen, topheavy)],
       file.path(OUT, "tab_obs_topheavy_b1twin.csv"))
cat("DONE\n")
