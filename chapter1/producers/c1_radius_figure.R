# ==============================================================================
# Combined circle-vs-square radius-sensitivity figure (annex).
# Overlays per-radius R2 / RMSE / MAE (MuSICA Tair @1m vs HOBO, hourly) for both
# footprint geometries; highlights the chosen 12.5 m radius (25 m extent).
#   Rscript c1_radius_figure.R
# Out: out_files/Chapter1/figures/FigAnnex_radius_sensitivity.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
tdir <- "out_files/Chapter1/tables"; fdir <- "out_files/Chapter1/figures"
dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

cir <- fread(file.path(tdir, "radius_metrics_circle_summary.csv"))[, geom := "Circle (radius r)"]
sqr <- fread(file.path(tdir, "radius_metrics_square_summary.csv"))[, geom := "Square (side 2r)"]
d <- rbind(cir, sqr)
pal <- c("Circle (radius r)" = "#2C7BB6", "Square (side 2r)" = "#D7191C")
CHOSEN <- 12.5

mk <- function(yvar, ylab) {
  ggplot(d, aes(radius, get(yvar), colour = geom)) +
    geom_vline(xintercept = CHOSEN, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = c(5,10,12.5,15,20,25,30,40,50)) +
    labs(x = "Footprint radius r (m)", y = ylab) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

g <- mk("R2", expression(R^2)) + mk("RMSE", "RMSE (°C)") + mk("MAE", "MAE (°C)") +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "top")
ggsave(file.path(fdir, "FigAnnex_radius_sensitivity.png"), g, width = 11, height = 4, dpi = 200, bg = "white")
cat("DONE ->", file.path(fdir, "FigAnnex_radius_sensitivity.png"), "\n")
