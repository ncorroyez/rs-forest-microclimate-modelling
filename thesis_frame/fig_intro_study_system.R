# Figure — the study system.
#
# Opens Section 1.7. The only figure of the introduction that carries measured
# values, and it is declared as such: it locates the three forests the thesis
# works in, and the sensor network one of them carries.
#
# (a) The three national state forests, on a metropolitan France outline.
#     Chapter 2 uses all three; Chapters 1 and 3 use Blois alone.
# (b) Blois, with the 53 temperature loggers whose records every simulation in
#     Chapters 1 and 3 is answerable to. A single colour: the structural
#     typology of those plots is a result of Chapter 1, not of this figure.
#
# Sources. Site coordinates: Chapter 2, Section 2.1. Logger positions:
# in_files/data_Blois_utm31n.geojson (60 plots), restricted to the 53 retained
# in the analysis by outputs/figures_pipeline/data/clusters.rds. That 53-plot
# membership is identical across the three pipeline lineages present in the
# repository, so the filter does not commit the figure to any of them.
# France outline: ~/Documents/Presentations/france.geojson.

suppressPackageStartupMessages({library(sf)})
source(here::here("thesis_frame", "_intro_fig_style.R"))

RM   <- here::here()
PRES <- path.expand("~/Documents/Presentations")

# ---- (a) metropolitan France ------------------------------------------------
fr <- st_read(file.path(PRES, "france.geojson"), quiet = TRUE)
fr <- st_make_valid(fr)
metro <- suppressWarnings(
  st_crop(fr, st_bbox(c(xmin = -5.5, ymin = 41.0, xmax = 10.0, ymax = 51.5),
                      crs = st_crs(fr))))
metro <- st_transform(st_union(metro), 2154)          # Lambert-93

sites <- data.frame(
  name = c("Mormal", "Blois", "Aigoual"),
  desc = c("lowland oak and beech", "lowland sessile oak", "montane beech"),
  lon  = c(3.74, 1.29, 3.52),
  lat  = c(50.20, 47.57, 44.12),
  # label placement: side of the point the text hangs on
  side = c("right", "left", "right"), stringsAsFactors = FALSE)
sp <- st_transform(st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), 2154)
sxy <- st_coordinates(sp)

# ---- (b) Blois and its logger network ---------------------------------------
plots <- st_read(file.path(RM, "in_files", "data_Blois_utm31n.geojson"),
                 quiet = TRUE)
keep  <- readRDS(file.path(RM, "outputs", "figures_pipeline", "data",
                           "clusters.rds"))$id_plot
plots <- plots[plots$id_plot %in% keep, ]
pxy   <- cbind(plots$coord_x_utm31n, plots$coord_y_utm31n)

draw_figure <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1.12))

  # ---- panel (a) ------------------------------------------------------------
  par(mar = c(0.4, 0.4, 2.2, 0.4), xaxs = "i", yaxs = "i")
  bb <- st_bbox(metro)
  padx <- 0.24 * (bb["xmax"] - bb["xmin"])
  plot.new()
  plot.window(xlim = c(bb["xmin"] - 0.02 * padx, bb["xmax"] + padx),
              ylim = c(bb["ymin"], bb["ymax"]), asp = 1)
  plot(metro, add = TRUE, col = "#F2F0EC", border = col_border, lwd = 0.7)
  points(sxy[, 1], sxy[, 2], pch = 21, bg = col_leaf, col = col_leaf_d,
         cex = 1.25, lwd = 0.9)
  dxl <- 0.030 * (bb["xmax"] - bb["xmin"])
  for (i in seq_len(nrow(sites))) {
    at <- if (sites$side[i] == "right") 4 else 2
    off <- if (at == 4) dxl else -dxl
    text(sxy[i, 1] + off, sxy[i, 2] + 0.028 * (bb["ymax"] - bb["ymin"]),
         sites$name[i], adj = c(if (at == 4) 0 else 1, 0.5), cex = 0.62,
         font = 2)
    text(sxy[i, 1] + off, sxy[i, 2] - 0.028 * (bb["ymax"] - bb["ymin"]),
         sites$desc[i], adj = c(if (at == 4) 0 else 1, 0.5), cex = 0.545,
         col = col_arrow)
  }
  mtext("(a)  three national state forests", side = 3, line = 0.7, adj = 0,
        cex = 0.66, font = 2)
  mtext("Chapter 2 uses all three; Chapters 1 and 3 use Blois", side = 3,
        line = -0.1, adj = 0, cex = 0.545, col = col_arrow)
  # a 100 km bar, in Lambert-93 metres
  x0 <- bb["xmin"] + 0.06 * (bb["xmax"] - bb["xmin"])
  y0 <- bb["ymin"] + 0.06 * (bb["ymax"] - bb["ymin"])
  segments(x0, y0, x0 + 1e5, y0, lwd = 1.4, col = col_txt <- "#3A3A3A")
  segments(c(x0, x0 + 1e5), y0 - 8e3, c(x0, x0 + 1e5), y0 + 8e3, lwd = 1.4)
  text(x0 + 5e4, y0 + 3.4e4, "100 km", cex = 0.545)

  # ---- panel (b) ------------------------------------------------------------
  par(mar = c(0.4, 0.4, 2.2, 0.4))
  rx <- range(pxy[, 1]); ry <- range(pxy[, 2])
  mx <- 0.10 * diff(rx); my <- 0.10 * diff(ry)
  plot.new()
  plot.window(xlim = c(rx[1] - mx, rx[2] + mx), ylim = c(ry[1] - my, ry[2] + my),
              asp = 1)
  points(pxy[, 1], pxy[, 2], pch = 21, bg = col_als, col = "#5D778E",
         cex = 0.85, lwd = 0.7)
  mtext("(b)  Blois: the 53 temperature loggers", side = 3, line = 0.7,
        adj = 0, cex = 0.66, font = 2)
  mtext(sprintf("a network spanning about %.0f by %.0f km", diff(range(pxy[, 1])) / 1000,
                diff(range(pxy[, 2])) / 1000), side = 3,
        line = -0.1, adj = 0, cex = 0.545, col = col_arrow)
  bx <- rx[1] - 0.4 * mx; by <- ry[1] - 0.4 * my
  segments(bx, by, bx + 2000, by, lwd = 1.4)
  segments(c(bx, bx + 2000), by - 130, c(bx, bx + 2000), by + 130, lwd = 1.4)
  text(bx + 1000, by + 480, "2 km", cex = 0.545)
  # north arrow
  nx <- rx[2] + 0.5 * mx; ny <- ry[2] - 0.10 * diff(ry)
  arrows(nx, ny, nx, ny + 0.09 * diff(ry), length = 0.05, lwd = 1.2)
  text(nx, ny + 0.115 * diff(ry), "N", cex = 0.60, font = 2)
}

intro_save(draw_figure, "Fig_intro_study_system", width = 6.4, height = 3.5)

# ---- sanity checks ----------------------------------------------------------
stopifnot(nrow(plots) == 53)                       # the analysed network, not 60
stopifnot(all(st_within(sp, st_union(metro), sparse = FALSE)[, 1]))
message(sprintf(
  "check: %d loggers plotted; all %d sites fall inside the France outline; Blois extent %.1f x %.1f km",
  nrow(plots), nrow(sites), diff(range(pxy[, 1])) / 1000,
  diff(range(pxy[, 2])) / 1000))
