# Figure — the two sensors do not read the same part of the canopy column.
#
# Illustrates Sections 1.4 and 1.5.2. Airborne LiDAR recovers the whole column
# because pulses travel to the ground and back, so gap fraction integrates over
# the full profile. An optical satellite sees reflected sunlight, which stops
# carrying information once enough leaf area lies above: below an effective
# optical depth, added foliage barely changes the signal.
#
# Saturation is shown STRUCTURALLY -- the sensed band keeps a comparable
# optical path while the canopy grows, so it covers a smaller and smaller share
# of the column. No saturation curve is plotted, because plotting one would
# invite the figure to be read as a result of Chapter 2.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

z <- seq(0, 1, length.out = 500)
dz <- z[2] - z[1]

# Two canopies of the same shape and height, differing only in leaf quantity.
shape   <- exp(-((z - 0.58)^2) / (2 * 0.20^2)) * (z > 0.08)
shape   <- shape / (sum(shape) * dz)
pad_mod <- shape * 1.0          # moderate canopy
pad_den <- shape * 2.6          # dense canopy, same shape and height

xmax <- max(pad_den) * 1.15

# Effective optical depth: the level below which the optical signal carries
# little further information. Defined here by a fixed cumulative leaf area
# measured downward from the canopy top -- the same optical path in both
# panels, which is precisely why it corresponds to a shallower geometric
# depth, and a smaller share of the column, in the denser canopy.
tau_star <- 0.85

depth_level <- function(pad) {
  cum <- rev(cumsum(rev(pad)) * dz)   # leaf area above each level
  i <- which(cum >= tau_star)
  if (!length(i)) return(min(z))
  z[max(i)]
}

panel <- function(pad, title, note, xoff, w, seen_label) {
  y0 <- 0.300; y1 <- 0.720
  bw <- w * 0.60
  ax <- xoff - 0.014
  px <- function(v) xoff + (v / xmax) * bw
  py <- function(v) y0 + v * (y1 - y0)

  zd <- depth_level(pad)
  keep <- pad > 0.01 * max(pad)
  zk <- z[keep]; pk <- pad[keep]

  # full profile
  polygon(c(px(0), px(pk), px(0)), c(py(min(zk)), py(zk), py(max(zk))),
          col = col_leaf, border = col_leaf_d, lwd = 0.9)

  # the part an optical satellite effectively senses, shaded on top
  up <- keep & z >= zd
  zu <- z[up]; pu <- pad[up]
  polygon(c(px(0), px(pu), px(0)), c(py(min(zu)), py(zu), py(max(zu))),
          col = col_s2, border = col_leaf_d, lwd = 0.9)
  segments(px(0), py(zd), xoff + bw * 1.32, py(zd), lty = 2, lwd = 1,
           col = "#8A6E43")

  segments(ax, py(0), xoff + bw * 1.70, py(0), col = col_ground, lwd = 1.4)
  arrows(ax, py(0), ax, py(1.06), length = 0.038, lwd = 1, col = col_arrow)

  # sunlight in and reflectance out: the return turns back at the depth below
  # which the signal no longer carries information
  xs <- xoff + bw * 1.10
  arrows(xs, py(1.14), xs, py(zd) + 0.006, length = 0.050, lwd = 1.1,
         col = col_s2)
  arrows(xs + bw * 0.16, py(zd) + 0.006, xs + bw * 0.16, py(1.14),
         length = 0.050, lwd = 1.1, col = col_s2)

  # a laser pulse: down to the ground and back, through the whole column
  xl <- xoff + bw * 1.45
  arrows(xl, py(1.14), xl, py(0.015), length = 0.050, lwd = 1.1, col = col_als)
  arrows(xl + bw * 0.16, py(0.015), xl + bw * 0.16, py(1.14), length = 0.050,
         lwd = 1.1, col = col_als)

  text(xoff + bw * 0.80, 0.905, title, cex = 0.64, font = 2)
  text(xoff + bw * 0.80, 0.866, note, cex = 0.535, col = col_arrow)
  text(xoff + bw * 0.80, py(0) - 0.042, seen_label, cex = 0.525, font = 3,
       col = "#8A6E43")
}

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.165, 1))

  w <- 0.34
  panel(pad_mod, "(a) moderate canopy",
        "the sensed band covers most of the profile", 0.075, w,
        "the satellite reads most of the leaf area")
  panel(pad_den, "(b) dense canopy",
        "the same optical path, a much smaller share", 0.560, w,
        "the satellite reads only the upper part")

  text(0.028, 0.50, "height above ground", srt = 90, cex = 0.54,
       col = col_arrow)

  # legend, on one row below the panels so it can never meet a panel title
  ly <- 0.198
  lg <- function(x, fill, lab, dashed = FALSE, line_col = NULL) {
    if (!is.null(line_col)) {
      if (dashed) segments(x, ly, x + 0.022, ly, col = line_col, lwd = 1.2,
                           lty = 2)
      else segments(x, ly, x + 0.022, ly, col = line_col, lwd = 1.6)
    } else {
      rect(x, ly - 0.014, x + 0.022, ly + 0.014, col = fill,
           border = col_leaf_d, lwd = 0.7)
    }
    text(x + 0.030, ly, lab, adj = c(0, 0.5), cex = 0.495)
  }
  lg(0.045, col_s2,   "leaf area the optical signal responds to")
  lg(0.330, col_leaf, "leaf area it does not reach")
  lg(0.530, NULL, "effective optical depth", dashed = TRUE,
     line_col = "#8A6E43")
  lg(0.725, NULL, "LiDAR pulse: ground and back", line_col = col_als)

}

intro_save(draw_figure, "Fig_intro_sensor_column", width = 6.6, height = 3.2)

# ---- what the drawing asserts, checked ------------------------------------
frac <- function(pad) {
  zd <- depth_level(pad)
  sum(pad[z >= zd]) / sum(pad)
}
message(sprintf("share of leaf area above the optical depth: moderate %.0f%%, dense %.0f%%",
                100 * frac(pad_mod), 100 * frac(pad_den)))
stopifnot(frac(pad_mod) > frac(pad_den))
