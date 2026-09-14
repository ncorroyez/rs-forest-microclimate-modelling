# Figure — the four coupled processes that set the sub-canopy temperature.
#
# Illustrates Section 1.2.1. Shortwave interception, longwave exchange,
# turbulent transfer and latent heat all act at the height where the foliage
# sits, which is why canopy structure is three-dimensional rather than scalar.
# Two smaller terms close the balance: soil heat storage, and the humidity that
# the same latent flux carries.
#
# This is the physics of buffering, independent of any model. It must NOT drift
# towards Fig_intro_musica, which shows what MuSICA discretises, nor towards the
# MuSICA process inventory reproduced from J. Ogee.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

z  <- seq(0, 1, length.out = 400)
dz <- z[2] - z[1]

GY <- 0.13                      # ground line
TY <- 0.88                      # canopy top
pad <- exp(-((z - 0.62)^2) / (2 * 0.18^2)) * (z > 0.10)
pad <- pad / max(pad)

zc  <- GY + (TY - GY) * z                       # canopy heights on the panel
cum <- rev(cumsum(rev(pad)) * dz)               # leaf area above each height
trans <- exp(-2.2 * cum / max(cum))             # shortwave transmittance
wind  <- exp(-2.6 * cum / max(cum))             # wind speed inside the canopy

# self-checks: the argument the figure makes
stopifnot(all(diff(trans) >= -1e-12))           # transmittance falls downward
stopifnot(trans[1] < 0.25 * trans[length(trans)])
stopifnot(all(diff(wind) >= -1e-12))

XL <- 0.165; WID <- 0.105                       # column footprint

draw_figure <- function() {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2)); on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")

  # --- ground and soil ------------------------------------------------------
  rect(0.035, GY - 0.075, 0.545, GY, col = "#EDE7DE", border = NA)
  segments(0.035, GY, 0.545, GY, col = col_ground, lwd = 1.6)
  text(0.048, GY - 0.040, "soil", cex = 0.50, col = col_ground, adj = c(0, 0.5))

  # --- (1) sunlight: a wedge that narrows as it is intercepted --------------
  xs <- 0.075; w0 <- 0.030
  polygon(c(xs - w0 * trans, rev(xs + w0 * trans)), c(zc, rev(zc)),
          col = "#EBD9A8", border = "#C79A3A", lwd = 0.6)
  polygon(c(xs - w0, xs + w0, xs + w0, xs - w0), c(TY, TY, 0.945, 0.945),
          col = "#EBD9A8", border = "#C79A3A", lwd = 0.6)
  arrows(xs, GY + 0.055, xs, GY + 0.010, length = 0.045, lwd = 1.0,
         col = "#A8801F")
  text(xs, 0.972, "sunlight", cex = 0.50, col = "#A8801F")

  # --- the canopy column ----------------------------------------------------
  polygon(c(rep(XL, length(z)), rev(XL + WID * pad)), c(zc, rev(zc)),
          col = col_leaf, border = col_leaf_d, lwd = 0.7)
  text(XL + WID * 0.45, TY + 0.050, "canopy", cex = 0.53, col = col_leaf_d,
       font = 2)

  # --- (4) latent heat: crown, understory, soil -----------------------------
  for (p4 in list(c(0.196, 0.70), c(0.213, 0.34), c(0.230, GY + 0.015))) {
    arrows(p4[1], p4[2], p4[1], p4[2] + 0.095, length = 0.032, lwd = 0.9,
           col = "#4E7E92", lty = 2)
  }

  # --- (2) longwave: emitted up to the sky and down to the floor ------------
  xw <- 0.300
  arrows(xw, TY - 0.015, xw, 0.945, length = 0.042, lwd = 1.0, col = "#9B8BA8")
  arrows(xw, 0.50, xw, GY + 0.012, length = 0.042, lwd = 1.0, col = "#9B8BA8")
  text(xw + 0.028, 0.972, "longwave", cex = 0.50, col = "#7E6E8C")

  # --- (3) wind profile, decaying into the canopy ---------------------------
  xu <- 0.360
  segments(xu, GY, xu, 0.945, col = col_border, lwd = 0.7, lty = 3)
  lines(xu + 0.135 * wind, zc, col = col_als, lwd = 1.5)
  for (yy in seq(0.19, 0.86, length.out = 6)) {
    i <- which.min(abs(zc - yy))
    arrows(xu, yy, xu + 0.135 * wind[i], yy, length = 0.028, lwd = 0.8,
           col = col_als)
  }
  text(xu + 0.068, 0.972, "wind", cex = 0.50, col = "#5E7B95")

  # --- (5) soil heat storage ------------------------------------------------
  arrows(0.215, GY - 0.012, 0.215, GY - 0.062, length = 0.028, lwd = 0.9,
         col = col_ground, code = 3)

  # --- callouts -------------------------------------------------------------
  LX <- 0.615
  items <- list(
    list(y = 0.930, n = "1", t = "Shortwave interception. Radiation is attenuated with cumulative leaf area, so less energy reaches the floor as the canopy thickens."),
    list(y = 0.745, n = "2", t = "Longwave exchange. Foliage emits downward and screens the understory from a cold sky, so it warms the understory at night."),
    list(y = 0.560, n = "3", t = "Turbulent transfer. Drag decays the wind with depth, insulating the sub-canopy air but also trapping the heat released below it."),
    list(y = 0.375, n = "4", t = "Latent heat. Energy consumed by evaporation, from crown, understory and soil, does not warm the air, and raises its humidity instead."),
    list(y = 0.190, n = "5", t = "Soil heat storage. The ground absorbs by day and releases at night, so soil and air offsets are correlated but not interchangeable.")
  )
  for (it in items) {
    points(LX - 0.022, it$y, pch = 21, bg = "white", col = col_border,
           cex = 1.45, lwd = 0.8)
    text(LX - 0.022, it$y, it$n, cex = 0.50, font = 2, col = col_arrow)
    wrap_text(it$t, width = 52, x = LX, y = it$y + 0.028, cex = 0.50,
              lh = 0.031, col = "black")
  }

  text(0.30, 0.035, paste("Every term acts at the height where the foliage sits,",
       "so a canopy is three-dimensional, not a single number."),
       cex = 0.53, font = 3, col = col_arrow, adj = c(0.5, 0.5))
}

intro_save(draw_figure, "Fig_intro_energy_balance", width = 7.4, height = 4.2)
