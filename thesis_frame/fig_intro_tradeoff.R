# Figure — the spatio-temporal trade-off between structural sensors.
#
# Illustrates Section 1.5.3. Positions each sensor on two axes: how much of the
# three-dimensional canopy it recovers, and how often it revisits. The corner
# a mechanistic microclimate model would want -- full structure, repeated -- is
# empty, and that emptiness is the third gap. Coverage is a third property, and
# it is stated in words under each sensor rather than encoded in a symbol size
# that would invite quantitative reading.
#
# Conceptual figure. Positions are qualitative; no axis carries tick values.

source(here::here("thesis_frame", "_intro_fig_style.R"))

# plot region
X0 <- 0.135; X1 <- 0.935
Y0 <- 0.300; Y1 <- 0.855

# GEDI is shown hollow and greyed. It has to appear, or the empty corner is open
# to the objection that spaceborne LiDAR already fills it; but it is a General
# Discussion perspective and is used in none of the three chapters, so it must
# not carry the same visual weight as the two sensors the thesis actually works
# with.
sensors <- list(
  list(x = 0.16, y = 0.86, col = col_als,  muted = FALSE, lab = "Airborne LiDAR",
       sub = c("full vertical profile", "one campaign, limited area")),
  list(x = 0.47, y = 0.60, col = col_gedi, muted = TRUE,  lab = "GEDI",
       sub = c("profile along sparse footprints, not wall-to-wall",
               "perspective only, used in no chapter")),
  list(x = 0.88, y = 0.17, col = col_s2,   muted = FALSE, lab = "Sentinel-2",
       sub = c("upper canopy only, saturating", "wall-to-wall, every few days"))
)

ux <- function(u) X0 + u * (X1 - X0)
uy <- function(u) Y0 + u * (Y1 - Y0)

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.18, 1))

  # the empty corner, drawn first so markers sit above it
  rect(ux(0.62), uy(0.62), ux(1.04), uy(1.04), col = "#F7F1F1",
       border = "#C08A8A", lwd = 0.9, lty = 2)
  text(ux(0.83), uy(0.93), "what a mechanistic model", cex = 0.53,
       col = "#9C5F5F", font = 2)
  text(ux(0.83), uy(0.855), "over a landscape would need", cex = 0.53,
       col = "#9C5F5F", font = 2)
  text(ux(0.83), uy(0.75), "no current sensor is here", cex = 0.53,
       col = "#9C5F5F", font = 3)

  # axes
  arrows(X0, Y0, X1 + 0.012, Y0, length = 0.05, lwd = 1.1, col = col_arrow)
  arrows(X0, Y0, X0, Y1 + 0.020, length = 0.05, lwd = 1.1, col = col_arrow)

  text((X0 + X1) / 2, Y0 - 0.088, "revisit frequency", cex = 0.60,
       font = 2, col = col_arrow)
  text(X0 + 0.030, Y0 - 0.048, "acquired once", cex = 0.52, adj = c(0, 0.5),
       col = col_arrow)
  text(X1, Y0 - 0.048, "every few days", cex = 0.52, adj = c(1, 0.5),
       col = col_arrow)

  text(X0 - 0.062, (Y0 + Y1) / 2, "structural information", srt = 90,
       cex = 0.60, font = 2, col = col_arrow)
  text(X0 - 0.026, Y0 + 0.030, "two-dimensional", srt = 90, cex = 0.52,
       adj = c(0, 0.5), col = col_arrow)
  text(X0 - 0.026, Y1 - 0.030, "full 3D profile", srt = 90, cex = 0.52,
       adj = c(1, 0.5), col = col_arrow)

  # sensors
  for (s in sensors) {
    x <- ux(s$x); y <- uy(s$y)
    points(x, y, pch = 21, bg = if (s$muted) "#FFFFFF" else s$col,
           col = if (s$muted) "#A9A9A9" else col_border,
           cex = if (s$muted) 1.9 else 2.3, lwd = if (s$muted) 1.2 else 0.9)
    if (s$muted)
      points(x, y, pch = 21, bg = s$col, col = NA, cex = 0.9)
    # label placed away from the axes so it never leaves the plot region
    adj <- if (s$x > 0.7) c(1, 0.5) else c(0, 0.5)
    dx  <- if (s$x > 0.7) -0.022 else 0.022
    lab_col <- if (s$muted) "#8A8A8A" else "black"
    text(x + dx, y + 0.030, s$lab, adj = adj, cex = if (s$muted) 0.545 else 0.60,
         font = 2, col = lab_col)
    for (i in seq_along(s$sub))
      text(x + dx, y - 0.004 - (i - 1) * 0.030, s$sub[i], adj = adj,
           cex = 0.485, font = if (s$muted && i == 2) 3 else 1,
           col = if (s$muted) "#8A8A8A" else col_arrow)
  }

}

intro_save(draw_figure, "Fig_intro_tradeoff", width = 6.6, height = 3.55)
