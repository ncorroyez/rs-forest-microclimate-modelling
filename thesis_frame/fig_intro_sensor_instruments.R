# Figure — the two sensors as instruments.
#
# Opens Section 1.4. The rest of the introduction compares what the two sensors
# see (Fig. on the canopy column) and what they cost in space and time (Fig. on
# the trade-off). This one compares what they physically do: one carries its own
# light and times its return, the other collects sunlight the canopy sends back.
# Everything downstream — the vertical profile on one side, the saturating
# two-dimensional signal on the other — follows from that single difference.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

rows <- data.frame(
  lab = c("energy source",
          "what is measured",
          "one acquisition returns",
          "what limits it"),
  als = c("its own laser pulse",
          "the time each return takes to come back",
          "many ranged returns per pulse, in three dimensions",
          "the flight plan, the pulse density, the cost"),
  s2  = c("sunlight the canopy reflects",
          "radiance in each spectral band",
          "one value per band over a 10 or 20 m pixel",
          "cloud cover, and the sun and view angles"),
  stringsAsFactors = FALSE)

col_als_d <- "#5D778E"
col_s2_d  <- "#8A6E43"

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))

  xl <- 0.030                       # label column, left edge
  xa <- 0.255; xb <- 0.615          # the two sensor columns
  cw <- 0.355                       # column width
  ytop <- 0.950; ysk <- 0.645       # sketch band: ysk .. ytop

  # ---- column headers and sketch panels -------------------------------------
  sketch <- function(x0, col, dark, title) {
    rect(x0, ysk, x0 + cw, ytop, col = col_panel, border = col_border,
         lwd = 0.7)
    text(x0 + cw / 2, ytop + 0.030, title, cex = 0.68, font = 2, col = dark)
  }
  sketch(xa, col_als, col_als_d, "airborne LiDAR")
  sketch(xb, col_s2,  col_s2_d,  "Sentinel-2")
  wrap_text("how one acquisition works", 16, xl, (ysk + ytop) / 2 + 0.028,
            cex = 0.545, lh = 0.032, font = 2, col = col_arrow)

  # a canopy silhouette, reused in both sketches
  crown <- function(x0, cx, base, top, rx, gy) {
    segments(x0 + cx, gy, x0 + cx, base + 0.020, col = "#7A6A50", lwd = 1.4)
    zz <- seq(base, top, length.out = 60)
    hw <- rx * sqrt(pmax(0, 1 - ((zz - (base + top) / 2) / ((top - base) / 2))^2))
    polygon(x0 + c(cx - hw, rev(cx + hw)), c(zz, rev(zz)),
            col = col_leaf, border = col_leaf_d, lwd = 0.6)
  }
  gy <- ysk + 0.075                 # ground line inside a sketch
  for (x0 in c(xa, xb)) {
    segments(x0 + 0.020, gy, x0 + cw - 0.020, gy, col = col_ground, lwd = 1.2)
    crown(x0, 0.095, gy + 0.026, gy + 0.092, 0.046, gy)
    crown(x0, 0.190, gy + 0.030, gy + 0.112, 0.053, gy)
    crown(x0, 0.283, gy + 0.024, gy + 0.086, 0.044, gy)
  }

  # (a) the aircraft carries its own light, and the returns come back ranged
  ax <- xa + 0.190; ay <- ytop - 0.040
  points(ax, ay, pch = 17, cex = 1.0, col = col_als_d)
  text(ax, ay + 0.024, "aircraft", cex = 0.455, col = col_als_d)
  gx <- xa + 0.190                        # the pulse goes straight down
  arrows(gx - 0.020, ay - 0.016, gx - 0.020, gy + 0.004, length = 0.028,
         lwd = 1, col = col_als_d)
  arrows(gx + 0.020, gy + 0.004, gx + 0.020, ay - 0.016, length = 0.028,
         lwd = 1, col = col_als_d)
  # each interception on the way down is a return with its own height
  for (rz in c(0.104, 0.078, 0.048, 0.004))
    points(gx - 0.020, gy + rz, pch = 16, cex = 0.5, col = col_als_d)
  text(xa + cw / 2, ysk + 0.014, "it ranges: every return has a height",
       cex = 0.455, font = 3, col = col_als_d)

  # (b) the sun supplies the light, the satellite collects what comes back up
  sx <- xb + 0.070; sy <- ytop - 0.042
  points(sx, sy, pch = 8, cex = 1.0, col = "#C9A200")
  text(sx, sy + 0.024, "sun", cex = 0.455, col = "#9A7E00")
  qx <- xb + 0.285
  points(qx, sy, pch = 15, cex = 0.8, col = col_s2_d)
  text(qx, sy + 0.024, "satellite", cex = 0.455, col = col_s2_d)
  arrows(sx + 0.016, sy - 0.018, xb + 0.166, gy + 0.126, length = 0.028,
         lwd = 1, col = "#C9A200")
  arrows(xb + 0.214, gy + 0.126, qx - 0.014, sy - 0.018, length = 0.028,
         lwd = 1, col = col_s2_d)
  # the footprint over which everything above it is averaged into one number
  fx0 <- xb + 0.075; fx1 <- xb + 0.305; fy <- gy - 0.010
  segments(fx0, fy, fx1, fy, col = col_s2_d, lwd = 1.2)
  segments(c(fx0, fx1), fy - 0.008, c(fx0, fx1), fy + 0.008, col = col_s2_d,
           lwd = 1.2)
  text((fx0 + fx1) / 2, fy - 0.024, "one pixel", cex = 0.455, col = col_s2_d)
  text(xb + cw / 2, ysk + 0.014, "it integrates: one value over the pixel",
       cex = 0.455, font = 3, col = col_s2_d)

  # ---- the attribute rows ---------------------------------------------------
  y <- 0.550; dy <- 0.116
  for (i in seq_len(nrow(rows))) {
    yy <- y - (i - 1) * dy
    if (i %% 2 == 1)
      rect(xl - 0.006, yy - dy * 0.46, xb + cw, yy + dy * 0.46,
           col = "#F4F2EE", border = NA)
    text(xl, yy, rows$lab[i], adj = c(0, 0.5), cex = 0.545, font = 2,
         col = col_arrow)
    wrap_text(rows$als[i], 30, xa + 0.012, yy + 0.024, cex = 0.545, lh = 0.030)
    wrap_text(rows$s2[i],  30, xb + 0.012, yy + 0.024, cex = 0.545, lh = 0.030)
  }
  segments(xa - 0.008, y + dy * 0.46, xa - 0.008, y - 3.54 * dy,
           col = "#DAD6CE", lwd = 0.8)
  segments(xb - 0.008, y + dy * 0.46, xb - 0.008, y - 3.54 * dy,
           col = "#DAD6CE", lwd = 0.8)

  # ---- the statement the figure exists to make ------------------------------
  text(0.5, 0.038,
       "one instrument is active and resolves height; the other is passive and resolves colour",
       cex = 0.575, font = 3, col = "#5A5A5A")
}

intro_save(draw_figure, "Fig_intro_sensor_instruments", width = 6.2, height = 4.0)

# ---- sanity check -----------------------------------------------------------
# The figure is a paired comparison: every attribute must be stated for both
# sensors, and no cell may be blank.
stopifnot(nrow(rows) == 4,
          all(nzchar(rows$lab)), all(nzchar(rows$als)), all(nzchar(rows$s2)))
message(sprintf("check: %d attributes, both sensors, no empty cell", nrow(rows)))
