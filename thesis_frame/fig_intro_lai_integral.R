# Figure — LAI is an integral, and integrals discard information.
#
# Illustrates Section 1.2.2: three canopies that share the same leaf area index
# but differ in the vertical arrangement of that leaf area and in the fraction
# of ground it covers. In (a) and (b) the drawn areas are equal by construction
# (checked at the end of this script). In (c) the drawn area is LARGER, because
# it is the density where foliage is present: the leaf area per unit GROUND is
# what is held equal. The panel and the caption both say so.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

# ---- three leaf area density profiles, normalised to the same integral -----
# z runs from 0 (ground) to 1 (top of the tallest canopy).
z <- seq(0, 1, length.out = 400)

# (a) shallow: foliage concentrated in a narrow upper crown layer
pad_a <- exp(-((z - 0.80)^2) / (2 * 0.055^2))
# (b) deep: foliage spread through a multi-layered column
pad_b <- exp(-((z - 0.55)^2) / (2 * 0.230^2)) * (z > 0.06)
# (c) same shape as (b) but present over only part of the ground
pad_c <- pad_b

LAI <- 1                                  # arbitrary common integral
area <- function(p) sum(p) * (z[2] - z[1])
pad_a <- pad_a / area(pad_a) * LAI
pad_b <- pad_b / area(pad_b) * LAI
cover_c <- 0.55                           # fraction of ground covered
pad_c <- pad_c / area(pad_c) * LAI / cover_c   # denser where present

xmax <- max(pad_a, pad_b, pad_c) * 1.12

panel <- function(pad, cover, title, subtitle, xoff, w) {
  # profile drawn in [xoff, xoff + bw] x [y0, y1]; the axis arrow sits to the
  # left of xoff so it can never run into the panel text.
  y0 <- 0.315; y1 <- 0.815
  bw <- w * 0.62
  ax <- xoff - 0.016
  px <- function(v) xoff + (v / xmax) * bw
  py <- function(v) y0 + v * (y1 - y0)

  # ground line
  segments(ax, py(0), xoff + bw * 1.06, py(0), col = col_ground, lwd = 1.4)

  # the profile, restricted to the height range where foliage is present, so
  # that the crown base is visible rather than a hairline down to the ground
  keep <- pad > 0.01 * max(pad)
  zk <- z[keep]; pk <- pad[keep]
  polygon(c(px(0), px(pk), px(0)),
          c(py(min(zk)), py(zk), py(max(zk))),
          col = col_leaf, border = col_leaf_d, lwd = 0.9)

  # horizontal extent of the canopy over the ground, drawn as a bar under the
  # ground line: this is the dimension LAI also integrates away
  rect(xoff, py(0) - 0.058, xoff + bw, py(0) - 0.028,
       col = "#EDEAE4", border = col_border, lwd = 0.7)
  rect(xoff, py(0) - 0.058, xoff + bw * cover, py(0) - 0.028,
       col = col_leaf, border = col_leaf_d, lwd = 0.7)
  text(xoff + bw / 2, py(0) - 0.082,
       if (cover < 1) "part of the ground only" else "the whole ground",
       cex = 0.50, col = col_arrow)
  # (c) is denser where foliage is present: its drawn area is therefore larger
  # than in (a) and (b). Say so, rather than let the eye read more leaf area.
  if (cover < 1)
    text(xoff + bw / 2, py(0) - 0.106, "(so locally denser)", cex = 0.50,
         font = 3, col = col_arrow)

  # height axis, no tick values
  arrows(ax, py(0), ax, py(1.06), length = 0.042, lwd = 1, col = col_arrow)

  # title and subtitle, centred on the profile, in a band reserved above it
  text(xoff + bw / 2, 0.960, title, cex = 0.66, font = 2)
  ln <- strwrap(subtitle, width = 27)
  for (i in seq_along(ln))
    text(xoff + bw / 2, 0.912 - (i - 1) * 0.030, ln[i], adj = c(0.5, 0.5),
         cex = 0.545)
}

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.16, 1))

  w <- 0.285
  panel(pad_a, 1.00, "(a) shallow",
        "leaf area in a narrow upper crown layer", 0.085, w)
  panel(pad_b, 1.00, "(b) deep",
        "the same leaf area through a deep column", 0.405, w)
  panel(pad_c, cover_c, "(c) clumped",
        "the same leaf area over part of the ground", 0.725, w)

  # axis labels
  text(0.034, 0.565, "height above ground", srt = 90, cex = 0.55,
       col = col_arrow)
  text(0.50, 0.185, "leaf area density  →   (shared axis)", cex = 0.55,
       col = col_arrow)

  # the statement the figure exists to make
}

intro_save(draw_figure, "Fig_intro_lai_integral", width = 6.6, height = 3.3)

# ---- sanity check: the three integrals really are equal --------------------
stopifnot(abs(area(pad_a) - LAI) < 1e-9,
          abs(area(pad_b) - LAI) < 1e-9,
          abs(area(pad_c) * cover_c - LAI) < 1e-9)
message("check: the three leaf area integrals are equal")
