# Figure — how a LiDAR leaf area index is made, and where its assumptions enter.
#
# Illustrates Section 1.4.1. Four steps, left to right: (a) the normalised point
# cloud; (b) the fraction of pulses still unintercepted at each height, read off
# the cloud layer by layer; (c) the Beer-Lambert inversion of that gap fraction
# into a plant area density, which is where the extinction coefficient k and the
# assumption of randomly distributed foliage enter; (d) the resulting profile,
# whose integral is the leaf area index.
#
# The point of the figure is that the number at the end is a model-based
# estimate resting on stated assumptions, not a direct measurement. The two
# assumptions are therefore annotated on the step that introduces them, not in
# the caption.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

# ---- a canopy, described once and reused by every panel ---------------------
z   <- seq(0, 1, length.out = 400)
pad <- exp(-((z - 0.62)^2) / (2 * 0.17^2)) * (z > 0.10)
pad <- pad / (sum(pad) * (z[2] - z[1]))          # unit integral

# cumulative interception from the top down, then the gap fraction that follows
k_ext <- 0.5
cum   <- rev(cumsum(rev(pad))) * (z[2] - z[1])   # PAI above height z
gap   <- exp(-k_ext * cum)                       # Beer-Lambert, top-down

# A point cloud drawn as four crowns over their stems, rather than as noise in a
# box: the vertical density still follows the profile above, but the eye has to
# recognise a canopy for panel (a) to mean anything.
set.seed(4)
crown <- data.frame(cx = c(0.19, 0.42, 0.63, 0.85),
                    zc = c(0.66, 0.60, 0.70, 0.57),   # centre of the crown
                    rz = c(0.30, 0.34, 0.27, 0.31),   # vertical half-extent
                    rx = c(0.11, 0.13, 0.10, 0.12))   # horizontal half-extent
np  <- 900
zp  <- sample(z, np, replace = TRUE, prob = pad)
ci  <- sample(nrow(crown), np, replace = TRUE)
# horizontal half-width of a crown at height z, zero outside its vertical extent
hw  <- pmax(0, 1 - ((zp - crown$zc[ci]) / crown$rz[ci])^2)
keep_p <- hw > 0
zp  <- zp[keep_p]; ci <- ci[keep_p]; hw <- hw[keep_p]
xp  <- crown$cx[ci] + runif(length(zp), -1, 1) * crown$rx[ci] * sqrt(hw)

# stems, and the ground returns that make the normalisation possible
ns  <- 26
zs  <- rep(seq(0.02, 1, length.out = ns), each = 4) *
       rep(crown$zc - 0.5 * crown$rz, times = ns)
xs_ <- rep(crown$cx, times = ns) + rnorm(4 * ns, 0, 0.004)
zg  <- runif(90, 0, 0.012); xg <- runif(90, 0.06, 0.94)

col_pulse <- "#8FA0AE"

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))

  y0 <- 0.235; h <- 0.480                        # common vertical band
  w  <- 0.196
  xs <- c(0.045, 0.283, 0.521, 0.759)
  PY <- function(v) y0 + v * h

  head <- function(x0, tag, title) {
    text(x0 + w / 2, y0 + h + 0.104, tag, cex = 0.62, font = 2)
    ln <- strwrap(title, width = 26)
    for (i in seq_along(ln))
      text(x0 + w / 2, y0 + h + 0.062 - (i - 1) * 0.030, ln[i], cex = 0.545,
           col = col_arrow)
  }
  ground <- function(x0) segments(x0, PY(0), x0 + w, PY(0), col = col_ground,
                                  lwd = 1.3)

  # ---- (a) the normalised point cloud ---------------------------------------
  head(xs[1], "(a)", "normalised point cloud")
  points(xs[1] + xs_ * w, PY(zs), pch = 16, cex = 0.14, col = "#8A7A5E")
  points(xs[1] + xp * w, PY(zp), pch = 16, cex = 0.16, col = "#6E8258")
  points(xs[1] + xg * w, PY(zg), pch = 16, cex = 0.16, col = col_ground)
  ground(xs[1])
  arrows(xs[1] - 0.016, PY(0), xs[1] - 0.016, PY(1.04), length = 0.034,
         lwd = 1, col = col_arrow)
  text(xs[1] - 0.036, PY(0.5), "height above ground", srt = 90, cex = 0.50,
       col = col_arrow)

  # ---- (b) gap fraction, layer by layer -------------------------------------
  head(xs[2], "(b)", "gap fraction at each height")
  # horizontal bars: the share of pulses still unintercepted at that height
  zb <- seq(0.03, 0.98, by = 0.052)
  gb <- approx(z, gap, zb)$y
  for (i in seq_along(zb))
    rect(xs[2], PY(zb[i]) - 0.010, xs[2] + gb[i] * w, PY(zb[i]) + 0.010,
         col = col_als, border = "#5D778E", lwd = 0.4)
  ground(xs[2])
  text(xs[2] + w / 2, PY(0) - 0.036, "fraction of pulses that get through",
       cex = 0.48, col = col_arrow)

  # ---- (c) the inversion ----------------------------------------------------
  head(xs[3], "(c)", "Beer-Lambert inversion")
  bx <- xs[3] + w / 2
  text(bx, PY(0.84), expression(P[gap](z) == e^{-k * PAI(z)}), cex = 0.72)
  text(bx, PY(0.71), expression(PAI(z) == -ln(P[gap](z)) / k), cex = 0.72)
  text(bx, PY(0.62), "then differentiate with height", cex = 0.475,
       col = col_arrow)
  # the two assumptions, on the step that introduces them
  rect(xs[3] + 0.004, PY(0.30), xs[3] + w - 0.004, PY(0.56),
       col = "#F6EFE4", border = "#C2A97F", lwd = 0.7)
  wrap_text("k fixes the leaf angle distribution at a nominal value",
            26, xs[3] + 0.014, PY(0.52), cex = 0.475, lh = 0.028,
            col = "#7A5F30")
  rect(xs[3] + 0.004, PY(0.02), xs[3] + w - 0.004, PY(0.26),
       col = "#F6EFE4", border = "#C2A97F", lwd = 0.7)
  wrap_text("foliage is assumed randomly placed within each layer",
            26, xs[3] + 0.014, PY(0.22), cex = 0.475, lh = 0.028,
            col = "#7A5F30")

  # ---- (d) the profile and its integral -------------------------------------
  head(xs[4], "(d)", "plant area density profile")
  pmax_ <- max(pad) * 1.10
  keep  <- pad > 0.005 * max(pad)
  zk <- z[keep]; pk <- pad[keep]
  polygon(xs[4] + c(0, pk / pmax_, 0) * w,
          PY(c(min(zk), zk, max(zk))),
          col = col_leaf, border = col_leaf_d, lwd = 0.9)
  ground(xs[4])
  text(xs[4] + 0.10 * w, PY(0.15), "LAI = the shaded area", cex = 0.50,
       font = 2, adj = c(0, 0.5), col = "#3F4C34")
  text(xs[4] + w / 2, PY(0) - 0.036, "leaf area per unit height", cex = 0.48,
       col = col_arrow)

  # ---- arrows between the steps ---------------------------------------------
  for (i in 1:3) {
    xa <- xs[i] + w + 0.008; xb <- xs[i + 1] - 0.008
    arrows(xa, PY(0.50), xb, PY(0.50), length = 0.045, lwd = 1.2,
           col = col_arrow)
  }

  # ---- the statement the figure exists to make ------------------------------
  text(0.5, 0.088,
       "every step is a modelling choice: the result is an estimate under stated assumptions, not a measurement",
       cex = 0.545, font = 3, col = "#5A5A5A")
}

intro_save(draw_figure, "Fig_intro_als_chain", width = 5.9, height = 3.15)

# ---- sanity checks ----------------------------------------------------------
# 1. the cumulative inversion is exact: reading the gap fraction back through
#    the Beer-Lambert law returns the plant area above each height.
pai_back <- -log(gap) / k_ext
stopifnot(max(abs(pai_back - cum)) < 1e-10)

# 2. differentiating that cumulative returns the profile drawn in panel (d),
#    to within the discretisation of the 400-point height grid.
dz  <- z[2] - z[1]
der <- rep(NA_real_, length(z))
der[2:(length(z) - 1)] <- -(pai_back[3:length(z)] - pai_back[1:(length(z) - 2)]) /
  (2 * dz)
ok  <- which(z > 0.15 & z < 0.95)
rel <- max(abs(der[ok] - pad[ok])) / max(pad)
stopifnot(rel < 0.01)

# 3. gap fraction increases with height and stays within [0, 1].
stopifnot(all(diff(gap) >= -1e-12), all(gap >= 0), all(gap <= 1))

message(sprintf(
  "check: cumulative inversion exact to %.1e; profile recovered to %.2f%% of peak; gap fraction in [%.2f, %.2f]",
  max(abs(pai_back - cum)), 100 * rel, min(gap), max(gap)))
