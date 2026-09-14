# Figure — how a Sentinel-2 leaf area index is made, and what limits it.
#
# Illustrates Section 1.4.2. (a) what the satellite records: one reflectance
# value per band, at two native resolutions, from above. (b) how that becomes a
# leaf area index: a radiative-transfer model is run over prescribed parameter
# ranges to build a simulated library, and a regressor trained on that library
# is applied to the observation. The prior ranges are the regularisation, and
# they are tuned primarily for crops. (c) why the inversion is ill-posed: two
# different canopies produce nearly the same spectrum. (d) why it saturates: the
# reflectance stops responding once leaf area has accumulated.
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

# ---- a schematic vegetation reflectance spectrum ----------------------------
lam <- seq(400, 2400, length.out = 600)
gauss <- function(x, m, s) exp(-((x - m)^2) / (2 * s^2))
spectrum <- function(nir = 0.46, cab = 1) {
  vis <- 0.030 + 0.055 * gauss(lam, 550, 40) / cab
  edge <- 1 / (1 + exp(-(lam - 715) / 16))
  swir <- 1 - 0.42 * gauss(lam, 1450, 95) - 0.60 * gauss(lam, 1930, 110) -
    0.30 * pmax(0, (lam - 1100) / 1300)
  r <- vis * (1 - edge) + nir * edge * swir
  pmax(r, 0.012)
}
ref <- spectrum()

# the two Sentinel-2 resolutions, and where the bands sit
bands <- data.frame(
  wl  = c(492, 560, 665, 704, 740, 783, 833, 865, 1610, 2190),
  res = c(10, 10, 10, 20, 20, 20, 10, 20, 20, 20))

# ---- the ill-posedness pair -------------------------------------------------
# more leaf area with less pigment, against less leaf area with more pigment
ref_a <- spectrum(nir = 0.470, cab = 0.86)
ref_b <- spectrum(nir = 0.455, cab = 1.18)

# ---- the saturating response ------------------------------------------------
lai_ax <- seq(0, 1, length.out = 200)
sat    <- 1 - exp(-3.6 * lai_ax)

col_lut <- "#C9C4BA"
col_ink <- "#3A3A3A"

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))

  y0 <- 0.300; h <- 0.430; w <- 0.196
  xs <- c(0.045, 0.283, 0.521, 0.759)
  PY <- function(v) y0 + v * h
  PX <- function(x0, v) x0 + v * w
  L  <- function(v) (v - 400) / 2000                # wavelength to [0, 1]

  head <- function(x0, tag, title) {
    text(x0 + w / 2, y0 + h + 0.112, tag, cex = 0.62, font = 2)
    ln <- strwrap(title, width = 27)
    for (i in seq_along(ln))
      text(x0 + w / 2, y0 + h + 0.070 - (i - 1) * 0.030, ln[i], cex = 0.545,
           col = col_arrow)
  }
  axes <- function(x0, xlab, ylab) {
    arrows(x0, PY(0), x0 + w, PY(0), length = 0.032, lwd = 1, col = col_arrow)
    arrows(x0, PY(0), x0, PY(1.02), length = 0.032, lwd = 1, col = col_arrow)
    text(x0 + w / 2, PY(0) - 0.040, xlab, cex = 0.475, col = col_arrow)
    text(x0 - 0.020, PY(0.5), ylab, srt = 90, cex = 0.475, col = col_arrow)
  }

  # ---- (a) what the satellite records ---------------------------------------
  head(xs[1], "(a)", "what the satellite records")
  ry <- function(r) PY(r / 0.55)
  lines(PX(xs[1], L(lam)), ry(ref), col = col_s2, lwd = 1.4)
  for (i in seq_len(nrow(bands))) {
    bx <- PX(xs[1], L(bands$wl[i])); by <- ry(spectrum()[which.min(abs(lam - bands$wl[i]))])
    points(bx, by, pch = if (bands$res[i] == 10) 19 else 21,
           bg = "white", col = "#8A6E43", cex = 0.45, lwd = 0.7)
  }
  axes(xs[1], "wavelength", "reflectance")
  points(PX(xs[1], 0.06), PY(0.93), pch = 19, col = "#8A6E43", cex = 0.45)
  text(PX(xs[1], 0.12), PY(0.93), "10 m band", cex = 0.455, adj = c(0, 0.5))
  points(PX(xs[1], 0.06), PY(0.80), pch = 21, bg = "white", col = "#8A6E43",
         cex = 0.45, lwd = 0.7)
  text(PX(xs[1], 0.12), PY(0.80), "20 m band", cex = 0.455, adj = c(0, 0.5))
  text(xs[1] + w / 2, PY(0) - 0.070, "one value per band, seen from above",
       cex = 0.455, font = 3, col = col_arrow)

  # ---- (b) the simulated library --------------------------------------------
  head(xs[2], "(b)", "inversion against a simulated library")
  set.seed(3)
  for (j in 1:26) {
    rj <- spectrum(nir = runif(1, 0.30, 0.56), cab = runif(1, 0.75, 1.35))
    lines(PX(xs[2], L(lam)), PY(rj / 0.55), col = col_lut, lwd = 0.5)
  }
  lines(PX(xs[2], L(lam)), PY(ref / 0.55), col = col_s2, lwd = 1.4)
  axes(xs[2], "wavelength", "reflectance")
  text(xs[2] + w / 2, PY(0) - 0.070,
       "a regressor is trained on the library", cex = 0.455, font = 3,
       col = col_arrow)
  rect(xs[2] - 0.004, PY(0) - 0.176, xs[2] + w + 0.004, PY(0) - 0.090,
       col = "#F6EFE4", border = "#C2A97F", lwd = 0.7)
  wrap_text("the prescribed parameter ranges are the prior, and they are tuned for crops",
            32, xs[2] + 0.008, PY(0) - 0.106, cex = 0.455, lh = 0.026,
            col = "#7A5F30")

  # ---- (c) the inversion is ill-posed ---------------------------------------
  head(xs[3], "(c)", "different canopies, one spectrum")
  lines(PX(xs[3], L(lam)), PY(ref_a / 0.55), col = "#7E9BB5", lwd = 1.5)
  lines(PX(xs[3], L(lam)), PY(ref_b / 0.55), col = "#B08A5A", lwd = 1.5,
        lty = 2)
  axes(xs[3], "wavelength", "reflectance")
  segments(PX(xs[3], 0.06), PY(0.97), PX(xs[3], 0.13), PY(0.97),
           col = "#7E9BB5", lwd = 1.5)
  text(PX(xs[3], 0.16), PY(0.97), "more leaf area, less pigment", cex = 0.44,
       adj = c(0, 0.5))
  segments(PX(xs[3], 0.06), PY(0.86), PX(xs[3], 0.13), PY(0.86),
           col = "#B08A5A", lwd = 1.5, lty = 2)
  text(PX(xs[3], 0.16), PY(0.86), "less leaf area, more pigment", cex = 0.44,
       adj = c(0, 0.5))
  text(xs[3] + w / 2, PY(0) - 0.070,
       "so the retrieval needs a prior to choose", cex = 0.455, font = 3,
       col = col_arrow)

  # ---- (d) the response saturates -------------------------------------------
  head(xs[4], "(d)", "and the response saturates")
  # the regime closed temperate canopies occupy, drawn first so the curve sits on top
  rect(PX(xs[4], 0.60), PY(0), PX(xs[4], 1), PY(1.02), col = "#EDE7DD",
       border = NA)
  lines(PX(xs[4], lai_ax), PY(sat * 0.92), col = col_s2, lwd = 1.6)
  axes(xs[4], "leaf area index", "reflectance")
  text(PX(xs[4], 0.80), PY(0.40), "closed canopy", cex = 0.455,
       col = "#7C6A4E")
  text(PX(xs[4], 0.80), PY(0.31), "sits here", cex = 0.455, col = "#7C6A4E")
  text(xs[4] + w / 2, PY(0) - 0.070,
       "differences stop reaching the sensor", cex = 0.455, font = 3,
       col = col_arrow)

  # ---- one arrow for the retrieval, a divider before the two limits ---------
  xd <- (xs[2] + w + xs[3] - 0.03) / 2
  segments(xd, PY(-0.10), xd, PY(1.16), col = "#D6D2CA", lwd = 0.9)
  text((xs[1] + xs[2] + w) / 2, y0 + h + 0.198, "how the retrieval works",
       cex = 0.575, font = 2, col = "#6E6E6E")
  text((xs[3] + xs[4] + w) / 2, y0 + h + 0.198, "what limits it",
       cex = 0.575, font = 2, col = "#6E6E6E")

  text(0.5, 0.052,
       "wall-to-wall and repeated, but two-dimensional: the depth of the canopy is what the signal cannot reach",
       cex = 0.545, font = 3, col = "#5A5A5A")
}

intro_save(draw_figure, "Fig_intro_prosail_chain", width = 5.9, height = 3.5)

# ---- sanity checks ----------------------------------------------------------
# 1. the two "different canopies" really do produce near-identical spectra:
#    their maximum separation is a small fraction of the spectrum's own range.
sep <- max(abs(ref_a - ref_b)) / diff(range(ref))
stopifnot(sep < 0.10)

# 2. the response really saturates: its slope in the dense half is a small
#    fraction of its slope in the open half.
d   <- diff(sat) / diff(lai_ax)
lo  <- mean(d[lai_ax[-1] < 0.3]); hi <- mean(d[lai_ax[-1] > 0.7])
stopifnot(hi < 0.10 * lo)

message(sprintf(
  "check: the two spectra differ by %.1f%% of the spectral range; slope falls to %.1f%% of its open-canopy value",
  100 * sep, 100 * hi / lo))
