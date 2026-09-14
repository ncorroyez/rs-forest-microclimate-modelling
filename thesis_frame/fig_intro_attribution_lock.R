# Figure — the attribution lock, and the design that opens it.
#
# Illustrates Section 1.5.1. Panel (a): in real stands the structural variables
# covary, so the observations occupy a narrow diagonal band. A regression fitted
# to that band estimates one direction only; the combinations off the band do
# not exist in the field, so no amount of extra sampling along the same gradient
# separates the two effects. Panel (b): the perturbation design moves ONE
# variable at a time, in small steps, around EACH stand's own measured canopy.
# Because the moves are axis-aligned and the cloud is diagonal, they reach
# combinations the field does not offer (checked at the end of this script).
#
# Conceptual figure. No measured values; no numeric tick labels.

source(here::here("thesis_frame", "_intro_fig_style.R"))

# ---- a collinear cloud of "real stands" -------------------------------------
set.seed(11)
n       <- 46
tt      <- seq(0.10, 0.90, length.out = n) + runif(n, -0.018, 0.018)
sd_perp <- 0.026                                # narrow: the variables covary
perp    <- rnorm(n, 0, sd_perp)

ang <- 42 * pi / 180
u   <- c(cos(ang), sin(ang))                    # along the band
v   <- c(-sin(ang), cos(ang))                   # across the band

ctr  <- c(0.50, 0.50)
half <- 0.33                                    # half-length of the band
sx   <- ctr[1] + (tt - 0.5) * 2 * half * u[1] + perp * v[1]
sy   <- ctr[2] + (tt - 0.5) * 2 * half * u[2] + perp * v[2]

# ---- the perturbation design ------------------------------------------------
# four stands spread along the band, each perturbed one variable at a time
ai   <- vapply(c(0.15, 0.40, 0.65, 0.88),
               function(a) which.min(abs(tt - a)), integer(1))
step <- 0.105                                   # one small native step

cross_pts <- function(x, y)
  data.frame(x = c(x - step, x + step, x, x),
             y = c(y, y, y - step, y + step))

col_band  <- "#EFECE6"   # the strip of combinations the field does offer
col_never <- "#DFDAD1"   # the rest of the plane: no stand is ever found there
col_fit   <- "#8A5F5F"

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))

  w <- 0.360; h <- 0.520; y0 <- 0.235
  x0a <- 0.105; x0b <- 0.560
  P <- function(x0) function(x, y) cbind(x0 + x * w, y0 + y * h)

  one_panel <- function(x0, title, subtitle) {
    rect(x0, y0, x0 + w, y0 + h, col = col_never, border = col_border,
         lwd = 0.7)
    text(x0 + w / 2, y0 + h + 0.098, title, cex = 0.66, font = 2)
    ln <- strwrap(subtitle, width = 40)
    for (i in seq_along(ln))
      text(x0 + w / 2, y0 + h + 0.056 - (i - 1) * 0.030, ln[i], cex = 0.545,
           col = col_arrow)
    # the strip the field actually samples, clipped to the panel
    Pn <- P(x0); L <- 1.4; b <- 2.3 * sd_perp
    cx <- ctr[1] + c(-L, L) * u[1]; cy <- ctr[2] + c(-L, L) * u[2]
    px <- c(cx + b * v[1], rev(cx) - b * v[1])
    py <- c(cy + b * v[2], rev(cy) - b * v[2])
    pp <- Pn(px, py)
    clip(x0, x0 + w, y0, y0 + h)
    polygon(pp[, 1], pp[, 2], col = col_band, border = NA)
    clip(0, 1, 0, 1)                            # release the clip region
    Pn
  }

  # ---- (a) what the field offers --------------------------------------------
  pa <- one_panel(x0a, "(a) what real stands offer",
                  "the two variables covary, so the observations occupy one direction")

  ends <- pa(ctr[1] + c(-half, half) * u[1], ctr[2] + c(-half, half) * u[2])
  segments(ends[1, 1], ends[1, 2], ends[2, 1], ends[2, 2], col = col_fit,
           lwd = 1.6)
  pt <- pa(sx, sy)
  points(pt[, 1], pt[, 2], pch = 21, bg = col_leaf, col = col_leaf_d,
         cex = 0.60, lwd = 0.6)

  d  <- 0.135
  a1 <- pa(ctr[1] - d * v[1], ctr[2] - d * v[2])
  a2 <- pa(ctr[1] + d * v[1], ctr[2] + d * v[2])
  arrows(a1[1, 1], a1[1, 2], a2[1, 1], a2[1, 2], length = 0.032, code = 3,
         lwd = 1.1, lty = 2, col = "#6E6E6E")
  text(ends[2, 1] + 0.006, ends[2, 2] - 0.020, "estimated", cex = 0.50,
       adj = c(0, 1), col = col_fit)
  text(a2[1, 1] - 0.010, a2[1, 2] + 0.004, "not estimated", cex = 0.50,
       adj = c(1, 0), col = "#6E6E6E")
  text(x0a + 0.80 * w, y0 + 0.10 * h, "no stand here", cex = 0.50, font = 3,
       adj = c(0.5, 0.5), col = "#7C766C")

  # ---- (b) what the design generates ---------------------------------------
  pb <- one_panel(x0b, "(b) what the perturbation design generates",
                  "one variable moved at a time, around each stand's own canopy")

  pt <- pb(sx, sy)
  points(pt[, 1], pt[, 2], pch = 21, bg = "#CFCCC4", col = "#AFACA4",
         cex = 0.48, lwd = 0.4)
  for (k in ai) {
    cp <- cross_pts(sx[k], sy[k])
    o  <- pb(sx[k], sy[k]); e <- pb(cp$x, cp$y)
    segments(o[1, 1], o[1, 2], e[, 1], e[, 2], col = col_als, lwd = 1)
    points(e[, 1], e[, 2], pch = 21, bg = col_als, col = "#5D778E",
           cex = 0.52, lwd = 0.5)
    points(o[1, 1], o[1, 2], pch = 21, bg = col_leaf, col = col_leaf_d,
           cex = 0.70, lwd = 0.7)
  }

  # ---- shared axes ----------------------------------------------------------
  for (x0 in c(x0a, x0b)) {
    arrows(x0, y0 - 0.014, x0 + w, y0 - 0.014, length = 0.036, lwd = 1,
           col = col_arrow)
    arrows(x0 - 0.014, y0, x0 - 0.014, y0 + h, length = 0.036, lwd = 1,
           col = col_arrow)
    text(x0 + w / 2, y0 - 0.052, "leaf quantity", cex = 0.55, col = col_arrow)
    text(x0 - 0.040, y0 + h / 2, "vertical arrangement", srt = 90, cex = 0.55,
         col = col_arrow)
  }

  # ---- legend ---------------------------------------------------------------
  ly <- 0.115
  points(0.108, ly, pch = 21, bg = col_leaf, col = col_leaf_d, cex = 0.7,
         lwd = 0.7)
  text(0.126, ly, "a measured stand", cex = 0.545, adj = c(0, 0.5))
  points(0.350, ly, pch = 21, bg = col_als, col = "#5D778E", cex = 0.7,
         lwd = 0.5)
  text(0.368, ly, "a simulated perturbation", cex = 0.545, adj = c(0, 0.5))
  rect(0.645, ly - 0.015, 0.668, ly + 0.015, col = col_band,
       border = col_border, lwd = 0.5)
  text(0.678, ly, "what the field offers", cex = 0.545, adj = c(0, 0.5))
  rect(0.845, ly - 0.015, 0.868, ly + 0.015, col = col_never,
       border = col_border, lwd = 0.5)
  text(0.878, ly, "what it does not", cex = 0.545, adj = c(0, 0.5))
}

intro_save(draw_figure, "Fig_intro_attribution_lock", width = 6.0, height = 3.35)

# ---- sanity checks ----------------------------------------------------------
# 1. the cloud really is collinear: its spread across the band is a small
#    fraction of its spread along it.
along  <- (sx - ctr[1]) * u[1] + (sy - ctr[2]) * u[2]
across <- (sx - ctr[1]) * v[1] + (sy - ctr[2]) * v[2]
stopifnot(sd(across) < 0.20 * sd(along))

# 2. the axis-aligned perturbations reach outside the band the stands occupy,
#    which is the whole claim of panel (b).
obs_half   <- max(abs(across))
tips       <- do.call(rbind, lapply(ai, function(k) cross_pts(sx[k], sy[k])))
tip_across <- (tips$x - ctr[1]) * v[1] + (tips$y - ctr[2]) * v[2]
stopifnot(all(tapply(abs(tip_across) > obs_half,
                     rep(seq_along(ai), each = 4), any)))

message(sprintf(
  "check: spread across/along = %.2f; %d of %d perturbation points outside the observed band",
  sd(across) / sd(along), sum(abs(tip_across) > obs_half), nrow(tips)))
