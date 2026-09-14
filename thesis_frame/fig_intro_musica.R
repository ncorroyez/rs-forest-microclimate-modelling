# Figure — why the canopy scheme has to be multilayer.
#
# Illustrates Section 1.3.2. The argument there is not what MuSICA computes but
# what it can ACCEPT: a vertical profile of leaf area density is a native input,
# so the arrangement of foliage can be contrasted against its quantity. A
# big-leaf scheme collapses the column into one effective surface and can only
# take a scalar, so the question of Section 1.2.2 cannot even be posed to it.
#
# Deliberately NOT a model-internals diagram, and deliberately not the Chapter 1
# pipeline figure: no sampling design, no typology, no results.
#
# Conceptual figure. No measured values anywhere.

source(here::here("thesis_frame", "_intro_fig_style.R"))

z <- seq(0, 1, length.out = 300)
pad <- exp(-((z - 0.60)^2) / (2 * 0.21^2)) * (z > 0.10)
pad <- pad / max(pad)

# canopy column geometry, shared by both panels
CY0 <- 0.335; CY1 <- 0.735          # ground to canopy top
SY0 <- 0.245; SY1 <- CY0            # soil block

col_air  <- "#DCE3E8"
col_soil <- "#E3DCD2"
col_sun  <- "#E8C77A"
col_shd  <- "#9AA88F"

# ---- panel (a): multilayer -------------------------------------------------
panel_multi <- function(x0) {
  px <- function(u) x0 + u * 0.085          # profile occupies [x0, x0+0.085]
  cx0 <- x0 + 0.125; cx1 <- cx0 + 0.155     # the canopy column

  # the LiDAR-derived profile, arriving as forcing
  k <- pad > 0.01
  polygon(c(px(0), px(pad[k]), px(0)),
          c(CY0 + min(z[k]) * (CY1 - CY0), CY0 + z[k] * (CY1 - CY0),
            CY0 + max(z[k]) * (CY1 - CY0)),
          col = col_leaf, border = col_leaf_d, lwd = 0.9)
  text(px(0.5), CY1 + 0.038, "leaf area density", cex = 0.50, col = col_arrow)
  text(px(0.5), CY1 + 0.012, "profile (from ALS)", cex = 0.50, col = col_arrow)
  arrows(px(1.10), (CY0 + CY1) / 2, cx0 - 0.008, (CY0 + CY1) / 2,
         length = 0.055, lwd = 1.3, col = col_leaf_d)

  # the discretised canopy: each layer carries sunlit and shaded foliage
  n <- 8
  ys <- seq(CY0, CY1, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(cx0, ys[i], cx1, ys[i + 1], col = "#FFFFFF", border = col_border,
         lwd = 0.6)
    # local leaf area sets how much foliage the layer holds
    u <- approx(z, pad, xout = (mean(ys[i:(i + 1)]) - CY0) / (CY1 - CY0))$y
    w <- (cx1 - cx0) * 0.40 * u
    rect(cx0 + 0.006, ys[i] + 0.004, cx0 + 0.006 + w, ys[i + 1] - 0.004,
         col = col_sun, border = NA)
    rect(cx1 - 0.006 - w, ys[i] + 0.004, cx1 - 0.006, ys[i + 1] - 0.004,
         col = col_shd, border = NA)
  }
  rect(cx0, SY0, cx1, SY1, col = col_soil, border = col_border, lwd = 0.7)
  text((cx0 + cx1) / 2, (SY0 + SY1) / 2, "soil", cex = 0.50, col = col_ground)

  list(cx0 = cx0, cx1 = cx1)
}

# ---- panel (b): big-leaf ---------------------------------------------------
panel_big <- function(x0) {
  px <- function(u) x0 + u * 0.085
  cx0 <- x0 + 0.125; cx1 <- cx0 + 0.155

  # the profile is not an admissible input: it is reduced to one number
  rect(px(0.05), (CY0 + CY1) / 2 - 0.030, px(0.95), (CY0 + CY1) / 2 + 0.030,
       col = "#EDEAE4", border = col_border, lwd = 0.8)
  text(px(0.5), (CY0 + CY1) / 2, "LAI", cex = 0.62, font = 2)
  text(px(0.5), CY1 + 0.038, "one scalar", cex = 0.50, col = col_arrow)
  text(px(0.5), CY1 + 0.012, "(the profile is lost)", cex = 0.50,
       col = "#9C5F5F", font = 3)
  arrows(px(1.10), (CY0 + CY1) / 2, cx0 - 0.008, (CY0 + CY1) / 2,
         length = 0.055, lwd = 1.3, col = col_border)

  # one effective surface
  rect(cx0, CY0, cx1, CY1, col = "#FBFBFA", border = col_border, lwd = 0.6,
       lty = 3)
  yb <- CY0 + 0.62 * (CY1 - CY0)
  rect(cx0, yb - 0.028, cx1, yb + 0.028, col = col_leaf, border = col_leaf_d,
       lwd = 0.9)
  text((cx0 + cx1) / 2, yb, "one effective surface", cex = 0.485)
  rect(cx0, SY0, cx1, SY1, col = col_soil, border = col_border, lwd = 0.7)
  text((cx0 + cx1) / 2, (SY0 + SY1) / 2, "soil", cex = 0.50, col = col_ground)

  list(cx0 = cx0, cx1 = cx1)
}

# ---- shared decoration: forcing above, sub-canopy air below ----------------
decorate <- function(g, note, note_col) {
  cxm <- (g$cx0 + g$cx1) / 2
  rect(g$cx0, CY1 + 0.055, g$cx1, CY1 + 0.100, col = col_air,
       border = col_border, lwd = 0.7)
  text(cxm, CY1 + 0.0775, "above-canopy forcing", cex = 0.485)
  arrows(cxm, CY1 + 0.050, cxm, CY1 + 0.008, length = 0.045, lwd = 1.1,
         col = col_arrow)
  # what comes out, at logger height
  arrows(g$cx1 + 0.012, CY0 + 0.045, g$cx1 + 0.048, CY0 + 0.045,
         length = 0.045, lwd = 1.1, col = col_arrow)
  text(g$cx1 + 0.052, CY0 + 0.045, expression(Delta*T[max]), adj = c(0, 0.5),
       cex = 0.56)
  text(cxm, SY0 - 0.038, note, cex = 0.495, font = 3, col = note_col)
}

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.138, 1))

  ga <- panel_multi(0.075)
  gb <- panel_big(0.565)
  decorate(ga, "the arrangement can be changed at fixed leaf area", col_arrow)
  decorate(gb, "only the quantity can be changed", "#9C5F5F")

  text(0.075 + 0.125 + 0.078, 0.930, "(a) multilayer scheme", cex = 0.64,
       font = 2)
  text(0.075 + 0.125 + 0.078, 0.892, "profile is a native input", cex = 0.535,
       col = col_arrow)
  text(0.565 + 0.125 + 0.078, 0.930, "(b) big-leaf scheme", cex = 0.64,
       font = 2)
  text(0.565 + 0.125 + 0.078, 0.892, "profile cannot be accepted", cex = 0.535,
       col = col_arrow)

  # legend for the two leaf classes
  ly <- 0.168
  rect(0.205, ly - 0.013, 0.227, ly + 0.013, col = col_sun, border = NA)
  text(0.235, ly, "sunlit foliage", adj = c(0, 0.5), cex = 0.495)
  rect(0.395, ly - 0.013, 0.417, ly + 0.013, col = col_shd, border = NA)
  text(0.425, ly, "shaded foliage", adj = c(0, 0.5), cex = 0.495)
  text(0.575, ly, "resolved layer by layer in (a) only", adj = c(0, 0.5),
       cex = 0.495, col = col_arrow)
}

intro_save(draw_figure, "Fig_intro_musica", width = 6.6, height = 3.8)
