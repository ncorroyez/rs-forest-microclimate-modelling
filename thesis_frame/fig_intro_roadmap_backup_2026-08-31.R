# Figure 1.4 — conceptual roadmap for the General Introduction.
# Maps the three scientific gaps (Section 1.5) onto the three research
# questions (Section 1.6) and the three research chapters (Section 1.7).
# Base R only, no dependencies. Outputs a vector PDF and a raster PNG.

out_dir <- here::here("thesis_frame", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- palette: three neutral greys, one accent per column ------------------
col_gap    <- "#E8E4DC"
col_rq     <- "#DCE3E8"
col_chap   <- "#E2E8DC"
col_border <- "#7A7A7A"
col_arrow  <- "#5A5A5A"

# ---- helper: draw a box with wrapped, top-aligned text -------------------
box_text <- function(x0, y0, x1, y1, title, body, fill,
                     cex_title = 0.62, cex_body = 0.55, wrap = 34) {
  rect(x0, y0, x1, y1, col = fill, border = col_border, lwd = 0.8)
  pad <- 0.012
  lines_body <- strwrap(body, width = wrap)
  n <- length(lines_body) + 1L
  # vertical layout: title then body, block centered in the box
  lh <- 0.030
  top <- (y0 + y1) / 2 + (n * lh) / 2 - lh / 2
  text(x0 + pad, top, title, adj = c(0, 0.5), cex = cex_title, font = 2)
  for (i in seq_along(lines_body)) {
    text(x0 + pad, top - i * lh, lines_body[i], adj = c(0, 0.5), cex = cex_body)
  }
}

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # column geometry
  cx <- list(gap  = c(0.035, 0.315),
             rq   = c(0.360, 0.640),
             chap = c(0.685, 0.965))
  # row geometry (top of each row, bottom of each row)
  ry <- list(r1 = c(0.660, 0.880),
             r2 = c(0.395, 0.615),
             r3 = c(0.130, 0.350))

  # ---- column headers ----
  hdr <- 0.925
  text(mean(cx$gap),  hdr, "Scientific gap (1.5)",      cex = 0.68, font = 2)
  text(mean(cx$rq),   hdr, "Research question (1.6)",   cex = 0.68, font = 2)
  text(mean(cx$chap), hdr, "Chapter and approach (1.7)", cex = 0.68, font = 2)

  gaps <- list(
    c("1.5.1  Attribution",
      "Canopy structural variables covary in real stands, so a correlative model cannot assign buffering to any one of them."),
    c("1.5.2  Inter-sensor consistency",
      "LiDAR and Sentinel-2 both report a leaf area index, but they integrate different parts of the canopy column."),
    c("1.5.3  Spatio-temporal trade-off",
      "No single sensor delivers three-dimensional structure and temporal dynamics at the same time.")
  )
  rqs <- list(
    c("RQ1",
      "Which dimension of canopy structure governs sub-canopy ΔTmax buffering: leaf quantity, or vertical arrangement?"),
    c("RQ2",
      "Can optical LAI be reconciled with LiDAR-derived structural LAI, or does optical saturation decouple them?"),
    c("RQ3",
      "Which remotely sensed structure should force the model, and in which canopy regime?")
  )
  chaps <- list(
    c("Chapter 1",
      "Controlled trait perturbation in MuSICA around each plot's own measured canopy; validated against 53 loggers.  ALS only, plot scale."),
    c("Chapter 2",
      "Sequential correction: effective optical depth, then PROSAIL priors, then horizontal heterogeneity.  ALS + S2, three sites, 10 m grid."),
    c("Chapter 3",
      "Competing LAI scenarios forced into MuSICA, adjudicated on four independent axes, stratified by canopy density.")
  )
  fills <- list(col_gap, col_rq, col_chap)

  for (i in 1:3) {
    r <- ry[[i]]
    box_text(cx$gap[1],  r[1], cx$gap[2],  r[2], gaps[[i]][1],  gaps[[i]][2],  col_gap)
    box_text(cx$rq[1],   r[1], cx$rq[2],   r[2], rqs[[i]][1],   rqs[[i]][2],   col_rq)
    box_text(cx$chap[1], r[1], cx$chap[2], r[2], chaps[[i]][1], chaps[[i]][2], col_chap)

    ym <- mean(r)
    arrows(cx$gap[2] + 0.006, ym, cx$rq[1] - 0.006, ym,
           length = 0.055, lwd = 1.1, col = col_arrow)
    arrows(cx$rq[2] + 0.006, ym, cx$chap[1] - 0.006, ym,
           length = 0.055, lwd = 1.1, col = col_arrow)
  }

  # ---- vertical dependency arrows, in the white space between rows ----
  for (i in 1:2) {
    arrows(mean(cx$gap), ry[[i]][1] - 0.006, mean(cx$gap), ry[[i + 1]][2] + 0.006,
           length = 0.05, lwd = 1.1, col = col_arrow)
  }

  # ---- the ordering argument, along the left edge ----
  text(0.016, mean(c(ry$r1[2], ry$r3[1])),
       "the gaps must be closed in this order",
       srt = 90, cex = 0.55, col = col_arrow)

  # ---- red thread, bottom band ----
  rect(0.035, 0.020, 0.965, 0.098, col = "#F4F2EE", border = col_border, lwd = 0.8)
  txt <- paste("Until we know which structural dimension the microclimate responds to, we cannot say what a remote sensor needs to measure;",
               "until we know what part of the canopy the optical sensor sees, we cannot say whether it measures that thing.")
  ln <- strwrap(txt, width = 118)
  for (i in seq_along(ln)) {
    text(0.5, 0.075 - (i - 1) * 0.028, ln[i], cex = 0.56, font = 3)
  }
}

W <- 7.4; H <- 5.0

cairo_pdf(file.path(out_dir, "Fig_intro_roadmap.pdf"), width = W, height = H)
draw_figure()
dev.off()

png(file.path(out_dir, "Fig_intro_roadmap.png"), width = W, height = H,
    units = "in", res = 300, type = "cairo")
draw_figure()
dev.off()

message("written: ", file.path(out_dir, "Fig_intro_roadmap.pdf"))
