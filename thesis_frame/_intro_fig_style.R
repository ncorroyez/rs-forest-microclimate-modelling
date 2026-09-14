# Shared style for the General Introduction figures.
# Sourced by fig_intro_*.R. Base R only, no dependencies.
#
# All introduction figures are conceptual. They carry no measured values, and
# no axis in them is given numeric tick labels: the shapes illustrate an
# argument, they do not report a result.

intro_out_dir <- function() {
  d <- here::here("thesis_frame", "figures")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# ---- palette (matches Fig_I4_roadmap) --------------------------------------
col_border <- "#7A7A7A"
col_arrow  <- "#5A5A5A"
col_panel  <- "#F4F2EE"   # panel background
col_leaf   <- "#B9C6AE"   # foliage / leaf area density
col_leaf_d <- "#8FA382"   # foliage, darker edge
col_als    <- "#7E9BB5"   # airborne LiDAR
col_s2     <- "#C9A87C"   # Sentinel-2
col_gedi   <- "#A99BB5"   # GEDI
col_ground <- "#6B6257"

# ---- helper: soft text wrap at a given character width ---------------------
wrap_text <- function(txt, width, x, y, cex, lh, adj = c(0, 1), font = 1,
                      col = "black") {
  ln <- strwrap(txt, width = width)
  for (i in seq_along(ln)) {
    text(x, y - (i - 1) * lh, ln[i], adj = adj, cex = cex, font = font,
         col = col)
  }
  invisible(length(ln))
}

# ---- helper: write both a vector PDF and a 300 dpi PNG ---------------------
intro_save <- function(draw, name, width, height) {
  d <- intro_out_dir()
  cairo_pdf(file.path(d, paste0(name, ".pdf")), width = width, height = height)
  draw()
  dev.off()
  png(file.path(d, paste0(name, ".png")), width = width, height = height,
      units = "in", res = 300, type = "cairo")
  draw()
  dev.off()
  message("written: ", file.path(d, paste0(name, ".pdf")))
}
