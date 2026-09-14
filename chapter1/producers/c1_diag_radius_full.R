# Diagnostic (supervisor #5/#6 extended, 2026-09-11): the FULL footprint-radius
# sweep, everything as a function of radius, per archetype, as a shift from the
# 20 m reference. Two MuSICA responses (ΔTmax, buffering slope) and the four
# structural inputs the footprint moves (one-sided LAI, fCover, Hmax, LAD
# centroid). All 7 radii on disk (5,10,12.5,15,20,25,50 m); 20 m = reference.
# Pure post-processing, CHS41-Rmerge / no wind. ~400 cLHS plots per radius.
# In : radius_clhs_by_archetype.csv (ΔTmax, slope) ; radius_clhs_vars.csv (structure)
# Out: out_files/Chapter1/figures/diag_2026-09-10/diag_radius_full_vs_radius.png
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
REF <- 20
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#A6D96A", P4="#1A9850", ALL="grey30")
# --- MuSICA responses: already a shift vs the 20 m reference ---
M <- fread("out_files/Chapter1/tables/radius_clhs_by_archetype.csv")
mus <- rbind(
  M[, .(group, radius, quantity = "ΔT_max vs 20 m (°C)", val = d_dtmax, lo = dt_lo, hi = dt_hi)],
  M[, .(group, radius, quantity = "buffering slope vs 20 m (–)", val = d_slope, lo = sl_lo, hi = sl_hi)])
# --- structural inputs: per-plot shift vs the plot's own 20 m value, then median per archetype + ALL ---
V <- fread("out_files/Chapter1/tables/radius_clhs_vars.csv")
V0 <- V[radius == REF, .(id_plot, LAI0 = LAI, H0 = Hmax, f0 = fCover, c0 = lad_centroid_rel)]
D <- merge(V, V0, by = "id_plot")
D[, `:=`(dLAI = LAI - LAI0, dfCover = fCover - f0, dHmax = Hmax - H0, dcent = lad_centroid_rel - c0)]
shift <- function(col, label){
  a <- D[, .(val = median(get(col), na.rm = TRUE)), by = .(group = P, radius)]
  b <- D[, .(group = "ALL", val = median(get(col), na.rm = TRUE)), by = radius]
  x <- rbind(a, b); x[, `:=`(quantity = label, lo = NA_real_, hi = NA_real_)]; x
}
str <- rbindlist(list(
  shift("dLAI",    "one-sided LAI vs 20 m"),
  shift("dfCover", "fractional cover vs 20 m"),
  shift("dHmax",   "H_max vs 20 m (m)"),
  shift("dcent",   "LAD centroid (rel.) vs 20 m")))
A <- rbind(mus, str, use.names = TRUE)
qlev <- c("ΔT_max vs 20 m (°C)","buffering slope vs 20 m (–)","one-sided LAI vs 20 m",
          "fractional cover vs 20 m","H_max vs 20 m (m)","LAD centroid (rel.) vs 20 m")
A[, quantity := factor(quantity, levels = qlev)]
A[, group := factor(group, levels = c("P1","P2","P3","P4","ALL"))]
g <- ggplot(A, aes(radius, val, colour = group, fill = group)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey65") +
  geom_vline(xintercept = REF, linetype = "dotted", colour = "grey65") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .10, colour = NA) +
  geom_line(linewidth = .55) + geom_point(size = 1.3) +
  scale_colour_manual(values = PAL, name = NULL) + scale_fill_manual(values = PAL, guide = "none") +
  facet_wrap(~quantity, scales = "free_y", ncol = 2) +
  labs(x = "footprint radius (m)", y = "shift from the 20 m reference",
       subtitle = "All 7 radii, 400 cLHS plots. Top row = MuSICA responses (with 95% CI); rows below = the structural inputs the footprint moves.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold"))
ggsave("out_files/Chapter1/figures/diag_2026-09-10/diag_radius_full_vs_radius.png",
       g, width = 9.5, height = 8, dpi = 300, bg = "white")
# console: pooled median |shift| vs 20 m, per quantity, at the extreme radii
cat("=== pooled (ALL) shift vs 20 m at 5 m and 50 m ===\n")
print(A[group == "ALL" & radius %in% c(5,50), .(quantity, radius, val = round(val,3))][order(quantity, radius)])
cat("\nradii covered:", paste(sort(unique(A$radius)), collapse = ", "), "m (20 m = reference)\n")
cat("DONE -> diag_radius_full_vs_radius.png\n")
