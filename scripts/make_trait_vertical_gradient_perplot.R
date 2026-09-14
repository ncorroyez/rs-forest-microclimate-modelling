# ==============================================================================
# Within-canopy VERTICAL GRADIENT response to each trait, WITHOUT a common baseline.
#
# WHY THIS EXISTS. The shipped Figure B3 (scripts/make_trait_vertical_gradient.R)
# builds every increment as a step from ONE common canopy, the 400-plot design mean
# (one-sided LAI 3.42, Hmax 23 m). That makes each increment's size depend on how far
# an archetype's measured trait sits from that baseline, not on how sensitive the
# archetype is: the profile appears to lead in P2 and P3 largely because their measured
# leaf areas (3.73, 3.67) sit almost ON the baseline, leaving the LAI step near zero,
# whereas the same step strips two thirds of the canopy in P1. It is also computed on
# the 53 logger pixels only (n = 8, 12, 13, 20 per archetype).
#
# This script rebuilds the same four panels on the ON-MANIFOLD design instead. Each of
# the 400 cLHS plots is perturbed around ITS OWN measured canopy by the fixed native
# steps of Section 2.6, exactly as Figure 4 does, so:
#   - there is no baseline and no distance-from-baseline artifact;
#   - n rises from 53 to 400 plots (100 per archetype);
#   - the vertical panels become directly commensurable with Figure 4.
#
# Pure re-extraction from the archived design NetCDFs. It does NOT run MuSICA.
#
# Reads : out_files/Chapter1/nc_sensitivity_perplot_units/S####_{base,LAIp,LAIm,
#         fCovp,fCovm,Hmaxp,Hmaxm,unifLAD}.nc     (400 plots x 8 configurations)
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds  (archetypes)
# Writes: outputs/figures_chap1/fig_trait_vertical_gradient_perplot.{png,pdf}
#         outputs/figures_chap1/tab_trait_vertical_gradient_perplot.csv
#   Rscript scripts/make_trait_vertical_gradient_perplot.R [n_plots_per_archetype]
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(data.table); library(ggplot2); library(patchwork)
  library(lubridate); library(parallel)
})
setwd(here::here())
source("scripts/_article_style.R"); source("R/cluster_relabel.R")

NC   <- "out_files/Chapter1/nc_sensitivity_perplot_units"
OUT  <- "outputs/figures_chap1"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MONTHS <- 6:9; HOURS <- 10:16          # same window as Fig. B2 / B3
P_HPA  <- 1013
args <- commandArgs(trailingOnly = TRUE)
NPER <- if (length(args)) as.integer(args[1]) else 100L   # plots per archetype (100 = all)

S <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
S[, pid := sprintf("S%04d", .I)][, P := factor(relabel_cluster(Cluster), levels = c("P1","P2","P3","P4"))]
set.seed(1); S <- S[, .SD[seq_len(min(.N, NPER))], by = P]
cat(sprintf("plots: %d (%s)\n", nrow(S), paste(S[, .N, by = P][order(P)]$N, collapse = "/")))

#' Tetens saturation vapour pressure (hPa)
#' @param Tc air temperature, degrees C
#' @return saturation vapour pressure, hPa
esat <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))   # Tetens, hPa

# midday-summer mean vertical profile of the four variables, one run
#' Midday-summer mean vertical profile of T, wind, RH and VPD for one design run
#' @param f path to a design NetCDF
#' @return data.table(lev, rh_rel, Tair, Wind, RH, VPD), or NULL if missing or truncated
prof <- function(f) {
  if (!file.exists(f) || file.size(f) < 1e5) return(NULL)
  nc <- try(nc_open(f), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value
  tv <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC") + ncvar_get(nc, "time") * 3600
  k  <- month(tv) %in% MONTHS & hour(tv) %in% HOURS
  if (!any(k)) return(NULL)
  Tk <- ncvar_get(nc, "Tair_z")[, k, drop = FALSE] - 273.15
  q  <- ncvar_get(nc, "wair_z")[, k, drop = FALSE]
  wd <- ncvar_get(nc, "wind_z")[, k, drop = FALSE]
  e  <- q * P_HPA / (0.622 + q)                       # vapor pressure, hPa
  es <- esat(Tk)
  rh <- 100 * e / es; rh[rh > 100] <- 100      # pmin() drops the matrix dims; keep them
  data.table(lev = seq_len(nrow(Tk)),
             rh_rel = ncvar_get(nc, "relative_height"),
             Tair = rowMeans(Tk), Wind = rowMeans(wd),
             RH = rowMeans(rh), VPD = rowMeans((es - e) / 10))   # kPa
}

STEPS <- c(dLAI = "LAIp", dHmax = "Hmaxp", dfCover = "fCovp", dLAD = "unifLAD")
#' Per-trait profile increments for one plot, base minus or plus the fixed native step
#' @param i row index in the sampled design table S
#' @return long data.table(pid, P, step, z, Tair, Wind, RH, VPD), or NULL if base unreadable
one <- function(i) {
  pid <- S$pid[i]
  b <- prof(file.path(NC, paste0(pid, "_base.nc"))); if (is.null(b)) return(NULL)
  rbindlist(lapply(names(STEPS), function(st) {
    p <- prof(file.path(NC, sprintf("%s_%s.nc", pid, STEPS[[st]])))
    if (is.null(p) || nrow(p) != nrow(b)) return(NULL)
    # SIGN CONVENTION: for the scalar traits the increment is (perturbed - base), the
    # response to ADDING the step. For the profile the chapter's contrast is
    # real-minus-uniform, so the uniform run is subtracted FROM the base.
    d <- if (st == "dLAD") b[, .(Tair, Wind, RH, VPD)] - p[, .(Tair, Wind, RH, VPD)]
         else               p[, .(Tair, Wind, RH, VPD)] - b[, .(Tair, Wind, RH, VPD)]
    cbind(data.table(pid = pid, P = S$P[i], step = st, z = b$rh_rel), d)
  }))
}
R <- rbindlist(mclapply(seq_len(nrow(S)), function(i) tryCatch(one(i), error = function(e) NULL),
                        mc.cores = max(1L, detectCores() - 2L)), fill = TRUE)
cat(sprintf("plots with a complete set: %d\n", uniqueN(R$pid)))

M <- melt(R, id.vars = c("pid","P","step","z"), variable.name = "var", value.name = "d")
# Each plot has its OWN relative-height vector, so averaging by raw z would mix a
# different set of plots at every level and produce spaghetti above the canopy top.
# Interpolate every plot onto one common grid first, then average.
ZG <- seq(0, 1.4, by = 0.05)   # capped: above 1.45 the per-plot grids stop covering, n collapses to 1
M <- M[, {
  o <- order(z)
  .(z = ZG, d = approx(z[o], d[o], xout = ZG, rule = 1)$y)
}, by = .(pid, P, step, var)]
A <- M[is.finite(d), .(d = mean(d, na.rm = TRUE), n = uniqueN(pid)), by = .(P, step, var, z)]
A[, var := factor(var, levels = c("Tair","Wind","RH","VPD"),
     labels = c("ΔAir temperature (°C)","ΔWind (m s⁻¹)","ΔRH (%)","ΔVPD (kPa)"))]
A[, step := factor(step, levels = names(STEPS), labels = c("ΔLAI","ΔHmax","ΔfCover","ΔLAD"))]
fwrite(A, file.path(OUT, "tab_trait_vertical_gradient_perplot.csv"))

PAL <- c("ΔLAI"="#2C7BB6","ΔHmax"="#E8A33D","ΔfCover"="#1A9850","ΔLAD"="#D46A9F")
# facet_grid shares the x scale down a column, which would squash VPD against RH.
# Build one row per variable and stack, so every variable keeps its own x range.
#' One variable's row of four archetype panels, on its own x scale
#' @param v variable level to draw
#' @param xlab x axis label
#' @param show_strip whether to draw the archetype strip labels (top row only)
#' @return a ggplot object
row <- function(v, xlab, show_strip) {
  ggplot(A[var == v], aes(d, z, colour = step)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
    geom_path(linewidth = 0.85) +
    facet_wrap(~ P, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = PAL, name = "Variable step") +
    coord_cartesian(ylim = c(0, 1.4)) +
    labs(x = xlab, y = NULL) + theme_article(11) +
    theme(legend.position = "none",
          strip.text = if (show_strip) element_text(face = "bold") else element_blank())
}
LV <- levels(A$var)
p <- (row(LV[1], LV[1], TRUE) / row(LV[2], LV[2], FALSE) /
      row(LV[3], LV[3], FALSE) / row(LV[4], LV[4], FALSE)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p <- patchwork::wrap_elements(p) +
  labs(tag = "Height  z / h_canopy  (1 = canopy top)") +
  theme(plot.tag = element_text(size = 11, angle = 90), plot.tag.position = "left")
ggsave_article(file.path(OUT, "fig_trait_vertical_gradient_perplot"), p, 11.5, 10)

cat("\n=== sub-canopy (z < 1) mean |effect| per archetype, largest trait first ===\n")
sub <- M[z < 1, .(m = mean(abs(d), na.rm = TRUE)), by = .(P, step, var)]
for (v in unique(sub$var)) {
  cat("\n", as.character(v), "\n")
  print(dcast(sub[var == v], step ~ P, value.var = "m"), digits = 2)
}
cat("\nDONE\n")
