# ==============================================================================
# DIRECT model-to-model comparison of MuSICA v3.2.0 vs v3.2.3 (yoyo / ABL_flag
# 'iter') -- WITHOUT going through the HOBO observations. For each of the 53
# plots both binaries are run on the same canopy + same forcing; here we pair
# their OUTPUT metrics plot-by-plot to isolate what the ABL coupling alone does
# to the simulated microclimate.
#   Run from repo root:
#   Rscript Chapitre1/comparaison_versions/scripts/02_analyses_microclimat/cmp_v320_v323_direct.R
# Out: figures/doc_media/image_directcmp.png  +  tables/tab_v320_v323_direct.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
DIR <- "Chapitre1/comparaison_versions"
TBL <- file.path(DIR, "tables"); MED <- file.path(DIR, "figures/doc_media")
PAL <- c(P1 = "#D9A441", P2 = "#7FBC41", P3 = "#3690C0", P4 = "#1B7837")

cl <- as.data.table(readRDS("outputs/figures_pipeline_z05/data/clusters.rds"))  # id_plot, Cluster (P1..P4)

# --- assemble paired (v320, v323) per-plot metrics from existing tables --------
rj <- fread(file.path(TBL, "recap_JJAS_plotdata.csv"))   # slope/equil/dTmaxS/dTmaxH _v320/_v323
ex <- fread(file.path(TBL, "compare_binaries_extra_axes.csv"))
ex[src == "v3.2.0", src := "v320"][src == "v3.2.3 iter", src := "v323"]
exw <- dcast(ex[src %in% c("v320", "v323")], id_plot ~ src,
             value.var = c("slope_day", "slope_night", "DTR", "dTmin", "bias"))
d <- merge(rj, exw, by = "id_plot")
d <- merge(d, cl, by = "id_plot", all.x = TRUE)
setnames(d, "Cluster", "P")

# metric registry: key, v320 col, v323 col, label, unit
M <- list(
  list("dTmaxH", "dTmaxH_v320", "dTmaxH_v323", "ΔTmax, hot days", "°C"),
  list("dTmaxS", "dTmaxS_v320", "dTmaxS_v323", "ΔTmax, summer",   "°C"),
  list("slope",  "slope_v320",  "slope_v323",  "slope (all hours)", ""),
  list("slopeN", "slope_night_v320", "slope_night_v323", "slope (night)", ""),
  list("dTmin",  "dTmin_v320",  "dTmin_v323",  "ΔTmin (night minima)", "°C"),
  list("DTR",    "DTR_v320",    "DTR_v323",    "DTR (diurnal range)", "°C")
)

# --- summary table over ALL axes ----------------------------------------------
allM <- c(M, list(
  list("slopeD", "slope_day_v320", "slope_day_v323", "slope (day)", ""),
  list("equil",  "equil_v320",  "equil_v323",  "equilibrium (unstable, slope~1)", "°C"),
  list("bias",   "bias_v320",   "bias_v323",   "warm bias vs macro", "°C")))
tab <- rbindlist(lapply(allM, function(m) {
  v0 <- d[[m[[2]]]]; v3 <- d[[m[[3]]]]; dd <- v3 - v0
  tt <- t.test(v3, v0, paired = TRUE)
  data.table(metric = m[[4]], unit = m[[5]],
             mean_v320 = mean(v0, na.rm = TRUE), mean_v323 = mean(v3, na.rm = TRUE),
             delta = mean(dd, na.rm = TRUE),
             ci_lo = unname(tt$conf.int[1]), ci_hi = unname(tt$conf.int[2]),
             r_models = cor(v0, v3, use = "complete.obs"), p_paired = tt$p.value)
}))
tab[, (names(tab)[3:9]) := lapply(.SD, round, 3), .SDcols = 3:9]
fwrite(tab, file.path(TBL, "tab_v320_v323_direct.csv"))
cat("=== Direct v3.2.0 vs v3.2.3 (model-to-model, paired over 53 plots) ===\n")
print(tab)

# --- 6-panel scatter v320 (x) vs v323 (y), 1:1 line ---------------------------
panel <- function(m) {
  v0 <- d[[m[[2]]]]; v3 <- d[[m[[3]]]]
  rng <- range(c(v0, v3), na.rm = TRUE)
  sub <- sprintf("Δ = %+.2f%s   r = %.2f", mean(v3 - v0, na.rm = TRUE),
                 ifelse(m[[5]] == "", "", paste0(" ", m[[5]])), cor(v0, v3, use = "complete.obs"))
  ggplot(d, aes(.data[[m[[2]]]], .data[[m[[3]]]], fill = P)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(shape = 21, colour = "black", size = 2, stroke = .3, alpha = .9) +
    scale_fill_manual(values = PAL, name = "Archetype") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = sprintf("v3.2.0  (%s)", m[[4]]), y = sprintf("v3.2.3  (%s)", m[[4]]),
         title = m[[4]], subtitle = sub) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"))
}
fig <- wrap_plots(lapply(M, panel), ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Direct comparison v3.2.0 vs v3.2.3 (yoyo) — model to model, 53 plots, no observations",
    subtitle = "Each point is one plot, same canopy and forcing run with both binaries. Dashed line = identity (v3.2.3 = v3.2.0). Points below the line: v3.2.3 lower. The yoyo damps the diurnal cycle: it cools daytime maxima and compresses DTR while raising night minima; models stay rank-correlated.",
    theme = theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 8.5, colour = "grey35"))) &
  theme(legend.position = "bottom")
ggsave(file.path(MED, "image_directcmp.png"), fig, width = 11, height = 7.6, dpi = 200, bg = "white")
cat(sprintf("\nDONE -> %s (n = %d plots)\n", file.path(MED, "image_directcmp.png"), nrow(d)))
