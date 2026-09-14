# ==============================================================================
# (1) CR robustness to the lower/upper split threshold (0.10/0.25/0.50*Hmax).
# (2) Hmax contrast: deciduous Blois (this study) vs Bouwen's sparse pine
#     (she reported hmax a POOR predictor; we test it here).
# Reuses per-plot slopes from tab_canopy_ratio_perplot.csv; only CR is recomputed.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(sf); library(terra); library(tidyverse)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/validation.R")); source(here::here("R/cluster_relabel.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
PP  <- fread(file.path(OUT, "tab_canopy_ratio_perplot.csv"))   # id_plot, LAI, VCI, fCover, Hmax, Cluster, log_slope_obs, log_slope_sim
if ("CR" %in% names(PP)) PP[, CR := NULL]   # drop old CR; recomputed per threshold

# ---- Rebuild LAD profiles (for CR only) ------------------------------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                            CFG$ids_to_remove, hobo_buffer_mode = "buffer25"))
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
df_hobo$id_plot <- hobo_pts$id_plot
lad_cols <- grep("^LAD_Layer_", names(df_hobo), value = TRUE)
z_lad    <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

cr_at <- function(split) {
  sapply(seq_len(nrow(df_hobo)), function(i) {
    d <- as.numeric(df_hobo[i, ..lad_cols]); d[is.na(d)] <- 0
    hmax <- df_hobo$Hmax[i]
    if (hmax <= 0 || sum(d) <= 0) return(NA_real_)
    hs <- split * hmax
    lo <- sum(d[z_lad <= hs]); up <- sum(d[z_lad > hs])
    if (up <= 0) return(NA_real_); lo / up
  })
}

# ---- (1) Robustness across thresholds --------------------------------------
cli_h1("(1) CR robustness to split threshold")
rob <- rbindlist(lapply(c(0.10, 0.25, 0.50), function(s) {
  cr <- data.table(id_plot = df_hobo$id_plot, CR = cr_at(s))
  d  <- merge(PP, cr, by = "id_plot")
  d  <- d[is.finite(CR) & is.finite(log_slope_obs)]
  ct  <- suppressWarnings(cor.test(d$CR, d$log_slope_obs, method = "spearman", exact = FALSE))
  m0  <- lm(log_slope_obs ~ LAI, d); m1 <- lm(log_slope_obs ~ LAI + CR, d)
  an  <- anova(m0, m1)
  data.table(split = sprintf("%.2f*Hmax", s),
             CR_med = median(d$CR), CR_share_gt1 = mean(d$CR > 1),
             rho_obs = unname(ct$estimate), p_obs = ct$p.value,
             dR2_beyond_LAI = summary(m1)$r.squared - summary(m0)$r.squared,
             p_CR_beyond_LAI = an$`Pr(>F)`[2])
}))
print(rob)
fwrite(rob, file.path(OUT, "tab_cr_robustness_threshold.csv"))

# ---- (2) Hmax contrast vs Bouwen pine --------------------------------------
cli_h1("(2) Hmax as a predictor — deciduous Blois vs Bouwen sparse pine")
d <- PP[is.finite(Hmax) & is.finite(log_slope_obs)]
for (v in c("Hmax", "LAI")) {
  ct <- suppressWarnings(cor.test(d[[v]], d$log_slope_obs, method = "spearman", exact = FALSE))
  cli_alert("{v}: rho={round(unname(ct$estimate),2)}  p={signif(ct$p.value,2)}")
}
# Is Hmax still predictive once LAI is known? (her claim: hmax redundant/poor)
m_lai      <- lm(log_slope_obs ~ LAI, d)
m_lai_hmax <- lm(log_slope_obs ~ LAI + Hmax, d)
an <- anova(m_lai, m_lai_hmax)
cli_alert("Hmax beyond LAI: dR2={round(summary(m_lai_hmax)$r.squared - summary(m_lai)$r.squared,4)}  p={signif(an$`Pr(>F)`[2],3)}")
# Spearman Hmax~LAI (in pine she said the usual taller=more-leaf link 'did not hold')
ct_hl <- suppressWarnings(cor.test(d$Hmax, d$LAI, method = "spearman", exact = FALSE))
cli_alert("Hmax~LAI coupling (Blois oak): rho={round(unname(ct_hl$estimate),2)}  p={signif(ct_hl$p.value,2)}")

hmax_tab <- data.table(
  quantity = c("rho(Hmax, logslope_obs)", "rho(LAI, logslope_obs)",
               "Hmax added beyond LAI (dR2)", "p(Hmax beyond LAI)",
               "rho(Hmax, LAI) coupling"),
  value = c(round(cor(d$Hmax, d$log_slope_obs, method="spearman"),3),
            round(cor(d$LAI,  d$log_slope_obs, method="spearman"),3),
            round(summary(m_lai_hmax)$r.squared - summary(m_lai)$r.squared,4),
            signif(an$`Pr(>F)`[2],3),
            round(unname(ct_hl$estimate),3)))
print(hmax_tab)
fwrite(hmax_tab, file.path(OUT, "tab_hmax_contrast_vs_pine.csv"))
