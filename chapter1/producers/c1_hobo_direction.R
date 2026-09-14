# ==============================================================================
# Does the MuSICA sensitivity go in the SAME DIRECTION as the HOBO field data?
# (Supervisor: "sensi MuSICA est-ce que les mesures HOBO vont dans le même sens ?")
#
# MuSICA sensitivity is ceteris-paribus; the 53-logger field relationship is
# CONFOUNDED (LAI, fCover, Hmax co-vary across plots). So the honest comparison
# is the SIGN of the response, not its magnitude (advisor, 2026-06-18):
#  * SIM sign  = sign of the per-cluster LOCAL ΔTmax sensitivity (tab_sensitivity).
#  * OBS sign  = sign of Spearman corr between each plot's observed mean ΔTmax
#                (53 HOBO, z05 branch) and its LiDAR variable (LAI, Hmax, fCover, VCI).
# Agreement = same sign. Reported pooled (53 plots) and per cluster where n allows.
#   Rscript c1_hobo_direction.R
# Out: out_files/Chapter1/tables/tab_hobo_direction.csv  + FigSh_hobo_direction.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(sf); library(terra)
  library(dplyr); library(tidyr); library(rmusica); library(musica.tools); source("R/cluster_relabel.R")
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))                              # CFG, build_hobo_inputs, load_lidar_rasters

# ---- per-HOBO observed ΔTmax (z05 branch) -----------------------------------
# ALIGNED CLOCK (2026-07-29). The observations formerly came from the z05 branch's
# ref_validation.rds, which predates the clock alignment: the loggers are on UTC and
# the forcing on solar time (UTC+1), so those values run ~0.3 degC warm (mean -0.730
# against -1.030) and this script OVERWROTE the canonical table with them on every
# run. They now come from the aligned validation (macro_ref_obs, see
# R/dtmax_convention.R and c1_hobo_native20_validation_regen.R), so the model-free
# analyses share one observational convention with the validation.
V   <- readRDS("outputs/figures_pipeline_z05/data/ref_validation.rds")   # still used for structure below
.AL <- as.data.table(data.table::fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv"))
obs <- .AL[, .(id_plot, dTmax_obs = obs_dt)]

# ---- per-HOBO LiDAR structural variables (rebuild as in pipeline 03z) --------
z05     <- fread(sprintf("in_files/lad_z05/Blois_lad_%s_r25.csv", "z05"))
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
hob_r   <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
hp      <- sf::st_read(CFG$hobo_geojson, quiet = TRUE)
hp      <- hp[!hp$id_plot %in% CFG$ids_to_remove, ]
hob_r$id_plot <- hp$id_plot
fcov    <- hob_r[, .(id_plot, fCover, VCI = if ("VCI" %in% names(hob_r)) VCI else NA_real_)]
str_h   <- merge(z05[, .(id_plot, LAI, Hmax)], fcov, by = "id_plot")   # LAI one-sided
str_h[fCover < 0.5, fCover := 0.5]

# attach the CANONICAL cluster of each HOBO plot (the same assignment Fig 5 uses).
# WARNING: do NOT re-derive by nearest cLHS-archetype centroid. The cLHS sample LAI is
# TWO-SIDED while the logger LAI (z05) is ONE-SIDED, so standardizing logger LAI against
# cLHS centroids makes every logger look artificially sparse and EMPTIES the dense
# archetype (P4=0) — a units bug. clusters.rds already carries P1..P4 labels.
clu <- as.data.table(readRDS("outputs/figures_pipeline_z05/data/clusters.rds"))[, .(id_plot, P = as.character(Cluster))]
VARS <- c("LAI","Hmax","fCover","VCI")

obs_sl <- .AL[, .(id_plot, slope_obs = obs_sl)]                       # aligned clock, as above
D <- merge(obs, str_h, by = "id_plot"); D <- merge(D, obs_sl, by = "id_plot", all.x = TRUE)
D <- merge(D, clu, by = "id_plot")
cat(sprintf("matched %d HOBO plots\n", nrow(D))); print(table(D$P))

# --- VCI from the wall-to-wall NORMALIZED 10 m raster (the old vci_res_10_m.tif was on
#     ABSOLUTE elevation, compressed ~0.68×; recomputed normalized → vci_res_10_m_norm_Blois.tif).
#     Same 10 m product as the cLHS sample, for consistency. ---
vrf <- "in_files/vci_res_10_m_norm_Blois.tif"
if (file.exists(vrf)) {
  co <- sf::st_coordinates(sf::st_geometry(hp))
  vcir <- data.table(id_plot = hp$id_plot, VCI_r = as.numeric(terra::extract(terra::rast(vrf), co[, 1:2])[, 1]))
  D <- merge(D, vcir, by = "id_plot", all.x = TRUE)
  D[is.finite(VCI_r), VCI := VCI_r][, VCI_r := NULL]
  cat(sprintf("VCI (raster 10 m normalisé) : %d/%d plots, médiane %.3f\n",
              sum(is.finite(D$VCI)), nrow(D), median(D$VCI, na.rm = TRUE)))
}
# DO NOT WRITE THE CANONICAL TABLE (2026-07-29). This script rebuilds the structural
# variables from the z05 pipeline branch, which differ materially from the canonical
# native20 ones (Hmax by up to 25 m, fCover mean 0.76 vs 0.94) AND re-derives P, so
# running it used to silently overwrite tab_hobo_perplot_cluster.csv, changing the
# archetype labels and cascading into every figure. The canonical table is now
# maintained by c1_align_cluster_obs.R, which touches ONLY the observed columns.
fwrite(D, "out_files/Chapter1/tables/tab_hobo_direction_z05branch.csv")

# ---- Does the MEASURED microclimate appear in the clustering? ----------------
# (supervisor: "est-ce que le microclimat mesuré se retrouve dans le clustering")
kw_t <- suppressWarnings(kruskal.test(dTmax_obs ~ P, D))
kw_s <- suppressWarnings(kruskal.test(slope_obs ~ P, D[is.finite(slope_obs)]))
cat(sprintf("\nObserved microclimate by cluster — Kruskal-Wallis:\n  dTmax_obs: chi2=%.1f p=%.2g | slope_obs: chi2=%.1f p=%.2g\n",
            kw_t$statistic, kw_t$p.value, kw_s$statistic, kw_s$p.value))
print(D[, .(n=.N, dTmax_obs=round(mean(dTmax_obs,na.rm=TRUE),2), slope_obs=round(mean(slope_obs,na.rm=TRUE),3)), by=P][order(P)])
MLAB <- c(dTmax_obs="Observed~Delta*T[max]~('°C')", slope_obs="Observed~micro/macro~slope")
mb <- melt(D, id.vars="P", measure.vars=names(MLAB), variable.name="metric", value.name="val")[is.finite(val)]
mb[, metric := factor(metric, levels=names(MLAB), labels=MLAB)]
set.seed(1)
pmc <- ggplot(mb, aes(P, val, fill=P, colour=P)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.4) +
  geom_boxplot(width=0.6, alpha=0.35, outlier.shape=NA, linewidth=0.6) +
  geom_jitter(width=0.12, size=1.1, alpha=0.5) +
  facet_wrap(~metric, scales="free_y", labeller=label_parsed) +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") + scale_colour_manual(values=PAL_CLUSTER, guide="none") +
  labs(x=NULL, y=NULL, title="Measured microclimate (53 HOBO) by structural cluster") +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11, face="bold"))
ggsave("out_files/Chapter1/figures/FigSh_microclimate_by_cluster.png", pmc, width=8.5, height=4.4, dpi=200, bg="white")
cat("saved FigSh_microclimate_by_cluster.png\n")

# ---- observed sign: Spearman(dTmax_obs, variable) ---------------------------
obs_sign <- rbindlist(lapply(VARS, function(v) {
  ok <- is.finite(D$dTmax_obs) & is.finite(D[[v]])
  ct <- suppressWarnings(cor.test(D$dTmax_obs[ok], D[[v]][ok], method = "spearman"))
  data.table(variable = v, rho_obs = unname(ct$estimate), p_obs = ct$p.value, n = sum(ok)) }))

# ---- simulated sign: pooled local ΔTmax sensitivity (mean over clusters) -----
sens <- fread("out_files/Chapter1/tables/tab_sensitivity_percluster.csv")[scope == "local"]
sens[variable == "LAD", variable := "VCI"]                  # shape term: sim "LAD" <-> obs "VCI"
sim_sign <- sens[, .(sens_sim = mean(Tmax_all, na.rm = TRUE)), by = variable]

cmp <- merge(obs_sign, sim_sign, by = "variable")
# VCI field correlation is confounded by collinearity with density (LAI), unlike
# the ceteris-paribus sim term -> flag it as not cleanly sign-comparable.
cmp[, note := ifelse(variable == "VCI", "collinear w/ LAI in field", "ceteris-paribus comparable")]
cmp[, `:=`(sign_obs = sign(rho_obs), sign_sim = sign(sens_sim))]
cmp[, agree := ifelse(sign_obs == sign_sim, "YES", "NO")]
cmp[, variable := factor(variable, levels = VARS)]
setorder(cmp, variable)
fwrite(cmp, "out_files/Chapter1/tables/tab_hobo_direction.csv")
cat("\n=== SIGN agreement MuSICA sensitivity vs HOBO field correlation ===\n")
print(cmp[, .(variable, rho_obs = round(rho_obs,2), p_obs = round(p_obs,3),
              sens_sim = round(sens_sim,2), agree, note)])

# ---- figure: observed scatter dTmax_obs vs each variable, sign annotated -----
DL <- melt(D, id.vars = c("id_plot","dTmax_obs","P"), measure.vars = VARS,
           variable.name = "variable", value.name = "x")
DL[, variable := factor(variable, levels = VARS)]
labs_agree <- cmp[, .(variable, lab = sprintf("rho[obs]==%.2f~~(%s)", rho_obs, agree))]
p <- ggplot(DL, aes(x, dTmax_obs)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.7) +
  geom_point(aes(colour = P), size = 1.8, alpha = 0.8) +
  geom_text(data = labs_agree, aes(x = -Inf, y = Inf, label = lab), parse = TRUE,
            hjust = -0.08, vjust = 1.3, size = 3.4, inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = PAL_CLUSTER, name = "Cluster") +
  labs(x = "LiDAR variable (per HOBO plot)", y = "Observed mean ΔTmax (°C, 53 HOBO)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave("out_files/Chapter1/figures/FigSh_hobo_direction.png", p, width = 11, height = 4.4, dpi = 200, bg = "white")
cat("\nDONE -> tab_hobo_direction.csv + FigSh_hobo_direction.png\n")
