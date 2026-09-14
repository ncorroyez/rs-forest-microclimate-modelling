# ==============================================================================
# PIPELINE STAGE 08 — Typology figures: FPCA (LAD-shape modes) + clusters.
# Reuses the cached FPCA (stage 01) and the FROZEN cLHS sample (clusters).
# Each figure is wrapped in tryCatch so one failure doesn't abort the stage.
# Outputs (OUT_FIG): fig_fpca_{harmonics,loadings,gcv}.{png,pdf},
#   fig_clusters_lad_profiles.*, fig_kmeans_elbow.*, fig_archetype_profiles.*
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
suppressMessages({ library(terra); library(fda); library(patchwork) })
cli_h1("STAGE 08 — typology (FPCA + clusters)")

save_fig <- function(p, name, w=10, h=7) {
  if (is.null(p)) { cli_alert_warning("{name}: NULL — skipped"); return(invisible()) }
  ggsave(file.path(PIPE$OUT_FIG, paste0(name,".png")), p, width=w, height=h, dpi=300, bg="white")
  ggsave(file.path(PIPE$OUT_FIG, paste0(name,".pdf")), p, width=w, height=h, device=cairo_pdf)
  cli_alert_success("Saved {name}")
}
try_fig <- function(expr, name, ...) {
  p <- tryCatch(expr, error=function(e){ cli_alert_warning("{name} FAILED: {e$message}"); NULL })
  save_fig(p, name, ...)
}

# ---- FPCA (from stage 01 cache) ---------------------------------------------
fc <- readRDS(file.path(PIPE$OUT_DATA, "forest_fpca.rds"))
if (is.null(fc$fpca)) cli_alert_warning("No FPCA in cache — rerun 01_data.") else {
  try_fig(plot_fpc_harmonics(fc$fpca, NULL, fc$z_break, fc$forest$Hmax), "fig_fpca_harmonics", 13, 6)
  try_fig(plot_fpc_loadings(fc$fpca), "fig_fpca_loadings", 11, 6)
  try_fig(fc$fpca$gcv_plot, "fig_fpca_gcv", 8, 6)
  cli_alert("FPCA varprop (FPC1-3): {paste(round(100*fc$fpca$varprop[1:3],1), collapse='/')} %")
}

# ---- Clusters (FROZEN sample, relabelled P1..P4) ----------------------------
PAL_CLU <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")
samp <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))
samp[, Cluster := relabel_cluster(Cluster)]
samp[, Cluster := factor(Cluster, levels=paste0("P",1:4))]

# LAD profiles per cluster — cluster palette, NO legend, NO in-plot stats text
try_fig({
  lad_cols <- grep("^LAD_Layer_", names(samp), value=TRUE)
  samp[, pid := .I]
  long <- melt(samp, id.vars=c("pid","Cluster"), measure.vars=lad_cols,
               variable.name="h", value.name="density")
  long[, height := as.numeric(gsub("LAD_Layer_","",h))]
  meanp <- long[, .(density=mean(density,na.rm=TRUE)), by=.(Cluster,height)]
  ggplot() +
    geom_path(data=long, aes(density, height, group=pid, colour=Cluster), alpha=0.06) +
    geom_path(data=meanp, aes(density, height, colour=Cluster), linewidth=1.3) +
    facet_wrap(~Cluster, nrow=1) +
    scale_colour_manual(values=PAL_CLU) +
    labs(x=expression(LAD~(m^2~m^{-3})), y="Height above ground (m)") +
    theme_bw(base_size=17) +
    theme(legend.position="none", panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold", size=18))
}, "fig_clusters_lad_profiles", 13, 6)

# EXTERNAL table: per-cluster mean [min–max] of structural metrics
mr <- function(v) sprintf("%.2f [%.2f–%.2f]", mean(v), min(v), max(v))
tab_struct <- samp[, .(n=.N,
                       LAI    = mr(LAI),
                       fCover = mr(fCover),
                       Hmax   = mr(Hmax),
                       VCI    = if ("VCI" %in% names(samp)) mr(VCI) else NA_character_),
                   by=Cluster][order(Cluster)]
fwrite(tab_struct, file.path(PIPE$OUT_TAB, "tab_cluster_structure.csv"))
cli_alert("Cluster structure table:"); print(tab_struct)

# k-means elbow (descriptive) on structural vars
try_fig({
  sc <- scale(as.data.frame(samp)[, c("LAI","Hmax","fCover")])
  el <- optimise_k_elbow(sc, k_range=1:10); el$plot
}, "fig_kmeans_elbow", 7, 5)

# archetype mean-shape profiles (one synthetic plot per cluster)
try_fig({
  arch <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "archetypes.rds")))
  arch[, Cluster := relabel_cluster(Cluster)]
  plot_cluster_all_profiles(as.data.frame(arch))
}, "fig_archetype_profiles", 12, 7)

# ---- cluster map (reuse existing self-contained script) ---------------------
tryCatch({ source(here::here("scripts/make_map_clusters_blois.R"))
           cli_alert_success("Cluster map via make_map_clusters_blois.R") },
         error=function(e) cli_alert_warning("cluster map skipped: {e$message}"))

cli_alert_success("Stage 08 done.")
