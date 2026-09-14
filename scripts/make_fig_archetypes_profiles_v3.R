# ==============================================================================
# fig_archetypes_profiles V3 — Hmax-coherent mean profiles.
#
# FIX (supervisor 2026-06-18): V2 averaged LAD at FIXED absolute height across
# plots with DIFFERENT Hmax, so tall plots smeared density up to ~40 m and the
# mean profile extended far above the cluster's mean Hmax ("Hmax bizarre").
#
# V3: rescale EACH plot's profile to the cluster-mean Hmax via the .lad_rescale
# homothety (stretch/compress the height axis, renormalize density so each plot's
# LAI is preserved), THEN average on the common grid. The mean profile now tops
# out exactly at the cluster mean Hmax, the mean SHAPE is undistorted, and
# "varying Hmax stretches the profile" is built in (height is relative).
# Display convention unchanged: raw lidR LAI/LAD = stored /2 (one-sided).
#   Rscript scripts/make_fig_archetypes_profiles_v3.R
# Out: outputs/figs_MEB2026_final/fig_archetypes_profiles_v3.{png,pdf}
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(cli); library(tidyverse) })
OUT <- here::here("outputs/figs_MEB2026_final")
source(here::here("R/cluster_relabel.R"))
source(here::here("R/lad.R"))                       # .lad_rescale (homothety)

cli_h1("fig_archetypes_profiles V3 (Hmax-coherent, relative-height averaging)")
df <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
df <- df[is.finite(Hmax) & Hmax > 0]
lad_cols <- grep("^LAD_Layer_", names(df), value = TRUE)
z_real   <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

# --- per cluster: rescale every plot to the cluster-mean Hmax, then average ----
profiles <- list(); stats <- list()
for (cl in c("1","2","3","4")) {
  sub <- df[as.character(Cluster) == cl]
  hmax_bar <- mean(sub$Hmax, na.rm = TRUE)
  grid <- seq_len(ceiling(hmax_bar))                # common absolute grid, top = mean Hmax
  M <- matrix(0, nrow = nrow(sub), ncol = length(grid))
  for (i in seq_len(nrow(sub))) {
    d  <- as.numeric(sub[i, ..lad_cols]) / 2; d[is.na(d)] <- 0   # raw lidR PAD (display)
    lai_i <- sum(d)
    r  <- .lad_rescale(z_real, d, as.numeric(sub$Hmax[i]), hmax_bar, lai_i)  # stretch to mean Hmax
    M[i, ] <- approx(r$height, r$density, grid, rule = 2)$y
  }
  profiles[[cl]] <- data.table(Cluster_old = cl, height = grid,
                               lad_mean = colMeans(M, na.rm = TRUE),
                               lad_sd   = apply(M, 2, sd, na.rm = TRUE))
  stats[[cl]] <- data.table(Cluster_old = cl,
    LAI_mean = mean(sub$LAI)/2, LAI_sd = sd(sub$LAI)/2,
    Hmax_mean = hmax_bar, Hmax_sd = sd(sub$Hmax),
    fCover_mean = mean(sub$fCover), fCover_sd = sd(sub$fCover))
}
prof_long <- rbindlist(profiles); prof_long[, Cluster := relabel_cluster(Cluster_old)]
stats <- rbindlist(stats); stats[, Cluster := relabel_cluster(Cluster_old)]; setorder(stats, Cluster)
cli_alert("Cluster stats (display LAI = stored/2). Profile top now = mean Hmax:"); print(stats[, .(Cluster, LAI_mean=round(LAI_mean,2), Hmax_mean=round(Hmax_mean,1), fCover_mean=round(fCover_mean,2))])

stats[, lab_LAI    := sprintf("LAI == %.2f %%+-%% %.2f", LAI_mean, LAI_sd)]
stats[, lab_Hmax   := sprintf("H[max] == %.1f %%+-%% %.1f * ' m'", Hmax_mean, Hmax_sd)]
stats[, lab_fCover := sprintf("fCover == %.2f %%+-%% %.2f", fCover_mean, fCover_sd)]

max_lad <- max(prof_long$lad_mean + prof_long$lad_sd, na.rm = TRUE) * 1.05
ymax    <- max(stats$Hmax_mean) * 1.12
x_ann   <- max_lad * 0.32
ANN <- 6.5
p <- ggplot(prof_long, aes(x = lad_mean, y = height, colour = Cluster, fill = Cluster)) +
  geom_ribbon(aes(xmin = pmax(0, lad_mean - lad_sd), xmax = lad_mean + lad_sd), alpha = 0.25, colour = NA) +
  geom_path(linewidth = 1.4) +
  # dashed marker at each cluster's mean Hmax (top of the profile)
  geom_segment(data = stats, aes(x = 0, xend = max_lad*0.9, y = Hmax_mean, yend = Hmax_mean, colour = Cluster),
               linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = stats, aes(x = x_ann, y = ymax*0.96, label = lab_LAI, colour = Cluster),
            parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1, size = ANN, fontface = "bold", show.legend = FALSE) +
  geom_text(data = stats, aes(x = x_ann, y = ymax*0.86, label = lab_Hmax, colour = Cluster),
            parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1, size = ANN, fontface = "bold", show.legend = FALSE) +
  geom_text(data = stats, aes(x = x_ann, y = ymax*0.76, label = lab_fCover, colour = Cluster),
            parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1, size = ANN, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Cluster, nrow = 1L) +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL, drop = FALSE) +
  scale_fill_manual(values = PAL_CLUSTER, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(0, max_lad), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, ymax)) +
  labs(x = bquote("LAD" ~ "(m"^"2"~"m"^"-3"*")"), y = "Height (m)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 22),
        axis.title = element_text(face = "bold", size = 22), axis.text = element_text(size = 16),
        legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(OUT, "fig_archetypes_profiles_v3.png"), p, width = 18, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_archetypes_profiles_v3.pdf"), p, width = 18, height = 7, device = cairo_pdf)
cli_alert_success("Saved fig_archetypes_profiles_v3 (PNG + PDF)")
