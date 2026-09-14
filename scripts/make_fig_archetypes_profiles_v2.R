# ==============================================================================
# fig_archetypes_profiles — V2 cLHS (LAI ×2 corrigé), labels P1..P4, palette
# harmonisée avec les autres figures MEB :
#   P1 = #1B9E77 (vert)
#   P2 = #D95F02 (orange)
#   P3 = #7570B3 (violet)
#   P4 = #E7298A (rose)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(tidyverse)
})

OUT <- here::here("outputs/figs_MEB2026_final")

source(here::here("R/cluster_relabel.R"))

cli_h1("fig_archetypes_profiles (V2, clusters ordonnes par LAI croissant)")
df_sample <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))

# Floor fCover at 0.5 already done in V2 ; keep as-is
lad_cols <- grep("^LAD_Layer_", names(df_sample), value = TRUE)
z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

# Build LAD profile mean ± sd per cluster (old IDs).
# Same convention as LAI : V2 stores LAD_Layer_* × 2 (PAD→LAI). For DISPLAY we
# divide by 2 to recover the raw lidR PAD profile, consistent with the LAI label.
profiles <- list()
for (cl in c("1","2","3","4")) {
  sub <- df_sample[as.character(Cluster) == cl, ..lad_cols]
  sub_mat <- as.matrix(sub) / 2     # raw lidR PAD
  sub_mat[is.na(sub_mat)] <- 0
  profiles[[cl]] <- data.table(
    Cluster_old = cl,
    height = z_breaks,
    lad_mean = colMeans(sub_mat),
    lad_sd   = apply(sub_mat, 2, sd))
}
prof_long <- rbindlist(profiles)
prof_long[, Cluster := relabel_cluster(Cluster_old)]

# Trait stats per cluster.
# V2 sample stores LAI ×2 (PAD→LAI spherical leaf angle correction) for use
# inside MuSICA. For DISPLAY we want the raw lidR PAD value : LAI / 2.
# Hmax and fCover are unaffected.
stats <- df_sample[, .(LAI_mean = mean(LAI, na.rm=TRUE) / 2,
                          LAI_sd   = sd(LAI, na.rm=TRUE) / 2,
                          Hmax_mean = mean(Hmax, na.rm=TRUE),
                          Hmax_sd   = sd(Hmax, na.rm=TRUE),
                          fCover_mean = mean(fCover, na.rm=TRUE),
                          fCover_sd   = sd(fCover, na.rm=TRUE)),
                      by = Cluster]
stats[, Cluster := relabel_cluster(Cluster)]
setorder(stats, Cluster)
cli_alert("Cluster stats (raw lidR LAI = V2/2, for display) :"); print(stats)

# Plotmath labels (3 lines per cluster)
stats[, lab_LAI    := sprintf("LAI == %.2f %%+-%% %.2f", LAI_mean, LAI_sd)]
stats[, lab_Hmax   := sprintf("H[max] == %.1f %%+-%% %.1f * ' m'",
                                Hmax_mean, Hmax_sd)]
stats[, lab_fCover := sprintf("fCover == %.2f %%+-%% %.2f",
                                fCover_mean, fCover_sd)]

max_lad <- max(prof_long$lad_mean + prof_long$lad_sd, na.rm = TRUE) * 1.05
x_ann   <- max_lad * 0.30   # annotation x position
ANNOT_SIZE <- 6.5

p <- ggplot(prof_long, aes(x = lad_mean, y = height,
                            colour = Cluster, fill = Cluster)) +
  geom_ribbon(aes(xmin = pmax(0, lad_mean - lad_sd),
                    xmax = lad_mean + lad_sd),
                alpha = 0.25, colour = NA) +
  geom_path(linewidth = 1.4) +
  geom_text(data = stats, aes(x = x_ann, y = 38, label = lab_LAI,
                                   colour = Cluster),
             parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 0,
             size = ANNOT_SIZE, fontface = "bold", show.legend = FALSE) +
  geom_text(data = stats, aes(x = x_ann, y = 33, label = lab_Hmax,
                                   colour = Cluster),
             parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 0,
             size = ANNOT_SIZE, fontface = "bold", show.legend = FALSE) +
  geom_text(data = stats, aes(x = x_ann, y = 28, label = lab_fCover,
                                   colour = Cluster),
             parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 0,
             size = ANNOT_SIZE, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Cluster, nrow = 1L) +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL, drop = FALSE) +
  scale_fill_manual(values = PAL_CLUSTER, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(0, max_lad),
                      expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 42),
                      breaks = c(0, 10, 20, 30, 40)) +
  labs(x = bquote("LAD" ~ "(m"^"2"~"m"^"-3"*")"),
       y = "Height (m)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text  = element_text(face = "bold", size = 22),
        axis.title  = element_text(face = "bold", size = 22),
        axis.text   = element_text(size = 16),
        legend.position = "none",
        panel.grid.minor = element_blank())

# Slide-quality outputs (PNG 600 dpi + PDF vector)
ggsave(file.path(OUT, "fig_archetypes_profiles.png"),
        p, width = 18, height = 7, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_archetypes_profiles.pdf"),
        p, width = 18, height = 7, device = cairo_pdf)
cli_alert_success("Saved fig_archetypes_profiles (PNG 600dpi + PDF)")

# Also dump the palette in a small CSV so I have a single source of truth
fwrite(data.table(Cluster = names(PAL_CLUSTER), Hex = PAL_CLUSTER),
        file.path(OUT, "cluster_palette.csv"))
cli_alert_success("Saved cluster_palette.csv (single source of truth)")
