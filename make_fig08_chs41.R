#!/usr/bin/env Rscript
# ==============================================================================
# fig08_ch3_stratified.png — CSI3 rapport asymmetric verdict figure
# Redraw on the CHS41 hybrid forcing numbers (2026-08-24 rerun).
# Scenarios : STATIC_ALS (LiDAR direct) vs STATIC_S2_ATBD (Sentinel-2 direct).
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

TAB <- "/home/corroyez/Documents/z_Example_rmusica_31012025/out_files/Chapter3_CHS41/tables/figI_all_stratified.csv"
OUT <- "/home/corroyez/Documents/Admin/ED_GAIA/CSI3/figures/fig08_ch3_stratified.png"

d <- fread(TAB)
pick <- d[scenario %in% c("STATIC_ALS", "STATIC_S2_ATBD")]

df <- rbindlist(list(
  data.table(sensor="LiDAR",      stratum="Open canopy",   n=28, R2=pick[scenario=="STATIC_ALS",       dT_open]),
  data.table(sensor="LiDAR",      stratum="Dense canopy",  n=25, R2=pick[scenario=="STATIC_ALS",       dT_dense]),
  data.table(sensor="LiDAR",      stratum="Pooled (all 53)", n=53, R2=pick[scenario=="STATIC_ALS",     dT_pool]),
  data.table(sensor="Sentinel-2", stratum="Open canopy",   n=28, R2=pick[scenario=="STATIC_S2_ATBD",   dT_open]),
  data.table(sensor="Sentinel-2", stratum="Dense canopy",  n=25, R2=pick[scenario=="STATIC_S2_ATBD",   dT_dense]),
  data.table(sensor="Sentinel-2", stratum="Pooled (all 53)", n=53, R2=pick[scenario=="STATIC_S2_ATBD", dT_pool])
))
df[, stratum := factor(stratum, levels = c("Open canopy", "Dense canopy", "Pooled (all 53)"))]
df[, sensor  := factor(sensor,  levels = c("LiDAR", "Sentinel-2"))]
df[, lab := sprintf("%.2f", R2)]

cols <- c(LiDAR = "#1F6FA8", `Sentinel-2` = "#D2691E")

g <- ggplot(df, aes(x = stratum, y = R2, fill = sensor)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.62, colour = NA) +
  geom_text(aes(label = lab),
            position = position_dodge(width = 0.75),
            vjust = -0.35, size = 4.2, colour = "grey15") +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "LiDAR ranks dense; Sentinel-2's open-canopy edge is narrow",
    subtitle = "Amplitude ΔTmax (R²) — split at LAI 3.86, CHS41 station-based hybrid forcing",
    y        = "Between-plot R² (sim vs observed)",
    x        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle    = element_text(size = 11, colour = "grey35", hjust = 0,
                                    margin = margin(b = 12)),
    legend.position  = "top",
    legend.text      = element_text(size = 12),
    legend.key.width = unit(1.2, "lines"),
    axis.title.y     = element_text(size = 11, margin = margin(r = 8)),
    axis.text.x      = element_text(size = 11, colour = "grey20"),
    axis.text.y      = element_text(size = 10, colour = "grey40"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92"),
    plot.margin      = margin(14, 18, 12, 14)
  )

ggsave(OUT, g, width = 9.5, height = 4.4, dpi = 300, bg = "white")
cat("WROTE ", OUT, "\n")
