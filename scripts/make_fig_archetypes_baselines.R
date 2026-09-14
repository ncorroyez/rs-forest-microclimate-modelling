# ==============================================================================
# fig_archetypes_baselines — per archetype, the mean REAL profile (solid) and the
# TWO uniform LAD baselines used in the attribution:
#   - per-cluster baseline (MAIN): uniform at that archetype's mean LAI & Hmax (dashed)
#   - global baseline (COMPLEMENTARY): uniform at the landscape mean LAI & Hmax,
#     identical in every panel (dotted).
# Makes explicit which baseline each scheme uses per profile (Section 2.4 / 2.6).
#   Rscript scripts/make_fig_archetypes_baselines.R
# Out: outputs/figs_MEB2026_final/fig_archetypes_baselines.{png,pdf}
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(cli); library(tidyverse) })
OUT <- here::here("outputs/figs_MEB2026_final"); source(here::here("R/cluster_relabel.R"))

df <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
lad_cols <- grep("^LAD_Layer_", names(df), value = TRUE)
z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

# --- mean REAL profile per cluster (display PAD = stored/2) --------------------
real <- rbindlist(lapply(c("1","2","3","4"), function(cl) {
  m <- as.matrix(df[as.character(Cluster) == cl, ..lad_cols]) / 2; m[is.na(m)] <- 0
  data.table(Cluster_old = cl, height = z_breaks, lad = colMeans(m), lad_sd = apply(m, 2, sd))
}))
real[, Cluster := relabel_cluster(Cluster_old)][, type := "Real (mean)"]

box <- function(dens, hm) data.table(height = c(0, hm, hm), lad = c(dens, dens, 0))

# --- per-cluster baseline (MAIN) : uniform at each archetype's mean LAI & Hmax -
hmax <- df[, .(Hmax = mean(Hmax, na.rm = TRUE)), by = Cluster][, Cluster := relabel_cluster(Cluster)]
perclu <- rbindlist(lapply(unique(real$Cluster), function(cl) {
  tot <- sum(real[Cluster == cl, lad]); hm <- hmax[Cluster == cl, Hmax]
  cbind(Cluster = cl, box(tot / hm, hm), type = "Per-cluster baseline (main)")
}))

# --- global baseline (COMPLEMENTARY) : uniform at landscape mean, same everywhere
gm   <- colMeans(as.matrix(df[, ..lad_cols]) / 2, na.rm = TRUE); gm[is.na(gm)] <- 0
gtot <- sum(gm); ghmax <- mean(df$Hmax, na.rm = TRUE); gdens <- gtot / ghmax
glob <- rbindlist(lapply(unique(real$Cluster), function(cl)
  cbind(Cluster = cl, box(gdens, ghmax), type = "Global baseline (complementary)")))
cli_alert("Global baseline: density {round(gdens,3)} m2 m-3 up to Hmax {round(ghmax,1)} m (same in every panel)")

allp <- rbind(real[, .(Cluster, height, lad, type)], perclu, glob)
lev  <- c("Real (mean)", "Per-cluster baseline (main)", "Global baseline (complementary)")
allp[, type := factor(type, levels = lev)]
max_lad <- max(real$lad + real$lad_sd, na.rm = TRUE) * 1.05

p <- ggplot(allp, aes(x = lad, y = height, colour = Cluster)) +
  geom_ribbon(data = real, aes(y = height, xmin = pmax(0, lad - lad_sd), xmax = lad + lad_sd, fill = Cluster),
              alpha = 0.20, colour = NA, inherit.aes = FALSE, orientation = "y") +
  geom_path(aes(linetype = type), linewidth = 1.25) +
  facet_wrap(~ Cluster, nrow = 1L) +
  scale_colour_manual(values = PAL_CLUSTER, guide = "none", drop = FALSE) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none", drop = FALSE) +
  scale_linetype_manual(values = c("Real (mean)" = "solid",
                                   "Per-cluster baseline (main)" = "dashed",
                                   "Global baseline (complementary)" = "dotted"), name = NULL) +
  scale_x_continuous(limits = c(0, max_lad), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 42), breaks = c(0, 10, 20, 30, 40)) +
  labs(x = bquote("LAD" ~ "(m"^"2" ~ "m"^"-3" * ")"), y = "Height (m)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 22),
        axis.title = element_text(face = "bold", size = 22), axis.text = element_text(size = 16),
        legend.position = "bottom", legend.text = element_text(size = 17),
        panel.grid.minor = element_blank()) +
  guides(linetype = guide_legend(nrow = 1, override.aes = list(colour = "grey25")))

ggsave(file.path(OUT, "fig_archetypes_baselines.png"), p, width = 18, height = 7.6, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_archetypes_baselines.pdf"), p, width = 18, height = 7.6, device = cairo_pdf)
cli_alert_success("Saved fig_archetypes_baselines (PNG 600dpi + PDF)")
