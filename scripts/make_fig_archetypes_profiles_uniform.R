# ==============================================================================
# fig_archetypes_profiles_uniform — like fig_archetypes_profiles (V2) but adds,
# per archetype, the UNIFORM LAD baseline used in the Shapley attribution:
# a vertically constant profile with the SAME total leaf area (LAI) and SAME
# Hmax as the mean real profile, only the shape flattened (dashed). Shows what
# the "LAD = uniform" baseline (Section 2.4) means visually.
#   Rscript scripts/make_fig_archetypes_profiles_uniform.R
# Out: outputs/figs_MEB2026_final/fig_archetypes_profiles_uniform.{png,pdf}
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(cli); library(tidyverse) })
OUT <- here::here("outputs/figs_MEB2026_final"); source(here::here("R/cluster_relabel.R"))

df <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
lad_cols <- grep("^LAD_Layer_", names(df), value = TRUE)
z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))   # 1.5, 2.5, ... (1 m bins)

# --- mean REAL profile per cluster (display = raw lidR PAD = stored/2) ---------
real <- rbindlist(lapply(c("1","2","3","4"), function(cl) {
  m <- as.matrix(df[as.character(Cluster) == cl, ..lad_cols]) / 2; m[is.na(m)] <- 0
  data.table(Cluster_old = cl, height = z_breaks, lad = colMeans(m), lad_sd = apply(m, 2, sd))
}))
real[, Cluster := relabel_cluster(Cluster_old)][, type := "Real (mean)"]

# --- UNIFORM baseline per cluster : same integrated PAD, same Hmax ------------
# density = (total integrated PAD over the mean profile) / Hmax_mean, constant
# from ground to Hmax_mean (then 0). Drawn as the right+top edge of the box.
hmax  <- df[, .(Hmax = mean(Hmax, na.rm = TRUE)), by = Cluster][, Cluster := relabel_cluster(Cluster)]
unif  <- rbindlist(lapply(unique(real$Cluster), function(cl) {
  tot <- sum(real[Cluster == cl, lad])                 # integral (1 m bins -> sum)
  hm  <- hmax[Cluster == cl, Hmax]; dens <- tot / hm
  data.table(Cluster = cl,
             height = c(0, hm, hm), lad = c(dens, dens, 0), lad_sd = NA_real_,
             type = "Uniform baseline")
}))
allp <- rbind(real[, .(Cluster, height, lad, lad_sd, type)], unif)
allp[, type := factor(type, levels = c("Real (mean)", "Uniform baseline"))]

max_lad <- max(real$lad + real$lad_sd, na.rm = TRUE) * 1.05

p <- ggplot(allp, aes(x = lad, y = height, colour = Cluster)) +
  geom_ribbon(data = real, aes(y = height, xmin = pmax(0, lad - lad_sd), xmax = lad + lad_sd, fill = Cluster),
              alpha = 0.22, colour = NA, inherit.aes = FALSE, orientation = "y") +
  geom_path(aes(linetype = type), linewidth = 1.3) +
  facet_wrap(~ Cluster, nrow = 1L) +
  scale_colour_manual(values = PAL_CLUSTER, guide = "none", drop = FALSE) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none", drop = FALSE) +
  scale_linetype_manual(values = c("Real (mean)" = "solid", "Uniform baseline" = "dashed"), name = NULL) +
  scale_x_continuous(limits = c(0, max_lad), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 42), breaks = c(0, 10, 20, 30, 40)) +
  labs(x = bquote("LAD" ~ "(m"^"2" ~ "m"^"-3" * ")"), y = "Height (m)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 22),
        axis.title = element_text(face = "bold", size = 22), axis.text = element_text(size = 16),
        legend.position = "bottom", legend.text = element_text(size = 18),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "fig_archetypes_profiles_uniform.png"), p, width = 18, height = 7.4, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_archetypes_profiles_uniform.pdf"), p, width = 18, height = 7.4, device = cairo_pdf)
cli_alert_success("Saved fig_archetypes_profiles_uniform (PNG 600dpi + PDF)")
