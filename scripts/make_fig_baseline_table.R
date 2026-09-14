# ==============================================================================
# fig_baseline_table — which baseline VALUE is used for each variable, per
# archetype P1..P4 (MAIN, per-cluster baseline = that archetype's mean) and for
# the GLOBAL baseline (complementary = landscape mean). Exact values as fed to
# MuSICA (c3_shapley_chunk.R::bl). LAD baseline is the vertically uniform profile.
#   Rscript scripts/make_fig_baseline_table.R
# Out: outputs/figs_MEB2026_final/fig_baseline_table.{png,pdf}
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
OUT <- here::here("outputs/figs_MEB2026_final"); if (!dir.exists(OUT)) OUT <- "outputs/figs_MEB2026_final"
df <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
df[, P := relabel_cluster(Cluster)]

# Baseline values per the pipeline. LAI DISPLAYED one-sided = stored/2 (the stored
# sample LAI is the to-be-doubled value MuSICA gets; Fig 1 shows the same /2).
blf <- function(d) data.table(LAI = mean(d$LAI, na.rm=TRUE) / 2,
                              Hmax = round(mean(floor(d$Hmax)+1, na.rm=TRUE)),
                              fCover = mean(d$fCover, na.rm=TRUE))
percl <- df[, blf(.SD), by = P][order(P)]
glob  <- cbind(P = "Global", blf(df))
tab   <- rbind(percl, glob)

m <- melt(tab, id.vars = "P", variable.name = "variable", value.name = "value")
m <- rbind(m, data.table(P = tab$P, variable = "LAD", value = NA_real_))   # LAD = uniform
m[, P := factor(P, levels = c("P1","P2","P3","P4","Global"))]
m[, variable := factor(variable, levels = c("LAD","fCover","Hmax","LAI"))] # LAI on top
m[, label := fifelse(variable == "LAD", "uniform",
              fifelse(variable == "Hmax", sprintf("%.0f m", value), sprintf("%.2f", value)))]
m[, fillv := (value - min(value, na.rm=TRUE)) / (max(value, na.rm=TRUE) - min(value, na.rm=TRUE) + 1e-9), by = variable]
m[, txt := fifelse(is.na(fillv), "grey15", fifelse(fillv > 0.55, "white", "grey15"))]

p <- ggplot(m, aes(P, variable)) +
  geom_tile(aes(fill = fillv), colour = "white", linewidth = 1.2) +
  geom_vline(xintercept = 4.5, linewidth = 1.1, colour = "grey40") +     # separate Global
  geom_text(aes(label = label, colour = txt), size = 4.6, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_viridis_c(option = "mako", direction = -1, na.value = "grey88", guide = "none") +
  scale_x_discrete(position = "top") +
  annotate("text", x = 2.5, y = 4.9, label = "per-cluster baseline (MAIN)", size = 4, fontface = "italic", colour = "grey25") +
  annotate("text", x = 5,   y = 4.9, label = "complementary", size = 3.6, fontface = "italic", colour = "grey25") +
  coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
  labs(x = NULL, y = NULL,
       caption = "LAI one-sided; Hmax in m; fCover fractional. LAD baseline = vertically uniform profile (Fig. 1c).") +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(), axis.text = element_text(face = "bold", size = 14),
        plot.caption = element_text(size = 9, colour = "grey40", hjust = 0))
ggsave(file.path(OUT, "fig_baseline_table.png"), p, width = 8.8, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_baseline_table.pdf"), p, width = 8.8, height = 4.6, device = cairo_pdf)
cat("=== baseline values ===\n"); print(tab); cat("DONE -> fig_baseline_table.png\n")
