# ==============================================================================
# By-cluster (P1->P4) version of the LOO Δv boxplots on ΔTmax and ΔVPDmax @1m,
# for whole period vs 10% hottest days. Reads the precomputed per-plot table.
# Two figures (one per metric); facet_grid(period ~ cluster).
#   Δv<0 = trait lowers the daily max (cools / dries-down) => buffers.
# Outputs: fig_loo_boxplot_tmax_by_cluster.{png,pdf}
#          fig_loo_boxplot_vpd_by_cluster.{png,pdf}
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(tidyverse); library(patchwork)
})
OUT  <- here::here("outputs/figs_MEB2026_final")
long <- fread(file.path(OUT, "tab_loo_deltav_tmax_vpd.csv"))

long[, period := factor(period, levels = c("All period","10% hottest days"))]
long[, Cluster := factor(Cluster, levels = paste0("P", 1:4))]
trait_lab <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
ord <- c("LAI","fCover","LAD","Hmax")
long[, trait_lab := factor(trait_lab[trait], levels = trait_lab[ord])]
trait_fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")

make_fig <- function(metric_key, ylab_expr, fname) {
  d  <- long[metric == metric_key]
  # shared clipped y across clusters (keeps magnitudes comparable; annotate clip)
  q  <- quantile(d$delta, c(0.04, 0.96), na.rm = TRUE)
  pad <- diff(q) * 0.10; yl <- c(q[1] - pad, q[2] + pad)
  nc <- sum(d$delta < yl[1] | d$delta > yl[2], na.rm = TRUE)
  p <- ggplot(d, aes(trait_lab, delta, fill = trait)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.4, colour = "grey20") +
    facet_grid(period ~ Cluster) +
    coord_cartesian(ylim = yl) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_fill_manual(values = trait_fill) +
    labs(x = NULL, y = ylab_expr,
         caption = sprintf("1 m (nair==1). Δv<0 = trait buffers. P1 sparse → P4 dense.%s",
                           if (nc > 0) sprintf(" %d amplification points clipped.", nc) else "")) +
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold", size = 17),
          axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
          plot.caption = element_text(size = 11, colour = "grey40"),
          legend.position = "none")
  ggsave(file.path(OUT, paste0(fname, ".png")), p, width = 16, height = 8.5, dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(fname, ".pdf")), p, width = 16, height = 8.5, device = cairo_pdf)
  cli_alert_success("Saved {fname} ; clipped={nc} ; ylim=[{round(yl[1],2)},{round(yl[2],2)}]")
}

cli_h1("By-cluster LOO boxplots: ΔTmax and ΔVPDmax")
make_fig("Delta*T[max]~(degree*C)",
         expression(Delta[v]^"LOO" ~ "on " * Delta*T[max] ~ "(" * degree * "C)"),
         "fig_loo_boxplot_tmax_by_cluster")
make_fig("Delta*VPD[max]~(kPa)",
         expression(Delta[v]^"LOO" ~ "on " * Delta*VPD[max] ~ "(kPa)"),
         "fig_loo_boxplot_vpd_by_cluster")

# medians table per cluster
cli_h2("Median Δv per cluster × metric × period × trait")
print(long[, .(median = round(median(delta, na.rm=TRUE),3)),
           by=.(metric, period, Cluster, trait)][order(metric, period, Cluster, median)])
