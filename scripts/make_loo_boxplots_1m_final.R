# ==============================================================================
# FINAL LOO Δv boxplots at FIXED 1 m (HOBO height) — replaces the nair==1
# versions for absolute ΔTmax / ΔVPDmax. Reads tab_loo_deltav_height_compare.csv
# (already holds both heights) and keeps the 1 m subset.
#   Δv<0 = trait lowers the daily max => buffers.
# Overwrites:
#   fig_loo_boxplot_tmax_vpd.{png,pdf}        (pooled, metric × period)
#   fig_loo_boxplot_tmax_by_cluster.{png,pdf} (period × cluster)
#   fig_loo_boxplot_vpd_by_cluster.{png,pdf}
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(tidyverse); library(patchwork)
})
OUT  <- here::here("outputs/figs_MEB2026_final")
long <- fread(file.path(OUT, "tab_loo_deltav_height_compare.csv"))
long <- long[height == "1 m (interp., HOBO)"]          # <-- fixed 1 m only
long[, period := factor(period, levels = c("All period","10% hottest"))]
long[, Cluster := factor(Cluster, levels = paste0("P",1:4))]

trait_lab  <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
ord <- c("LAI","fCover","LAD","Hmax")
long[, trait_lab := factor(trait_lab[trait], levels = trait_lab[ord])]
trait_fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")
CAP1M <- "Extracted at fixed 1 m (HOBO height). "

# ---- (1) pooled: metric × period -------------------------------------------
ylim_of <- function(d){ q<-quantile(d,c(0.02,0.98),na.rm=TRUE); p<-diff(q)*0.12; c(q[1]-p,q[2]+p) }
pooled_panel <- function(metric_key, ylab_expr){
  d <- long[metric == metric_key]; yl <- ylim_of(d$delta)
  nc <- sum(d$delta < yl[1] | d$delta > yl[2], na.rm=TRUE)
  ggplot(d, aes(trait_lab, delta, fill = trait)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85) +
    geom_jitter(width=0.12, size=1.0, alpha=0.35, colour="grey20") +
    facet_wrap(~ period) + coord_cartesian(ylim = yl) +
    scale_x_discrete(labels=function(x) parse(text=x)) +
    scale_fill_manual(values=trait_fill) +
    labs(x=NULL, y=ylab_expr,
         subtitle = if (nc>0) sprintf("(%d sparse-plot amplification points clipped)", nc) else NULL) +
    theme_bw(base_size=18) +
    theme(panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
          axis.text.x=element_text(face="bold", size=16),
          plot.subtitle=element_text(size=10, colour="grey45"), legend.position="none")
}
pT <- pooled_panel("Delta*T[max]~(degree*C)", expression(Delta[v]^"LOO" ~ "on " * Delta*T[max] ~ "(" * degree * "C)"))
pV <- pooled_panel("Delta*VPD[max]~(kPa)",    expression(Delta[v]^"LOO" ~ "on " * Delta*VPD[max] ~ "(kPa)"))
p_pool <- (pT / pV) + plot_annotation(
  caption = paste0(CAP1M, "Δv<0 = trait lowers the daily max (cools / dries-down) => buffers. fCover dominant & strengthens on hot days."),
  theme = theme(plot.caption=element_text(size=11, colour="grey40")))
ggsave(file.path(OUT,"fig_loo_boxplot_tmax_vpd.png"), p_pool, width=13, height=10, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_loo_boxplot_tmax_vpd.pdf"), p_pool, width=13, height=10, device=cairo_pdf)
cli_alert_success("Saved fig_loo_boxplot_tmax_vpd (1 m)")

# ---- (2) by-cluster: period × cluster, one per metric ----------------------
by_cluster <- function(metric_key, ylab_expr, fname){
  d <- long[metric == metric_key]
  q <- quantile(d$delta, c(0.04,0.96), na.rm=TRUE); pad <- diff(q)*0.10
  yl <- c(q[1]-pad, q[2]+pad); nc <- sum(d$delta < yl[1] | d$delta > yl[2], na.rm=TRUE)
  p <- ggplot(d, aes(trait_lab, delta, fill = trait)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_boxplot(width=0.62, outlier.shape=NA, alpha=0.85) +
    geom_jitter(width=0.12, size=0.9, alpha=0.4, colour="grey20") +
    facet_grid(period ~ Cluster) + coord_cartesian(ylim=yl) +
    scale_x_discrete(labels=function(x) parse(text=x)) +
    scale_fill_manual(values=trait_fill) +
    labs(x=NULL, y=ylab_expr,
         caption=sprintf("%sΔv<0 = trait buffers. P1 sparse → P4 dense.%s", CAP1M,
                         if (nc>0) sprintf(" %d amplification points clipped.", nc) else "")) +
    theme_bw(base_size=18) +
    theme(panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold", size=16),
          axis.text.x=element_text(size=14, angle=30, hjust=1),
          plot.caption=element_text(size=11, colour="grey40"), legend.position="none")
  ggsave(file.path(OUT, paste0(fname,".png")), p, width=16, height=8.5, dpi=300, bg="white")
  ggsave(file.path(OUT, paste0(fname,".pdf")), p, width=16, height=8.5, device=cairo_pdf)
  cli_alert_success("Saved {fname} (1 m)")
}
by_cluster("Delta*T[max]~(degree*C)", expression(Delta[v]^"LOO" ~ "on " * Delta*T[max] ~ "(" * degree * "C)"), "fig_loo_boxplot_tmax_by_cluster")
by_cluster("Delta*VPD[max]~(kPa)",    expression(Delta[v]^"LOO" ~ "on " * Delta*VPD[max] ~ "(kPa)"), "fig_loo_boxplot_vpd_by_cluster")

cli_h2("1 m median Δv per metric × period × trait (pooled)")
print(long[, .(median=round(median(delta,na.rm=TRUE),3)), by=.(metric,period,trait)][order(metric,period,median)])
