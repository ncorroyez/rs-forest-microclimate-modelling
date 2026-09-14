# ==============================================================================
# Do the 53 loggers cover the trait space the 400-plot design explores? Design cloud with
# the loggers overlaid. Loggers carry their native one-sided LAI, like the design.
#
# Reads : out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
#         in_files/data_Blois_utm31n.geojson / in_files_native20/*.tif
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv       (stage A2)
# Writes: out_files/Chapter1/figures/FigSh_plot_trait_coverage_native.png    (PNG only)
#   Rscript scripts/c1_hobo_coverage_native.R
# ==============================================================================
# D1 coverage — NATIVE, consistent ONE-SIDED LAI (the rds LAI is one-sided; run_musica_one
# doubles it for MuSICA, R/musica.R:69). The old script wrongly /2'd the cLHS ("two-sided")
# while leaving the loggers un-halved, putting them on different scales. Here both cLHS and
# loggers use their native one-sided LAI directly.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(terra); library(sf); source("R/cluster_relabel.R"); source("scripts/_article_style.R") })
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[, P := relabel_cluster(Cluster)]; samp[, LAI1 := LAI]                 # rds LAI is one-sided (no /2)
# native logger traits (one-sided), restricted to the 53 loggers with observations
g <- st_read("in_files/data_Blois_utm31n.geojson", quiet=TRUE); r0 <- rast("in_files_native20/lai_z1_res_10_m.tif"); g <- st_transform(g, crs(r0)); v <- vect(g)
H <- data.table(id_plot=g$id_plot, LAI1=terra::extract(r0,v)[,2], Hmax=terra::extract(rast("in_files_native20/max_res_10_m.tif"),v)[,2], fCover=terra::extract(rast("in_files_native20/fCover_res_10_m.tif"),v)[,2])
obs <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot)]
H <- H[id_plot %in% obs$id_plot & is.finite(LAI1)]
PAL <- PAL_CLUSTER   # house palette P1 red -> P4 green, as every other figure
#' One trait-space coverage panel: design cloud behind, the 53 loggers on top
#' @param xv name of the x trait column
#' @param xl x axis label
#' @param yv name of the y trait column
#' @param yl y axis label
#' @return a ggplot object
base <- function(xv,xl,yv,yl) ggplot() +
  geom_point(data=samp, aes(.data[[xv]], .data[[yv]], colour=P), alpha=.18, size=.9) +
  geom_point(data=H, aes(.data[[xv]], .data[[yv]]), shape=21, fill="black", colour="white", size=2.2, stroke=.3) +
  scale_colour_manual(values=PAL, name="Archetype") + labs(x=xl,y=yl) +
  theme_article(11)
pA <- base("LAI1","LAI (one-sided)","Hmax",expression(H[max]~(m))) 
pB <- base("LAI1","LAI (one-sided)","fCover","fCover") 
pC <- base("fCover","fCover","Hmax",expression(H[max]~(m))) 
nd <- nrow(H[LAI1>=4]); p4 <- round(mean(samp[P=="P4",LAI1]),1)
sub <- sprintf("%d sub-canopy plots (black); %d have one-sided LAI ≥ 4, the dense archetype P4 averaging %.1f. Faint points: the 400 cLHS design plots.", nrow(H), nd, p4)
fig <- (pA+pB+pC) + plot_layout(guides="collect") +
  plot_annotation(tag_levels=list(c("(a)","(b)","(c)")), subtitle=sub,
    theme=theme(plot.subtitle=element_text(size=8.5,colour="grey35"), plot.title=element_text(face="bold"))) & theme(legend.position="bottom")
ggsave("out_files/Chapter1/figures/FigSh_plot_trait_coverage_native.png", fig, width=15, height=4.8, dpi=200, bg="white")
cat(sprintf("DONE native coverage | n_loggers=%d, #one-sided LAI>=4 = %d, cLHS-P4 mean = %.1f\n", nrow(H), nd, p4))
