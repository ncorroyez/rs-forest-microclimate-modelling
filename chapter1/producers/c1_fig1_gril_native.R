# ==============================================================================
# Figure 1 — Gril-style: map (left) + per-archetype [side point cloud | LAD
# profile with Hmax/LAI/VCI] rows (right). Blois oak, P1..P4.

# Reads :
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
#         out_files/clusters/blois_clusters_P1P4_native20.tif
#         out_files/Chapter1/tables/typical_plots_per_cluster.csv
#         in_files/fig1_crosssection_cache.rds  (the LAS catalogue is read ONLY to rebuild it)
# Writes: outputs/figures_chap1/Fig1_typology_gril_native.png   (PNG only, plain ggsave)
#         in_files/fig1_crosssection_cache.rds   (only on a first run without the cache)
#   Rscript c1_fig1_gril_native.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(sf)
  library(png);library(grid);source("R/cluster_relabel.R");source("scripts/_article_style.R")})
dir.create("outputs/figures_chap1", showWarnings=FALSE, recursive=TRUE)
LV <- c("P1","P2","P3","P4")
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[, pid := sprintf("S%04d",.I)]; samp[, P := factor(relabel_cluster(Cluster), levels=LV)]
typ <- fread("out_files/Chapter1/tables/typical_plots_per_cluster.csv"); typ[,P:=factor(P,levels=LV)]

# ---- map (left) : ESRI satellite basemap + archetype dots + typical plots -----
suppressPackageStartupMessages(library(tidyterra))
clr <- terra::rast("out_files/clusters/blois_clusters_P1P4_native20.tif")
# Fill INTERIOR holes only (sub-canopy gaps/roads/water inside the massif, where a
# clustering trait was missing) by majority of neighbours, so the map is continuous;
# the exterior non-forest background is left untouched (forest outline unchanged).
na_m <- terra::ifel(is.na(clr), 1L, NA)
pat  <- terra::patches(na_m, directions=8)            # label connected NA regions
nr <- nrow(pat); nc <- ncol(pat); pv <- terra::values(pat)[,1]
bidx <- unique(c(1:nc, (nr-1)*nc + 1:nc, seq(1, by=nc, length.out=nr), seq(nc, by=nc, length.out=nr)))
ext_ids <- unique(stats::na.omit(pv[bidx]))           # NA patches touching the border = outside
hole <- terra::setValues(clr, as.integer(!is.na(pv) & !(pv %in% ext_ids)))  # 1 = interior hole
for (i in 1:10) {
  if (terra::global((hole==1) & is.na(clr), "sum", na.rm=TRUE)[[1]] == 0) break
  fm  <- terra::focal(clr, w=matrix(1,3,3), fun="modal", na.rm=TRUE, na.policy="only")
  clr <- terra::ifel((hole==1) & is.na(clr), fm, clr)
}
levels(clr) <- data.frame(id=1:4, archetype=c("P1","P2","P3","P4"))
xr <- c(terra::ext(clr)[1], terra::ext(clr)[2]); yr <- c(terra::ext(clr)[3], terra::ext(clr)[4])
pMap <- ggplot() +
  tidyterra::geom_spatraster(data=clr, aes(fill=archetype), maxcell=5e5) +
  geom_point(data=samp, aes(x,y), shape=21, fill="grey12", colour="white", size=1.5, stroke=0.5) +  # cLHS dots (dark w/ white ring)
  annotate("segment", x=xr[1]+0.04*diff(xr), xend=xr[1]+0.04*diff(xr)+1000, y=yr[1]+0.05*diff(yr), yend=yr[1]+0.05*diff(yr), linewidth=1.8, colour="grey15") +
  annotate("text", x=xr[1]+0.04*diff(xr)+500, y=yr[1]+0.05*diff(yr), label="1 km", vjust=-0.55, size=4.6, fontface="bold", colour="grey15") +
  annotate("segment", x=xr[2]-0.06*diff(xr), xend=xr[2]-0.06*diff(xr), y=yr[2]-0.15*diff(yr), yend=yr[2]-0.03*diff(yr), arrow=arrow(length=unit(4.5,"mm"),type="closed"), linewidth=1.6, colour="grey15") +
  annotate("text", x=xr[2]-0.06*diff(xr), y=yr[2]-0.01*diff(yr), label="N", fontface="bold", size=6, colour="grey15") +
  scale_fill_manual(values=PAL_CLUSTER, guide="none", na.translate=FALSE, na.value="transparent") +
  coord_sf(xlim=xr, ylim=yr, expand=FALSE, crs=32631, datum=32631) +
  labs(x=NULL, y=NULL) +
  theme_article(11) + theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank())

# ---- per-archetype rows: side cloud | LAD profile | metrics -------------------
# Cloud rendered as a 2D ggplot side view sharing the SAME height axis (0-40 m) as
# the profile, so the two are exactly aligned; coloured by height (vegetation ramp).
suppressPackageStartupMessages({library(ggtext)})
# REPRODUCIBILITY 2026-07-29. The four cross-section panels used to be clipped straight
# from the raw LAS catalogue on an EXTERNAL DRIVE, so Figure 1 could not be rebuilt
# unless /media/corroyez/MyPassport happened to be mounted. The clip is deterministic
# apart from the 30k subsample, so its OUTPUT is now cached inside the repo and the
# catalogue is read only to (re)build that cache. Delete the cache to force a re-clip.
LAS_DIR   <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/3-las_normalized_utm"
FIG1_CACHE <- "in_files/fig1_crosssection_cache.rds"
USE_CACHE  <- file.exists(FIG1_CACHE)
if (!USE_CACHE) {
  if (!dir.exists(LAS_DIR))
    stop("Figure 1: no cache at ", FIG1_CACHE, " and the LAS catalogue is not reachable at ",
         LAS_DIR, ".\n  Mount the drive once to build the cache, then the figure is self-contained.")
  suppressPackageStartupMessages(library(lidR))
  set.seed(42)                                   # the 30k point subsample below
  ctg <- readLAScatalog(LAS_DIR)
  sf::st_crs(ctg) <- 32631; lidR::opt_progress(ctg) <- FALSE
}
veg_ramp <- c("#6b4a2a","#3f5f2b","#4f8a2f","#79b23e","#bcd97e")   # trunk brown -> canopy green
# The representative plots are illustrative and predate the native cLHS sample, so
# rather than read a (mismatched) sample row, derive BOTH the cloud and its LAD
# profile from the SAME LAS clip (MacArthur-Horn, dz=1/k=0.5/z0=1, the pipeline recipe),
# guaranteeing cloud and profile describe the identical canopy.
prof_list <- if (USE_CACHE) readRDS(FIG1_CACHE) else lapply(LV, function(pcl){
  tr <- typ[as.character(P)==pcl][1]
  la <- tryCatch(clip_circle(ctg, tr$x, tr$y, 12), error=function(e)NULL)
  d <- data.table(z=numeric(), lad=numeric()); pc <- NULL
  if(!is.null(la)){
    dd <- as.data.table(la@data)[Z > -0.5 & Z < 40]
    lz <- as.data.table(lidR::LAD(dd$Z, dz=1, k=0.5, z0=1))[z <= tr$Hmax]
    if(nrow(lz) >= 4){ sm <- stats::spline(lz$z, lz$lad, n=160, method="natural")  # fluidify
                       d <- data.table(z=sm$x, lad=pmax(0, sm$y))[z <= tr$Hmax] }
    if(nrow(dd) > 30000) dd <- dd[sample(.N, 30000)]
    pc <- dd[,.(xr=X-tr$x, Z)]
  }
  list(pcl=pcl, tr=tr, d=d, pc=pc)
})
if (!USE_CACHE) { saveRDS(prof_list, FIG1_CACHE)
  cat(sprintf("Figure 1: cross-section cache written to %s (%.1f MB) — the drive is no longer needed\n",
              FIG1_CACHE, file.size(FIG1_CACHE)/1e6)) }
# ---- archetype MEAN LAD profiles (redesign 2026-08-28) ------------------------
# The profile panel used to show the LAD profile of the ONE typical plot whose
# point cloud sits beside it. It now shows the archetype MEAN profile: mean LAD at
# each 1 m layer across the archetype's 100 cLHS plots, with the interquartile
# range (25th-75th percentile) as a band. Above a plot's own canopy top its LAD is
# NA in the sample table; that is zero leaf area, so NA -> 0 before averaging.
# Truncation: each profile stops at the last layer still reached by at least half
# of the archetype's plots (Hmax >= layer centre), i.e. the archetype median Hmax;
# above that the mean would describe a shrinking minority of tall plots.
# The point cloud stays: it is ONE representative plot of the archetype, shown for
# shape only, whereas the curve is the archetype mean over its 100 plots. It carries
# no letter and is not located on the map: the two would invite the reader to read
# the mean curve as that single plot's profile. The annotated Hmax/LAI/VCI are now the
# ARCHETYPE MEANS over the same plots, not the example plot's values.
lad_cols <- grep("^LAD_Layer_", names(samp), value=TRUE)
lad_long <- melt(samp[, c("pid","P","Hmax", lad_cols), with=FALSE],
                 id.vars=c("pid","P","Hmax"), variable.name="layer", value.name="lad")
lad_long[, height := as.numeric(sub("LAD_Layer_", "", as.character(layer)))]
lad_long[is.na(lad), lad := 0]                       # no leaf area above the plot's top
prof_mean <- lad_long[, .(n_plots = .N, n_reach = sum(Hmax >= height),
                          lad_mean = mean(lad),
                          lad_q25 = as.numeric(quantile(lad, 0.25)),
                          lad_q75 = as.numeric(quantile(lad, 0.75))), by=.(P, height)]
prof_mean <- prof_mean[n_reach >= n_plots/2]         # truncate at the median Hmax
arch_stats <- samp[, .(n=.N, Hmax_m=mean(Hmax), LAI_m=mean(LAI), VCI_m=mean(VCI)), keyby=P]

XMAX <- max(prof_mean$lad_q75) * 1.05                # common PAD axis
rows <- lapply(prof_list, function(e){
  pcl <- e$pcl; tr <- e$tr
  dd <- prof_mean[P==pcl][order(height)]
  st <- arch_stats[P==pcl]
  # Floating point cloud: no axes/box visible, but the axis SPACE is kept (transparent
  # ink) so the panel matches the profile's and the two stay vertically aligned.
  cloud <- if(!is.null(e$pc)) ggplot(e$pc, aes(xr, Z, colour=Z)) + geom_point(size=0.14, alpha=0.5) +
      scale_colour_gradientn(colours=veg_ramp, guide="none") +
      coord_cartesian(ylim=c(0,40), xlim=c(-13,13)) + labs(x=expression(PAD~(m^2~m^-3)), y="Height (m)") +
      theme_article(10) + theme(axis.text=element_text(colour="transparent"),
        axis.title=element_text(colour="transparent"), axis.ticks=element_line(colour="transparent"),
        axis.line=element_blank(), panel.border=element_blank(), panel.background=element_blank(),
        plot.margin=margin(2,2,2,2)) else plot_spacer()
  # Metrics annotated INSIDE the profile (top-right corner, kept clear by the curve),
  # so the text column is dropped and cloud+profile can be enlarged.
  lab <- sprintf("**%s** (n = %d)<br>H<sub>max</sub> = %.1f m<br>LAI = %.1f<br>VCI = %.2f",
                 pcl, st$n, st$Hmax_m, st$LAI_m, st$VCI_m)
  gp <- ggplot(dd, aes(lad_mean, height)) +
    geom_ribbon(aes(y=height, xmin=lad_q25, xmax=lad_q75), orientation="y",
                fill=PAL_CLUSTER[[pcl]], alpha=0.25, colour=NA) +
    geom_path(colour=PAL_CLUSTER[[pcl]], linewidth=1.1) +
    ggtext::geom_richtext(data=data.frame(x=XMAX*0.99, y=39.5, lab=lab),
                     aes(x=x, y=y, label=lab), inherit.aes=FALSE, hjust=1, vjust=1, size=3.7,
                     fill=grDevices::adjustcolor("white",alpha=0.72), label.color=NA,
                     label.padding=unit(c(1,1,1,1),"pt"), colour=PAL_CLUSTER[[pcl]], lineheight=1.3) +
    labs(x=expression(PAD~(m^2~m^-3)), y=NULL) +
    coord_cartesian(ylim=c(0,40), xlim=c(0,XMAX)) +
    theme_article(10) + theme(plot.margin=margin(2,2,2,2))
  (cloud | gp) + plot_layout(widths=c(1.15, 1.0))
})
right <- rows[[1]] / rows[[2]] / rows[[3]] / rows[[4]]

# ---- bottom row: canopy-trait boxplots per archetype -------------------------
samp[, LAI1 := LAI]   # sample LAI is already one-sided (PAD); no /2
vn <- c(LAI1="LAI", Hmax="H<sub>max</sub> (m)", fCover="fCover", VCI="VCI")
#' Melt a trait table to long form for the raincloud row
#' @param dt table carrying P and the trait columns
#' @param cols named vector mapping column name to display label
#' @return long data.table(P, var, val) with non-finite values dropped
mk <- function(dt, cols){ m<-melt(dt[,c("P",names(cols)),with=F],id.vars="P",variable.name="var",value.name="val")
  m[,var:=factor(vn[as.character(var)],levels=vn)]; m[is.finite(val)] }
sl <- mk(samp, vn)
suppressPackageStartupMessages(library(ggdist))
# Clean coloured raincloud: filled half-violin (left) + jittered raw points (right),
# a subtle median tick. Coloured by archetype (P1..P4 red->green). The 400 cLHS
# design plots are the only sample shown: Chapter 1 rests on them alone.
pBox <- ggplot(sl, aes(P, val, fill=P, colour=P)) +
  ggdist::stat_slab(side="left", scale=0.62, justification=1.06, colour=NA, alpha=0.85,
                    normalize="groups", adjust=0.7) +
  geom_jitter(data=sl, aes(P,val,colour=P), inherit.aes=FALSE, size=0.5, alpha=0.45,
              position=position_jitter(width=0.07, height=0, seed=1)) +
  stat_summary(fun=median, geom="point", colour="grey15", fill="white", shape=21, size=1.9, stroke=0.6) +
  facet_wrap(~var, scales="free_y", nrow=1) +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") + scale_colour_manual(values=PAL_CLUSTER, guide="none") +
  coord_cartesian(clip="off") +
  labs(x=NULL, y=NULL, caption="Half-violin + colored points: 400 cLHS design plots, 100 per archetype") +
  theme_article(12) + theme(strip.text=ggtext::element_markdown(face="bold", size=12.5, margin=margin(b=4)),
    strip.background=element_blank(), panel.spacing=unit(1.4,"lines"),
    plot.caption=element_text(size=10.5, colour="grey25", hjust=0.5, margin=margin(t=6)),
    axis.line=element_line(colour="grey70", linewidth=0.3), panel.border=element_blank())

fig <- (pMap | patchwork::wrap_elements(right)) / pBox +
  plot_layout(heights=c(4, 1)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=17), plot.tag.position=c(0.01,0.99))
ggsave("outputs/figures_chap1/Fig1_typology_gril_native.png", fig, width=13.5, height=13, dpi=350, bg="white")
cat("DONE -> outputs/figures_chap1/Fig1_typology_gril_native.png\n")
