# ==============================================================================
# PIPELINE STAGE 06 — Shapley boxplots (pooled + per cluster), ΔTmax & ΔVPDmax.
# Outputs (OUT_FIG): fig_shapley_pooled.{png,pdf},
#   fig_shapley_tmax_by_cluster.{png,pdf}, fig_shapley_vpd_by_cluster.{png,pdf}
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
suppressMessages(library(patchwork))
cli_h1("STAGE 06 — Shapley boxplots")

phi <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "shapley_perplot.rds")))
phi[, Cluster := factor(Cluster, levels=paste0("P",1:4))]
mlabs <- c(Tmax="Delta*T[max]~(degree*C)", VPD="Delta*VPD[max]~(kPa)")
phi[, metric_lab := mlabs[metric]]
tl <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD"); ordv <- c("LAI","fCover","LAD","Hmax")
phi[, trait_lab := factor(tl[trait], levels=tl[ordv])]
fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")
yl_of <- function(d,q=c(.04,.96)){ z<-quantile(d,q,na.rm=TRUE); p<-diff(z)*.10; c(z[1]-p,z[2]+p) }

# SINGLE by-cluster figure with BOTH metrics: rows = metric×period (units differ →
# free-y per row, range set by geom_blank), cols = clusters P1..P4.
row_lvls <- c("Delta*T[max]*' (all, '*degree*'C)'", "Delta*T[max]*' (hot, '*degree*'C)'",
              "Delta*VPD[max]*' (all, kPa)'",       "Delta*VPD[max]*' (hot, kPa)'")
phi[, rowlab := fcase(
  metric=="Tmax" & period=="All period",  row_lvls[1],
  metric=="Tmax" & period=="10% hottest", row_lvls[2],
  metric=="VPD"  & period=="All period",  row_lvls[3],
  metric=="VPD"  & period=="10% hottest", row_lvls[4])]
phi[, rowlab := factor(rowlab, levels=row_lvls)]

# add an "All" column = all plots pooled, alongside P1..P4
phiC <- rbind(phi, copy(phi)[, Cluster := "All"])
clusters <- c(paste0("P",1:4), "All")
phiC[, Cluster := factor(Cluster, levels=clusters)]

rrng <- phiC[, { yl <- yl_of(phi, c(.03,.97)); .(ymin=yl[1], ymax=yl[2]) }, by=rowlab]
blankC <- rbindlist(lapply(clusters, function(cc)
  data.table(rowlab=rep(rrng$rowlab,2), phi=c(rrng$ymin,rrng$ymax), Cluster=cc)))
blankC[, `:=`(trait_lab=factor(tl["LAI"],levels=tl[ordv]),
              Cluster=factor(Cluster,levels=clusters), rowlab=factor(rowlab,levels=row_lvls))]
jitC <- merge(phiC, rrng, by="rowlab")[phi>=ymin & phi<=ymax]
nclipC <- nrow(phiC) - nrow(jitC)

p_clu <- ggplot(phiC, aes(trait_lab, phi, fill=trait)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_blank(data=blankC, aes(trait_lab, phi), inherit.aes=FALSE) +
  geom_boxplot(width=0.62, outlier.shape=NA, alpha=0.85) +
  geom_jitter(data=jitC, width=0.12, size=0.7, alpha=0.35, colour="grey20") +
  facet_grid(rowlab ~ Cluster, scales="free_y", labeller=labeller(rowlab=label_parsed)) +
  scale_x_discrete(labels=function(x) parse(text=x)) + scale_fill_manual(values=fill) +
  labs(x=NULL, y=expression(varphi[v]^"Shapley"~"(per-plot exact Shapley)"),
       caption=sprintf("Exact Shapley, 1 m, fCover baseline=mean. φ<0 = trait buffers. P1 sparse → P4 dense. %d amplification points clipped.", nclipC)) +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=13), axis.text.x=element_text(size=13,angle=30,hjust=1),
        plot.caption=element_text(size=10.5,colour="grey40"), legend.position="none")
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_by_cluster.png"), p_clu, width=18, height=11, dpi=300, bg="white")
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_by_cluster.pdf"), p_clu, width=18, height=11, device=cairo_pdf)
cli_alert_success("Saved fig_shapley_by_cluster (both metrics, P1-P4 + All, single figure)")

# SINGLE figure: facet_grid(metric × period), per-metric y-range via geom_blank
# (free_y), boxplots without outliers + jitter clipped to the per-metric range.
rng <- phi[, { yl <- yl_of(phi, c(.02,.98)); .(ymin=yl[1], ymax=yl[2]) }, by=.(metric, metric_lab)]
periods <- levels(phi$period)
blank_df <- rbindlist(lapply(periods, function(pp)
  data.table(metric_lab = rep(rng$metric_lab, 2),
             phi        = c(rng$ymin, rng$ymax),
             period     = pp)))
blank_df[, `:=`(trait_lab = factor(tl["LAI"], levels=tl[ordv]), period=factor(period, levels=periods))]
jit <- merge(phi, rng[, .(metric, ymin, ymax)], by="metric")[phi >= ymin & phi <= ymax]
nclip <- nrow(phi) - nrow(jit)

p_pool <- ggplot(phi, aes(trait_lab, phi, fill=trait)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_blank(data=blank_df, aes(trait_lab, phi), inherit.aes=FALSE) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85) +
  geom_jitter(data=jit, width=0.12, size=0.9, alpha=0.32, colour="grey20") +
  facet_grid(metric_lab ~ period, scales="free_y", labeller=labeller(metric_lab=label_parsed)) +
  scale_x_discrete(labels=function(x) parse(text=x)) + scale_fill_manual(values=fill) +
  labs(x=NULL, y=expression(varphi[v]^"Shapley"~"(per-plot exact Shapley)"),
       caption=sprintf("Exact Shapley, 1 m, fCover baseline=mean. φ<0 = trait buffers. Σφ = REF − baseline. %d amplification points clipped.", nclip)) +
  theme_bw(base_size=17) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=15), axis.text.x=element_text(face="bold", size=15),
        plot.caption=element_text(size=10.5, colour="grey40"), legend.position="none")
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_pooled.png"), p_pool, width=12, height=9, dpi=300, bg="white")
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_pooled.pdf"), p_pool, width=12, height=9, device=cairo_pdf)
cli_alert_success("Stage 06 done.")
