# ==============================================================================
# Merge + figures for the 6-metric per-cluster cLHS Shapley (c1_metrics_chunk.R).
# Metrics × periods: ΔTmax, micro/macro slope, ΔVPDmax × {all period, hot days}.
# Per-cluster baseline → importance = mean|φ| (raw, with units) AND relative share
# (|φ_trait| / Σ|φ|) so metrics with different units are comparable on one scale.
# Outputs: out_files/Chapter1/figures/FigSh_metrics_heatmap.png  (+ tables csv)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
pdir <- "out_files/Chapter1/tables/metrics_parts"
S <- rbindlist(lapply(list.files(pdir, "part_.*csv$", full.names=TRUE), fread))
S[, Cluster := relabel_cluster(Cluster)]
S[, c("var","period") := tstrsplit(metric, "_")]
S[, var := factor(var, levels=c("Tmax","slope","VPD"), labels=c("ΔTmax","micro/macro slope","ΔVPDmax"))]
S[, period := factor(period, levels=c("all","hot"), labels=c("all period","hot days"))]
Fv <- c("LAI","Hmax","fCover","LAD")

cat(sprintf("merged %d rows -> %d plots × 6 metric-periods (%s/cluster)\n",
            nrow(S), uniqueN(S$pid), paste(S[metric=="Tmax_all", .N, by=Cluster][order(Cluster)]$N, collapse="/")))

# additivity
S[, sumphi := LAI+Hmax+fCover+LAD][, fmb := full-base]
cat(sprintf("additivity max|Σφ-(full-base)| = %.2e  (must be ~0)\n", max(abs(S$sumphi-S$fmb), na.rm=TRUE)))

# mean |φ| per cluster × metric × trait (raw, with units)
imp <- melt(S, id.vars=c("Cluster","var","period"), measure.vars=Fv, variable.name="trait", value.name="phi")
agg <- imp[, .(absphi = mean(abs(phi), na.rm=TRUE)), by=.(Cluster, var, period, trait)]
agg[, share := absphi / sum(absphi), by=.(Cluster, var, period)]   # relative importance (rows sum to 1)
fwrite(dcast(agg, Cluster+var+period+trait ~ ., value.var=c("absphi","share")),
       "out_files/Chapter1/tables/Table_metrics_importance.csv")

cat("\n=== relative importance share (%) — leading trait per cluster × metric × period ===\n")
lead <- agg[, .SD[which.max(share)], by=.(Cluster, var, period)][order(var, period, Cluster)]
print(lead[, .(Cluster, var, period, trait, share=round(100*share))], nrow=99)

# --- heatmap: relative importance, x=trait, y=metric×period, facet=cluster -----
agg[, mp := factor(paste(var, period, sep=" — "),
      levels=rev(c("ΔTmax — all period","ΔTmax — hot days","micro/macro slope — all period",
                   "micro/macro slope — hot days","ΔVPDmax — all period","ΔVPDmax — hot days")))]
agg[, trait := factor(trait, levels=Fv)]
g <- ggplot(agg, aes(trait, mp, fill=share)) +
  geom_tile(colour="white", linewidth=0.6) +
  geom_text(aes(label=sprintf("%.0f", 100*share)), size=3) +
  facet_wrap(~ Cluster, nrow=1) +
  scale_fill_viridis_c(option="magma", direction=-1, name="relative\nimportance\n|φ| share (%)",
                       labels=function(x) round(100*x), limits=c(0, max(agg$share))) +
  labs(title="Within-archetype trait importance across metrics (per-point cLHS Shapley, per-cluster baseline)",
       subtitle="Share of total mean|φ| per trait. Density-dependent shift LAI→LAD from open (P1) to dense (P4) holds across ΔTmax, slope and ΔVPDmax, all period and hot days.",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(), strip.text=element_text(face="bold"),
        plot.subtitle=element_text(size=8.5), axis.text.x=element_text(angle=0))
figdir <- "out_files/Chapter1/figures"; dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(figdir,"FigSh_metrics_heatmap.png"), g, width=12, height=5.2, dpi=300, bg="white")

# --- 2nd version: ABSOLUTE mean|φ| in native units (free y per metric) --------
# Instead of the relative share, show how much each trait actually moves each
# metric (°C / dimensionless slope / kPa). Bars split all-period vs hot days.
agg[, var := factor(var, levels=c("ΔTmax","micro/macro slope","ΔVPDmax"))]
g2 <- ggplot(agg, aes(trait, absphi, fill=period)) +
  geom_col(position=position_dodge(width=0.72), width=0.66) +
  facet_grid(var ~ Cluster, scales="free_y", switch="y") +
  scale_fill_manual(values=c("all period"="#4575B4","hot days"="#D73027"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(title="Within-archetype trait effect on each metric (absolute mean|φ|, native units)",
       subtitle="Mean|φ| in native units (°C / dimensionless slope / kPa; free y per metric). The LAI→LAD shift from open (P1) to dense (P4) holds in magnitude, and hot days (red) sharpen every effect.",
       x=NULL, y=expression("mean |"*phi*"|  (native units of each metric)")) +
  theme_bw(base_size=11) +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"), strip.placement="outside",
        plot.subtitle=element_text(size=8), legend.position="bottom")
ggsave(file.path(figdir,"FigSh_metrics_absolute.png"), g2, width=12, height=6.4, dpi=300, bg="white")
cat("DONE -> out_files/Chapter1/figures/FigSh_metrics_absolute.png\n")

# --- raw mean|φ| table (with units) per metric×period -------------------------
cat("\n=== mean|φ| (raw units: °C for ΔTmax, dimensionless for slope, kPa for ΔVPDmax) by cluster ===\n")
raw <- dcast(agg, var+period+Cluster ~ trait, value.var="absphi")
print(raw[order(var,period,Cluster)], digits=2, nrow=99)
cat("\nDONE -> out_files/Chapter1/figures/FigSh_metrics_heatmap.png\n")
