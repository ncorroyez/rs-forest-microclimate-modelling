# Merge the 8 Shapley chunk parts -> aggregate per cluster + figure.
# Run from z_Example root:  Rscript c3_shapley_merge.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("Chapter3_config.R") })
pdir <- file.path(CFG_C3$out_dir, "tables", "shapley_parts")
parts <- list.files(pdir, pattern="part_.*\\.csv$", full.names=TRUE)
S <- rbindlist(lapply(parts, fread), fill=TRUE)
S <- S[is.finite(dTmax_full) & is.finite(dTmax_base)]
cat(sprintf("merged %d parts -> %d plots (%s/cluster)\n", length(parts), nrow(S),
            paste(S[,.N,by=Cluster][order(Cluster)]$N, collapse="/")))
Fv <- c("LAI","Hmax","fCover","LAD")
fwrite(S, file.path(CFG_C3$out_dir, "tables", "tab_shapley_clhs_perplot.csv"))

cat("\n=== mean Shapley (°C on ΔTmax, vs cluster baseline) per cluster ===\n")
print(S[, c(.(n=.N), lapply(.SD, function(x) round(mean(x,na.rm=T),3))), by=Cluster, .SDcols=Fv][order(Cluster)])
cat("\n=== within-cluster IMPORTANCE = mean |Shapley| per factor ===\n")
print(S[, lapply(.SD, function(x) round(mean(abs(x),na.rm=T),3)), by=Cluster, .SDcols=Fv][order(Cluster)])
cat("\n=== additivity (sum vs full-base) ===\n")
print(S[, .(sum=round(mean(LAI+Hmax+fCover+LAD,na.rm=T),3),
            fmb=round(mean(dTmax_full-dTmax_base,na.rm=T),3)), by=Cluster][order(Cluster)])

ml <- melt(S, id.vars=c("pid","Cluster"), measure.vars=Fv, variable.name="factor", value.name="shapley")
ml[, Cluster := paste0("Cluster ", Cluster)]
g <- ggplot(ml, aes(factor, shapley, fill=factor)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60") +
  geom_boxplot(outlier.size=0.4, alpha=0.85) + facet_wrap(~Cluster, nrow=1) +
  scale_fill_brewer(palette="Set2", guide="none") +
  labs(title="Per-point cLHS Shapley: drivers of within-cluster ΔTmax variation (2×LAI)",
       subtitle="100 cLHS plots/cluster; per-cluster baseline (each factor at its cluster mean); LAD baseline = uniform.\nEach plot decomposed vs its cluster archetype; boxes = within-cluster distribution of factor contributions.",
       x=NULL, y=expression(Shapley~(degree*C~"on"~Delta*T[max]))) +
  theme_minimal(base_size=11) + theme(plot.subtitle=element_text(size=8))
ggsave(file.path(CFG_C3$out_dir, "figures", "FigSh_shapley_clhs_boxplot_supp.png"), g, width=11, height=4.2, dpi=150)  # supplementary boxplot view; canonical Fig 2 = pipeline/11_clhs_attribution.R
cat("\nDONE merge\n")
