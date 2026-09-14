# Fig 3 (attribution) — CLEAN native version matching Table H1 exactly.
# The old c1_deltadelta_perplot.R double-normalised (×trait-SD on top of the
# already-per-SD metrics6 values), giving figures ~3.4x too small (P4 LAI −0.084
# instead of −0.288). Here we plot the SAME per-plot sensitivities the bootstrap
# (c1_attribution_bootstrap.R) aggregates for Table H1: per +1 within-archetype SD.
# Violins = per-plot distribution per archetype; white boxes = median + IQR.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("scripts/_article_style.R"); source("R/cluster_relabel.R") })
LV <- c("P1","P2","P3","P4")
M <- rbindlist(lapply(list.files("out_files/Chapter1/tables/metrics6_native20","part_.*csv$",full.names=TRUE), fread), fill=TRUE)
vci <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))[, .(pid=sprintf("S%04d",.I), VCI)]
M <- merge(M, vci, by="pid", all.x=TRUE)
M[, P := factor(relabel_cluster(Cluster), levels=LV)]
# per-plot sensitivities, EXACTLY as the bootstrap / Table H1
M[, `:=`(LAI=(LAI_add-LAI_rem)/2, fCover=(fCov_add-fCov_rem)/2, Hmax=(Hmax_add-Hmax_rem)/2, LAD=dT_LAD)]
# profile lever = RAW real-vs-uniform ΔTmax contrast (dT_LAD, set on the line above): a FULL
# real→uniform swap at fixed LAI/Hmax/fCover — an UPPER BOUND, not a ±1 SD move. NO (1-VCI)
# rescaling (committee 2026-07: a given VCI maps to infinitely many profiles, and this is not a
# per-SD quantity). The scalar levers stay per +1 within-archetype SD; footings differ by design.
# slope metric expressed on log(slope): Δlog(β) ≈ Δβ/β_base
issl <- grepl("slope", M$metric)
for (cc in c("LAI","fCover","Hmax","LAD")) M[issl & is.finite(base) & base>0, (cc) := get(cc)/base]
PERIOD <- Sys.getenv("PERIOD","all")
keep <- if (PERIOD=="hot") c(Tmax_hot="Delta*T[max]~(hot)", slope_hot="Delta*log(slope)~(hot)") else c(Tmax_all="Delta*T[max]", slope_all="Delta*log(slope)")
D <- M[metric %in% names(keep)]
L <- melt(D, id.vars=c("pid","P","metric"), measure.vars=c("LAI","fCover","Hmax","LAD"), variable.name="lever", value.name="s")
L[, lever := factor(lever, levels=c("LAI","fCover","Hmax","LAD"), labels=c("LAI","fCover","H[max]","vertical~profile"))]
L[, metriclab := factor(keep[metric], levels=unname(keep))]
cat("check native P4 LAI ΔTmax median (should be −0.288):", round(median(L[metric=="Tmax_all"&P=="P4"&lever=="LAI"]$s,na.rm=TRUE),3), "\n")
p <- ggplot(L[is.finite(s)], aes(lever, s, fill=P)) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.4) +
  geom_violin(aes(group=interaction(lever,P)), scale="width", width=0.85, colour=NA, alpha=0.75, position=position_dodge(width=0.9)) +
  geom_boxplot(aes(group=interaction(lever,P)), width=0.14, fill="white", outlier.size=0.25, position=position_dodge(width=0.9), linewidth=0.3) +
  facet_grid(metriclab ~ ., scales="free_y", labeller=label_parsed) +
  scale_fill_manual(values=PAL_CLUSTER, name=NULL) +
  scale_x_discrete(labels=scales::label_parse()) +
  labs(x=NULL, y="trait sensitivity (°C, or Δlog slope) — scalars: per +1 SD; profile: full real–uniform swap") +
  theme_article(12) + theme(legend.position="bottom")
OUTN <- paste0("out_files/Chapter1/figures/Fig3_attribution_native", if (PERIOD=="hot") "_hot" else "")
ggsave_article(OUTN, p, 9.5, 6.2)
cat("saved", OUTN, "\n")
