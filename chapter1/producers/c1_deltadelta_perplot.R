# ==============================================================================
# Fig 2 & 3 — ΔMÉTRIQUE / ΔVARIABLE (delta/delta), par archétype ET par placette.
# The "importance" figure the supervisors asked for: the actual per-unit slope of
# each metric wrt each trait (ΔTmax/ΔLAI, Δslope/ΔfCover, …), NOT a slope×spread
# bar. Metrics = ΔTmax and the micro–macro buffering slope (all-days); variables =
# LAI, fCover, Hmax (per native step) + LAD (real−uniform contrast).
#   Fig2 = median per archetype (points+line) ; Fig3 = per-plot distribution (box).
# Reads metrics6_<VER> (per-plot, per-metric per-unit effects; NO sims).
#   VER=v320 Rscript c1_deltadelta_perplot.R     (rerun: VER=iter)
# Out: FigSh_deltadelta_perarchetype.png (Fig2) + FigSh_deltadelta_perplot.png (Fig3)
#      + tab_deltadelta_perplot.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R"); source("scripts/_article_style.R")
})
VER <- Sys.getenv("VER", "v320")
PERIOD <- Sys.getenv("PERIOD", "all")   # all | hot  (hot = 10% hottest days, annex)
MET <- if (PERIOD == "hot") c(Tmax_hot = "Delta*T[max]~(hot)", slope_hot = "Delta*log(slope)~(hot)") else
                            c(Tmax_all = "Delta*T[max]",       slope_all = "Delta*log(slope)")
DIR <- sprintf("out_files/Chapter1/tables/metrics6_%s", VER)
M <- rbindlist(lapply(list.files(DIR, "part_.*csv$", full.names = TRUE), fread), fill = TRUE)
stopifnot(nrow(M) > 0)
M <- M[metric %in% names(MET)]
M[, P := relabel_cluster(Cluster)]
# corrected per-plot VCI (pid S%04d = cLHS sample row order; sample VCI is normalized-height)
.clhs <- if (VER == "native20") "out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds" else "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"
vci <- as.data.table(readRDS(.clhs))[, .(pid = sprintf("S%04d", .I), VCI)]
M <- merge(M, vci, by = "pid", all.x = TRUE)

# fCover boundaries (as in c1_sensitivity_perplot_figure): saturated add->0, floor excluded from remove
M[, fCov_add := ifelse(fCover >= 0.95, 0, fCov_add)]
M[, fCov_rem := ifelse(fCover <= 0.5001, NA_real_, fCov_rem)]

# per-plot signed Δmetric/Δvar (native step) — SEPARATE slope columns, trait cols kept for SD
M[, `:=`(sl_LAI    = (LAI_add  - LAI_rem )/2,   # per +1 one-sided LAI
         sl_fCover = (fCov_add - fCov_rem)/2,   # per +0.1 fCover
         sl_Hmax   = (Hmax_add - Hmax_rem)/2,   # per +5 m
         sl_LAD    = dT_LAD)]                    # real - uniform contrast (no per-unit)
# slope metric -> express on LOG(slope) (buffering convention, §2.5b): Δlog(β) ≈ Δβ / β_base
issl <- grepl("slope", M$metric)
for (cc in c("sl_LAI","sl_fCover","sl_Hmax","sl_LAD"))
  M[issl & is.finite(base) & base > 0, (cc) := get(cc) / base]

# OPTIONAL within-cluster normalisation (NORM=sd): express each slope per +1 SD of the
# trait WITHIN its own cluster (= native-step slope × cluster-SD / native-step) → traits
# comparable at the variation each REALLY shows in that cluster (committee #14). LAD kept
# as the real-uniform contrast. Trait value columns (LAI, Hmax, fCover) are two-sided cLHS.
NORM <- Sys.getenv("NORM", "native")
if (NORM == "sd") {
  M[, `:=`(sl_LAI    = sl_LAI    * sd(LAI/2,   na.rm = TRUE),        # SD of ONE-SIDED LAI per cluster
           sl_fCover = sl_fCover * (sd(fCover, na.rm = TRUE) / 0.1), # SD(fCover) in 0.1 steps
           sl_Hmax   = sl_Hmax   * (sd(Hmax,   na.rm = TRUE) / 5),   # SD(Hmax) in 5 m steps
           # LAD: per unit of structural complexity (1-VCI), then × within-cluster SD(1-VCI)=SD(VCI)
           sl_LAD    = (sl_LAD / pmax(1 - VCI, 0.05)) * sd(VCI, na.rm = TRUE)),
    by = Cluster]                                                    # WITHIN cluster
  VARLAB <- c(sl_LAI="per +1 SD LAI", sl_fCover="per +1 SD fCover", sl_Hmax="per +1 SD Hmax", sl_LAD="per +1 SD (1-VCI)")
} else {
  VARLAB <- c(sl_LAI="per +1 LAI", sl_fCover="per +0.1 fCover", sl_Hmax="per +5 m Hmax", sl_LAD="real - uniform LAD")
}
SUF <- paste0(if (NORM == "sd") "_norm" else "", if (PERIOD == "hot") "_hot" else "")
L <- melt(M, id.vars = c("pid","P","metric"), measure.vars = names(VARLAB),
          variable.name = "variable", value.name = "dd")
L <- L[is.finite(dd)]
L[, metric := factor(metric, levels = names(MET), labels = MET[names(MET)])]
L[, variable := factor(variable, levels = names(VARLAB), labels = VARLAB[names(VARLAB)])]

# table: median per archetype × metric × variable
tab <- L[, .(median = round(median(dd),4), q25 = round(quantile(dd,.25),4),
             q75 = round(quantile(dd,.75),4), n = .N), by = .(metric, variable, P)][order(metric, variable, P)]
fwrite(tab, sprintf("out_files/Chapter1/tables/tab_deltadelta_perplot%s.csv", SUF))
cat(sprintf("=== Δmetric/Δvar (VER=%s) — médianes par archétype ===\n", VER)); print(tab)

# ---- Fig 3 : per-plot distribution (box) by archetype ------------------------
p3 <- ggplot(L, aes(P, dd, fill = P)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_boxplot(width = 0.7, outlier.size = 0.4, outlier.alpha = 0.25, linewidth = 0.3) +
  facet_grid(metric ~ variable, scales = "free_y", labeller = labeller(metric = label_parsed, variable = label_value), switch = "y") +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("Per-plot Δmetric/Δtrait at each plot's own canopy (400 cLHS plots; VER=%s). Rows = metric, cols = trait. Box = archetype distribution.", VER)) +
  theme_article(11) + theme(strip.placement = "outside")
ggsave_article(sprintf("out_files/Chapter1/figures/FigSh_deltadelta_perplot%s", SUF), p3, 9.5, 5.2)

# ---- Fig 2 : median per archetype (points + line) ----------------------------
med <- L[, .(dd = median(dd)), by = .(metric, variable, P)]
p2 <- ggplot(med, aes(P, dd, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_line(colour = "grey60", linewidth = 0.5) +
  geom_point(aes(colour = P), size = 2.8) +
  facet_grid(metric ~ variable, scales = "free_y", labeller = labeller(metric = label_parsed, variable = label_value), switch = "y") +
  scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("Δmetric/Δtrait per archetype (median over the archetype's cLHS plots; VER=%s). Rows = metric, cols = trait.", VER)) +
  theme_article(11) + theme(strip.placement = "outside", legend.position = "bottom")
ggsave_article(sprintf("out_files/Chapter1/figures/FigSh_deltadelta_perarchetype%s", SUF), p2, 9.5, 5.2)
cat("DONE -> FigSh_deltadelta_perarchetype (Fig2) + FigSh_deltadelta_perplot (Fig3) + tab\n")
