# ==============================================================================
# Annex — LOCAL CURVATURE of the trait sensitivities (saturation probe), FREE.
# No new simulations: uses the per-plot _add / _rem effects already in metrics6.
#
# For each trait we already simulated BOTH sides of the plot's own canopy:
#   upper-side per-unit slope = _add        ( [f(x+h)-f(x)] / +h )
#   lower-side per-unit slope = -_rem        ( [f(x)-f(x-h)] / +h )
# Under LOCAL LINEARITY the two coincide (add = -rem). The GAP between them is
# the local curvature (∝ the 2nd finite difference add+rem) and, for a buffering
# trait, a flatter UPPER side (|upper| < |lower|) = SATURATION / diminishing
# returns as the trait grows. This is the within-cluster, per-trait version of
# the density-dependent saturation story — at zero extra sim cost.
#
# Scaling, log(slope) and boundary clipping are IDENTICAL to c1_deltadelta_perplot.R
# so the per-+1-SD axis matches Fig2/3.
#   VER=v320 PERIOD=all NORM=sd Rscript c1_sensitivity_asymmetry.R   (rerun: VER=iter)
# Out: FigSh_sensitivity_asymmetry[_norm][_hot].png  + tab_sensitivity_asymmetry*.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R"); source("scripts/_article_style.R")
})
VER    <- Sys.getenv("VER", "v320")
PERIOD <- Sys.getenv("PERIOD", "all")        # all | hot
NORM   <- Sys.getenv("NORM", "sd")           # sd (per +1 within-cluster SD) | native
MET <- if (PERIOD == "hot") c(Tmax_hot = "Delta*T[max]~(hot)", slope_hot = "Delta*log(slope)~(hot)") else
                            c(Tmax_all = "Delta*T[max]",       slope_all = "Delta*log(slope)")

DIR <- sprintf("out_files/Chapter1/tables/metrics6_%s", VER)
M <- rbindlist(lapply(list.files(DIR, "part_.*csv$", full.names = TRUE), fread), fill = TRUE)
stopifnot(nrow(M) > 0)
M <- M[metric %in% names(MET)]
M[, P := relabel_cluster(Cluster)]

# same boundary handling as Fig2/3: saturated fCover can't go up; floored can't go down
M[, fCov_add := ifelse(fCover >= 0.95,   0,        fCov_add)]
M[, fCov_rem := ifelse(fCover <= 0.5001, NA_real_, fCov_rem)]

# per-unit UPPER-side and LOWER-side slopes (native step), both signed "Δmetric per +1 unit"
M[, `:=`(up_LAI = LAI_add,  dn_LAI = -LAI_rem,
         up_fCover = fCov_add, dn_fCover = -fCov_rem,
         up_Hmax = Hmax_add, dn_Hmax = -Hmax_rem)]
# slope metric -> log(slope) (buffering convention): divide by base
issl <- grepl("slope", M$metric)
for (cc in c("up_LAI","dn_LAI","up_fCover","dn_fCover","up_Hmax","dn_Hmax"))
  M[issl & is.finite(base) & base > 0, (cc) := get(cc) / base]

# within-cluster per-+1-SD scaling (identical factors to c1_deltadelta_perplot.R NORM=sd)
if (NORM == "sd") {
  M[, `:=`(sLAI = sd(LAI/2, na.rm = TRUE),
           sFC  = sd(fCover, na.rm = TRUE) / 0.1,
           sHM  = sd(Hmax,   na.rm = TRUE) / 5), by = Cluster]
  M[, `:=`(up_LAI = up_LAI*sLAI, dn_LAI = dn_LAI*sLAI,
           up_fCover = up_fCover*sFC, dn_fCover = dn_fCover*sFC,
           up_Hmax = up_Hmax*sHM, dn_Hmax = dn_Hmax*sHM)]
  UNIT <- "per +1 SD"
} else UNIT <- "native step"
SUF <- paste0(if (NORM == "sd") "_norm" else "", if (PERIOD == "hot") "_hot" else "")

# long: trait x side
TRAITS <- c(LAI = "LAI", fCover = "fCover", Hmax = "H[max]")
L <- rbindlist(lapply(names(TRAITS), function(tr) rbindlist(list(
  M[, .(pid, P, metric, trait = tr, side = "upper (+1 SD)", val = get(paste0("up_", tr)))],
  M[, .(pid, P, metric, trait = tr, side = "lower (-1 SD)", val = get(paste0("dn_", tr)))]))))
L <- L[is.finite(val)]
L[, metric := factor(metric, levels = names(MET), labels = MET[names(MET)])]
L[, trait  := factor(trait,  levels = names(TRAITS), labels = TRAITS[names(TRAITS)])]
L[, side   := factor(side,   levels = c("lower (-1 SD)", "upper (+1 SD)"))]

# ---- curvature table: upper vs lower median + asymmetry index (upper - lower) ----
med <- L[, .(med = median(val), q25 = quantile(val,.25), q75 = quantile(val,.75), n = .N),
          by = .(metric, trait, P, side)]
asym <- dcast(med[, .(metric, trait, P, side, med)], metric + trait + P ~ side, value.var = "med")
setnames(asym, c("lower (-1 SD)","upper (+1 SD)"), c("lower","upper"), skip_absent = TRUE)
asym[, curvature := round(upper - lower, 4)]      # >0 with cooling trait => upper flatter => saturating
asym[, c("upper","lower") := .(round(upper,4), round(lower,4))]
fwrite(asym, sprintf("out_files/Chapter1/tables/tab_sensitivity_asymmetry%s.csv", SUF))
cat(sprintf("=== upper vs lower slope (%s ; VER=%s) ; curvature=upper-lower ===\n", UNIT, VER))
print(asym[order(metric, trait, P)])

# ---- figure: paired upper/lower median (+ IQR) per cluster, faceted by trait x metric ----
dodge <- position_dodge(width = 0.55)
p <- ggplot(med, aes(P, med, colour = side, group = side)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_linerange(aes(ymin = q25, ymax = q75), position = dodge, linewidth = 0.5, alpha = 0.6) +
  geom_point(position = dodge, size = 2.6) +
  facet_grid(metric ~ trait, scales = "free_y",
             labeller = labeller(metric = label_parsed, trait = label_parsed), switch = "y") +
  scale_colour_manual(values = c("lower (-1 SD)" = "#2C7FB8", "upper (+1 SD)" = "#E66101"), name = NULL) +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("Upper- vs lower-side sensitivity (%s; VER=%s). Coinciding = locally linear; gap = curvature (upper flatter = saturation).",
                          UNIT, VER)) +
  theme_article(11) + theme(strip.placement = "outside", legend.position = "bottom")
ggsave_article(sprintf("out_files/Chapter1/figures/FigSh_sensitivity_asymmetry%s", SUF), p, 9.5, 5.2)
cat("DONE -> FigSh_sensitivity_asymmetry", SUF, " + tab_sensitivity_asymmetry", SUF, "\n", sep="")
