# ==============================================================================
# Figure: per-plot marginal ΔTmax sensitivity DISTRIBUTIONS over the 400 cLHS
# plots (the distribution version of fig 22/23, which gave 4 cluster points).
# One marginal slope per plot per variable per direction; shown as per-cluster
# boxplots. Reads sensitivity_perplot/part_*.csv (NO sims).
#   Rscript c1_sensitivity_perplot_figure.R
# Out: FigSh_sensitivity_perplot.png + tab_sensitivity_perplot_summary.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R")
})
S <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot",
                                 "part_.*csv$", full.names = TRUE), fread), fill = TRUE)
stopifnot(nrow(S) > 0)
S[, P := relabel_cluster(Cluster)]
cat(sprintf("plots: %d (per cluster: %s)\n", nrow(S), paste(table(S$P), collapse="/")))

# ---- fCover boundary handling (mix of options 1 + 2) ------------------------
# UPPER saturation: fCover >= 0.95 cannot meaningfully gain cover -> the add-effect
# is 0 by saturation (option 2). This also removes the divide-by-tiny-step blow-up
# for near-saturated plots. LOWER bound: fCover at the cLHS floor (0.5) is an
# artificial sampling cap, NOT physical saturation, so removing below it is
# meaningless -> those plots are EXCLUDED from remove(-) (option 1), not coded 0.
SAT_HI <- 0.95; FLOOR_LO <- 0.5001
S[, fCov_add := ifelse(fCover >= SAT_HI, 0, fCov_add)]          # saturated add -> 0
S[, fCov_rem := ifelse(fCover <= FLOOR_LO, NA_real_, fCov_rem)] # floored remove -> excluded
cat(sprintf("fCover: %d plots saturated (>=%.2f, add coded 0); %d at floor (<=%.2f, remove excluded)\n",
            sum(S$fCover >= SAT_HI), SAT_HI, sum(S$fCover <= FLOOR_LO), 0.5))

L <- melt(S, id.vars = c("pid","P"),
          measure.vars = c("LAI_add","LAI_rem","Hmax_add","Hmax_rem","fCov_add","fCov_rem"),
          variable.name = "key", value.var = "dT")
L[, variable := sub("_(add|rem)$", "", key)]
L[, sign := ifelse(grepl("_add$", key), "add (+)", "remove (-)")]
L <- L[is.finite(value)]

# LAD: categorical contrast real vs uniform (not a +/- unit) — 4th panel, single box
LAD <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")
LAD[, P := relabel_cluster(Cluster)]
Llad <- LAD[is.finite(dT_LAD), .(pid, P, key="LAD_realunif", variable="LAD",
                                 sign="real vs uniform", value=dT_LAD)]
L <- rbind(L, Llad, fill=TRUE)

LABS <- c(LAI ="Delta*T[max]~per~1~LAI~unit~(degree*C)",
          fCov="Delta*T[max]~per~0.1~fCover~(degree*C)",
          Hmax="Delta*T[max]~per~5~m~H[max]~(degree*C)",
          LAD ="Delta*T[max]~real-uniform~LAD~(degree*C)")
L[, variable := factor(variable, levels=names(LABS), labels=LABS)]

# summary table (median per cluster x variable x direction)
summ <- L[, .(median=round(median(value),3), q25=round(quantile(value,.25),3),
              q75=round(quantile(value,.75),3), n=.N), by=.(P, variable=sub("~.*","",as.character(variable)), sign)]
fwrite(dcast(L[, .(median=median(value)), by=.(P, var=sub("_(add|rem)$","",key), sign)],
             P + var ~ sign, value.var="median"),
       "out_files/Chapter1/tables/tab_sensitivity_perplot_summary.csv")
cat("=== median ΔTmax per native unit, per cluster (°C) ===\n")
print(dcast(L[, .(median=round(median(value),3)), by=.(P, var=sub("_(add|rem)$","",key), sign)],
            P + var ~ sign, value.var="median"))

# N per box (after boundary handling), anchored at each facet's lower bound
NC <- L[, .(n=.N), by=.(variable, P, sign)]
NC <- merge(NC, L[, .(ytext = min(value, na.rm=TRUE)), by=variable], by="variable")

p <- ggplot(L, aes(P, value, fill=sign)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=0.4) +
  geom_boxplot(width=0.7, outlier.size=0.5, outlier.alpha=0.3, linewidth=0.3,
               position=position_dodge(width=0.8)) +
  geom_text(data=NC, aes(y=ytext, label=n, group=sign), vjust=1.6, size=2.5,
            colour="grey30", position=position_dodge(width=0.8)) +
  scale_y_continuous(expand=expansion(mult=c(0.10, 0.05))) +
  facet_wrap(~variable, scales="free_y", nrow=1, labeller=label_parsed) +
  scale_fill_manual(values=c("add (+)"="#2C7BB6","remove (-)"="#D7191C","real vs uniform"="#1A9850"), name=NULL) +
  labs(x=NULL, y=NULL,
       title="Per-plot marginal ΔTmax sensitivity (400 cLHS plots, distribution per archetype)",
       subtitle="Box = per-plot effect at each plot's own operating point (other traits fixed at real); N per box. fCover: saturated plots (≥0.95) coded 0 for add, floored (0.5) excluded from remove. LAD = real-vs-uniform profile contrast.") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
        legend.position="bottom", plot.title=element_text(size=11,face="bold"),
        plot.subtitle=element_text(size=8, colour="grey35"))
ggsave("out_files/Chapter1/figures/FigSh_sensitivity_perplot.png", p, width=14, height=5.0, dpi=200, bg="white")
cat("DONE -> FigSh_sensitivity_perplot.png + tab_sensitivity_perplot_summary.csv\n")
