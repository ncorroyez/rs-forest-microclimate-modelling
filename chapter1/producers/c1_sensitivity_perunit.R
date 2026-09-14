# ==============================================================================
# PER-UNIT marginal sensitivity (alternative to Shapley). Supervisor's intent for
# "ΔTmax/ΔLAI, ΔTmax/ΔHmax, ΔTmax/ΔfCover": not the range-normalized importance
# (Fig 15) and not Shapley (whose averaging of contributions may not be the right
# readout) — but the CONCRETE physical effect per ONE NATIVE UNIT of each variable,
# at each cluster's operating point:
#   °C per 1 (one-sided) LAI unit  |  °C per 1 m Hmax  |  °C per 0.1 fCover (10 pp)
# Local slope = Δmetric / Δvariable across the within-cluster p10-p90 window
# (reads the cached endpoints in tab_sensitivity_points.csv — NO new sims).
# LAI is stored 2x (two-sided); we report per ONE-SIDED unit (display convention),
# i.e. slope_per_stored x 2.
#   Rscript c1_sensitivity_perunit.R
# Out: out_files/Chapter1/figures/FigSh_sensitivity_perunit.png + tab_sensitivity_perunit.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R")
})
P <- fread("out_files/Chapter1/tables/tab_sensitivity_points.csv")
MET <- c("Tmax_all","VPD_all","slope_all")

perunit <- rbindlist(lapply(sort(unique(P$Cluster)), function(cid) {
  g <- function(t) P[Cluster == cid & tag == t]
  slope <- function(hi, lo, xcol, fac) {
    a <- g(hi); b <- g(lo); dx <- (a[[xcol]] - b[[xcol]])
    setNames(as.numeric(unlist(a[, ..MET]) - unlist(b[, ..MET])) / dx * fac, MET)
  }
  rbind(
    data.table(Cluster=cid, variable="LAI",    unit="per 1 LAI (one-sided)", as.data.table(as.list(slope("LAIloc_hi","LAIloc_lo","lai",  2)))),   # x2 -> one-sided
    data.table(Cluster=cid, variable="Hmax",   unit="per 1 m",               as.data.table(as.list(slope("Hmaxloc_hi","Hmaxloc_lo","hmax", 1)))),
    data.table(Cluster=cid, variable="fCover", unit="per 0.1 (10 pp)",       as.data.table(as.list(slope("fCovloc_hi","fCovloc_lo","fcover",0.1))))
  )
}))
perunit[, P := relabel_cluster(Cluster)]
fwrite(perunit, "out_files/Chapter1/tables/tab_sensitivity_perunit.csv")
cat("=== ΔTmax per native unit, per cluster (°C) ===\n")
print(dcast(perunit, P ~ variable, value.var="Tmax_all")[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])
cat("\n(VPD_all in kPa/unit, slope_all dimensionless/unit also in the CSV)\n")

LABS <- c(LAI="Delta*T[max]~per~1~LAI~unit~(degree*C)",
          Hmax="Delta*T[max]~per~1~m~H[max]~(degree*C)",
          fCover="Delta*T[max]~per~0.1~fCover~(degree*C)")
d <- copy(perunit); d[, variable := factor(variable, levels=names(LABS), labels=LABS)]
p <- ggplot(d, aes(P, Tmax_all, fill=P)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=0.4) +
  geom_col(width=0.7, colour="grey25", linewidth=0.2) +
  geom_text(aes(label=sprintf("%+.3g", Tmax_all), vjust=ifelse(Tmax_all>=0,-0.4,1.3)), size=3.2) +
  facet_wrap(~variable, scales="free_y", nrow=1, labeller=label_parsed) +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") +
  labs(x=NULL, y=NULL, title="Marginal ΔTmax sensitivity per native unit (local, at each cluster operating point)") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11, face="bold"))
ggsave("out_files/Chapter1/figures/FigSh_sensitivity_perunit.png", p, width=11, height=4.4, dpi=200, bg="white")
cat("DONE -> FigSh_sensitivity_perunit.png + tab_sensitivity_perunit.csv\n")
