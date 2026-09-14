# ==============================================================================
# What does LiDAR add to MuSICA, vs no LiDAR? (supervisor: "en quoi le lidar
# apporte des infos pertinentes pour MuSICA vs sans lidar ?")
#
# WITHOUT LiDAR, MuSICA has no per-plot structure: it must assume ONE generic
# canopy for the whole forest = the GLOBAL baseline coalition (0000_global =
# global-mean LAI/Hmax/fCover + uniform LAD) -> a single predicted ΔTmax for
# every plot (dTmax_base, constant).
# WITH LiDAR, each plot gets its real LAI/Hmax/fCover/LAD -> dTmax_full (per plot).
# Value of LiDAR (per plot) = dTmax_full - dTmax_base : the microclimate signal a
# structure-blind average canopy MISSES. Composition = the global Shapley φ
# (which LiDAR variable carries that information).
# Reads shapley_parts_global (LEGACY, after the binary standardization). No sims.
#   Rscript c1_lidar_value.R
# Out: out_files/Chapter1/figures/FigSh_lidar_value.png + tab_lidar_value.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); source("R/cluster_relabel.R")
})
SG <- rbindlist(lapply(list.files("out_files/Chapter1/tables/shapley_parts_global",
                                  "part_.*csv$", full.names = TRUE), fread), fill = TRUE)
stopifnot(nrow(SG) > 0)
SG[, P := relabel_cluster(Cluster)]
SG[, value := dTmax_full - dTmax_base]          # °C the LiDAR structure moves each plot vs the average canopy
base_val <- mean(SG$dTmax_base, na.rm = TRUE)   # the single "no-LiDAR" average-canopy prediction

Fv <- c("LAI","Hmax","fCover","LAD")
cat(sprintf("=== Value of LiDAR for MuSICA ===\n"))
cat(sprintf("'No-LiDAR' average-canopy ΔTmax (single value, all plots): %.2f °C\n", base_val))
cat(sprintf("'With-LiDAR' per-plot ΔTmax: range %.2f .. %.2f °C (span %.2f), SD %.2f\n",
            min(SG$dTmax_full,na.rm=TRUE), max(SG$dTmax_full,na.rm=TRUE),
            diff(range(SG$dTmax_full,na.rm=TRUE)), sd(SG$dTmax_full,na.rm=TRUE)))
summ <- SG[, .(n=.N, base=round(mean(dTmax_base),2), full_mean=round(mean(dTmax_full),2),
               value_mean=round(mean(value),2), value_absmean=round(mean(abs(value)),2),
               value_sd=round(sd(value),2)), by=P][order(P)]
print(summ)
phi <- SG[, lapply(.SD, function(x) mean(abs(x),na.rm=TRUE)), .SDcols=Fv]
cat("global mean|phi| (info carried by each LiDAR variable):\n"); print(round(phi,3))
fwrite(SG[, .(pid, P, dTmax_base, dTmax_full, value)], "out_files/Chapter1/tables/tab_lidar_value.csv")

# (a) per-plot WITH-LiDAR distribution vs the single NO-LiDAR average canopy
pa <- ggplot(SG, aes(P, dTmax_full, fill=P, colour=P)) +
  geom_hline(yintercept=base_val, linetype="dashed", colour="grey30", linewidth=0.7) +
  geom_boxplot(width=0.6, alpha=0.35, outlier.shape=NA, linewidth=0.6) +
  geom_jitter(width=0.12, size=0.7, alpha=0.35) +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") + scale_colour_manual(values=PAL_CLUSTER, guide="none") +
  labs(title="(a) LiDAR resolves a per-plot spread (dashed = no-LiDAR avg canopy)", x=NULL, y=expression("With-LiDAR "*Delta*T[max]~("°C"))) +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        plot.title=element_text(size=11,face="bold"))

dphi <- data.table(variable=factor(Fv, levels=Fv), phi=as.numeric(phi[1]))
pb <- ggplot(dphi, aes(variable, phi, fill=variable)) +
  geom_col(width=0.66, colour="grey25", linewidth=0.2) +
  scale_fill_manual(values=c(LAI="#2C7BB6",fCover="#5E3C99",Hmax="#E66101",LAD="#1A9850"), guide="none") +
  labs(title="(b) Which LiDAR variable carries the info", x=NULL, y=expression("global mean"~"|"*phi*"|"~("°C"))) +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        plot.title=element_text(size=11,face="bold"))

p <- pa + pb + plot_layout(widths=c(1.4,1)) +
  plot_annotation(caption=sprintf("Without LiDAR, MuSICA assumes one average canopy (%.1f °C for all plots). LiDAR structure moves each plot, resolving a %.1f °C span carried mostly by LAI.",
                                  base_val, diff(range(SG$dTmax_full,na.rm=TRUE)))) &
  theme(plot.caption=element_text(size=9, colour="grey40", hjust=0))
ggsave("out_files/Chapter1/figures/FigSh_lidar_value.png", p, width=11, height=4.8, dpi=200, bg="white")
cat("DONE -> FigSh_lidar_value.png + tab_lidar_value.csv\n")
