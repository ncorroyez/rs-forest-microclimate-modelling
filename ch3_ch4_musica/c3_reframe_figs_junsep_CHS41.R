# ==============================================================================
# c3_reframe_figs_junsep_CHS41.R — Chapter 3 figures, June-September window, on the
# CHS41-Rmerge / no-wind extraction (macro = the CHS 41 station). Same code as
# c3_reframe_figs_junsep.R; only the source tables and the output names differ.
# Built fresh from the junsep_full_* tables (the reframe pipeline scripts are gone).
#   Fig 1 premise    : coupling slope vs total LAI (obs + LiDAR-forced sim)
#   Fig 2 sensitivity: paired sim contrast (scenario - LiDAR fixe) by archetype, slope + ΔTmax
#   Fig 3 validation : ranking r (sim vs obs) by scenario x archetype, slope + ΔTmax
#   Fig 4 mechanism  : S2 slope error vs fraction of leaf area in top d_opt (P3+P4)
# Out: NC_Full/manuscripts/ch3/figures/Fig{1,2,3,4}_reframe_junsep_chs41.png
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
source("R/cluster_relabel.R")
TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
FIG <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures"; dir.create(FIG,showWarnings=FALSE,recursive=TRUE)
Zdir  <- "/home/corroyez/Documents/z_Example_rmusica_31012025/out_files/Chapter3_CHS41/tables"
Zdir0 <- "/home/corroyez/Documents/z_Example_rmusica_31012025/out_files/Chapter3/tables"  # LAI_ALS lives here
th <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(), legend.position="top")
SCEN <- c("LiDARfixe","S2seul","S2opt","Combinaison")
SCEN_EN <- c(LiDARfixe="LiDAR fixed", S2seul="Sentinel-2 alone",
             S2opt="Sentinel-2 \"opt\"", Combinaison="Combination")
scol <- c(LiDARfixe="#1A9850", S2seul="#D7191C", S2opt="#FDAE61", Combinaison="#2C7BB6")

PP  <- fread(file.path(Zdir,"junsep_chs41_perplot_macro-chs41.csv"))
# LAI_ALS is the leaf area MuSICA was forced with. The clustering LAI in
# reframe_seasonality_amplitude_RAW.csv differs on the open plots (up to x3) and
# must not be used here: the chapter reports the forcing LAI.
LAI <- unique(fread(file.path(Zdir0,"annual_notmasked_report.csv"))[,.(id_plot,LAI=LAI_ALS)])
FRT <- fread(file.path(TAB,"reframe_dopt_mechanism.csv"))[,.(id_plot,fr_top)]
BY  <- fread(file.path(Zdir,"junsep_chs41_by_archetype_macro-chs41.csv"))[P!="ALL"]
CON <- fread(file.path(Zdir,"junsep_chs41_paired_contrasts_macro-chs41.csv"))
W   <- merge(dcast(PP,id_plot+P+obs_slope_junsep~scenario,value.var="sim_slope_junsep"), LAI, by="id_plot")
W[,P:=factor(P,levels=paste0("P",1:4))]

# ---- Fig 1 : premise (slope and dTmax vs forcing LAI) ----
W2 <- merge(dcast(PP, id_plot + P + obs_dtmax_junsep ~ scenario,
                  value.var = "sim_dtmax_junsep"), LAI, by = "id_plot")
W2[, P := factor(P, levels = paste0("P", 1:4))]
premise_panel <- function(D, obs, sim, ylab) {
  d <- rbind(D[, .(P, LAI, y = get(obs), src = "Observed (HOBO)")],
             D[, .(P, LAI, y = get(sim), src = "Simulated (LiDAR)")])
  r_o <- cor(D[[obs]], D$LAI, use = "complete")
  r_s <- cor(D[[sim]], D$LAI, use = "complete")
  ggplot(d, aes(LAI, y, colour = P, shape = src)) +
    geom_point(size = 2, alpha = .85) +
    geom_smooth(aes(group = src, linetype = src), method = "lm", se = FALSE,
                colour = "grey30", linewidth = .6) +
    scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
    scale_shape_manual(values = c(16, 2), name = NULL) +
    scale_linetype_manual(values = c(1, 2), guide = "none") +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 3.2,
             label = sprintf("observed r = %.2f\nsimulated r = %.2f", r_o, r_s)) +
    labs(x = expression("LiDAR LAI (one-sided, m"^2*" m"^-2*")"), y = ylab) + th
}
g1 <- premise_panel(W,  "obs_slope_junsep", "LiDARfixe", "Coupling slope (T_micro ~ T_macro)") +
      premise_panel(W2, "obs_dtmax_junsep", "LiDARfixe", expression(Delta*T[max]*" ("*degree*"C)")) +
      plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave(file.path(FIG,"Fig1_reframe_junsep_chs41.png"), g1, width = 9.5, height = 4.6, dpi = 300, bg = "white")

# ---- Fig 2 : sensitivity (paired contrast vs LiDAR fixe) ----
C2 <- CON[grepl("S2seul|S2opt|Combinaison",contrast)]
C2[,scen:=sub("-LiDARfixe","",contrast)]; C2[,P:=factor(P,levels=paste0("P",1:4))]
C2[,scen:=factor(SCEN_EN[scen],levels=unname(SCEN_EN[c("S2seul","S2opt","Combinaison")]))]
p2 <- function(met,ylab){ ggplot(C2[metric==met],aes(P,mean,colour=scen,group=scen)) +
  geom_hline(yintercept=0,colour="grey45",linewidth=.4) +
  geom_linerange(aes(ymin=lo,ymax=hi),position=position_dodge(.5),linewidth=.6) +
  geom_point(position=position_dodge(.5),size=2.2) +
  scale_colour_manual(values=setNames(unname(scol[c("S2seul","S2opt","Combinaison")]),
                    unname(SCEN_EN[c("S2seul","S2opt","Combinaison")])),name=NULL) +
  labs(x=NULL,y=ylab) + th }
g2 <- p2("slope","Δ coupling slope vs LiDAR") + p2("dtmax","Δ ΔT_max vs LiDAR (°C)") +
  plot_layout(guides="collect") & theme(legend.position="top")
ggsave(file.path(FIG,"Fig2_reframe_junsep_chs41.png"),g2,width=9,height=4.4,dpi=300,bg="white")

# ---- Fig 3 : validation (ranking r) ----
B <- BY[,.(scenario,P,metric,r)]; B[,P:=factor(P,levels=paste0("P",1:4))]
B[,scenario:=factor(SCEN_EN[scenario],levels=unname(SCEN_EN[SCEN]))]
p3 <- function(met,ylab){ ggplot(B[metric==met],aes(P,r,fill=scenario)) +
  geom_hline(yintercept=0,colour="grey45",linewidth=.4) +
  geom_col(position=position_dodge(.8),width=.7) +
  scale_fill_manual(values=setNames(unname(scol[SCEN]),unname(SCEN_EN[SCEN])),name=NULL) + coord_cartesian(ylim=c(-0.4,0.9)) +
  labs(x=NULL,y=ylab) + th }
g3 <- p3("slope","Ranking r, coupling slope") + p3("dtmax","Ranking r, ΔT_max") +
  plot_layout(guides="collect") & theme(legend.position="top")
ggsave(file.path(FIG,"Fig3_reframe_junsep_chs41.png"),g3,width=9.5,height=4.4,dpi=300,bg="white")

# ---- Fig 4 : mechanism (S2 slope error vs fr_top, P3+P4) ----
M <- merge(W[,.(id_plot,P,S2seul,obs_slope_junsep)],FRT,by="id_plot")
M <- M[P%in%c("P3","P4")&is.finite(fr_top)]; M[,err:=abs(S2seul-obs_slope_junsep)]
r4 <- cor(M$fr_top,M$err,use="complete")
g4 <- ggplot(M,aes(fr_top,err,colour=P)) + geom_point(size=2.4,alpha=.85) +
  geom_smooth(method="lm",se=FALSE,colour="grey30",linewidth=.6) +
  scale_colour_manual(values=PAL_CLUSTER[c("P3","P4")],name=NULL) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.5,size=3.4,label=sprintf("r = %.2f",r4)) +
  labs(x="Fraction of leaf area in top d_opt",y="|Sentinel-2 slope error|") + th
ggsave(file.path(FIG,"Fig4_reframe_junsep_chs41.png"),g4,width=6,height=4.6,dpi=300,bg="white")

cat("DONE: Fig1-4_reframe_junsep.png in",FIG,"\n")
