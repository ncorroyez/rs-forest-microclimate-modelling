# ==============================================================================
# Forward validation, PER-CLUSTER addition order + an "All" row (global order).
# Each archetype adds its LiDAR traits cumulatively in ITS OWN |φ| importance
# order (from the per-cluster cLHS Shapley, ΔTmax all period); the "All" row uses
# the order pooled over all 400 points. Each cumulative step's simulated mean
# ΔTmax (z05, 1 m) is validated against the 53 HOBO of that cluster (r, RMSE).
# Shows that LAD only helps the fit where it leads (dense P4), LAI in the open.
# Run:  PIPE_BRANCH=z05 Rscript c1_forward_percluster.R
# Out:  out_files/Chapter1/figures/FigSh_forward_percluster.png  (+ table)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
if (!exists("PIPE")) source("pipeline/00_config.R")
coal <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))
clu  <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))      # id_plot, Cluster (P1..P4)
V    <- readRDS(file.path(PIPE$OUT_DATA, "ref_validation.rds"))
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq)); mean_macro <- mean(dm$Tmax_macro, na.rm=TRUE)
obs_pp <- as.data.table(V$obs_daily)[, .(Delta_obs = mean(Delta_obs, na.rm=TRUE)), by=id_plot]

# ---- per-cluster + pooled trait order from SENSITIVITY importance (non-Shapley) ----
# Add-order = descending sensitivity-based importance per archetype (and pooled
# "All"); table from c1_importance_sensitivity.R. LAI leads every archetype; the
# vertical profile (LAD) is second in the dense end (P3, P4).
IMP <- fread("out_files/Chapter1/tables/tab_importance_sensitivity.csv")
Fv <- c("LAI","Hmax","fCover","LAD")
ord_of <- function(r) Fv[order(-unlist(r[, ..Fv]))]
orders <- lapply(setNames(paste0("P",1:4), paste0("P",1:4)), function(cl) ord_of(IMP[P==cl]))
orders[["All"]] <- ord_of(IMP[P=="All"])
cat("=== addition order per row ===\n"); for(nm in names(orders)) cat(sprintf("  %-3s : %s\n", nm, paste(orders[[nm]], collapse=" > ")))

pos <- c(LAI=1L, Hmax=2L, fCover=3L, LAD=4L); mk <- function(a){ b<-rep("0",4L); b[a]<-"1"; paste(b,collapse="") }
build_row <- function(rowname, order, ids){
  active <- integer(0); steps <- list(list(bit="0000", add="Baseline", step=0L))
  for(k in seq_along(order)){ active <- c(active, pos[[order[k]]])
    steps[[k+1L]] <- list(bit=mk(active), step=k,
                          add=if(k==length(order)) sprintf("+%s (REF)", order[k]) else sprintf("+%s", order[k])) }
  rbindlist(lapply(steps, function(s){
    d <- coal[bit==s$bit & id_plot %in% ids, .(id_plot, Delta_sim = Tmax_all - mean_macro)]
    d <- merge(d, obs_pp, by="id_plot"); d <- merge(d, clu[, .(id_plot, pt_clu=Cluster)], by="id_plot")
    d[, `:=`(rowlab=rowname, step=s$step, add=s$add)]; d }))
}
rows <- rbindlist(c(
  lapply(paste0("P",1:4), function(cl) build_row(cl, orders[[cl]], clu[Cluster==cl, id_plot])),
  list(build_row("All", orders[["All"]], clu$id_plot)) ))
rows[, rowlab := factor(rowlab, levels=c(paste0("P",1:4),"All"))]

ann <- rows[, { ok <- is.finite(Delta_sim)&is.finite(Delta_obs)
  .(r=if(sum(ok)>2) cor(Delta_sim[ok],Delta_obs[ok]) else NA_real_,
    RMSE=sqrt(mean((Delta_sim[ok]-Delta_obs[ok])^2, na.rm=TRUE)),
    MAE=mean(abs(Delta_sim[ok]-Delta_obs[ok]), na.rm=TRUE), n=sum(ok)) }, by=.(rowlab, step)]
ann <- merge(ann, unique(rows[,.(rowlab,step,add)]), by=c("rowlab","step"))
fwrite(ann[order(rowlab,step)], file.path(PIPE$OUT_TAB,"tab_forward_percluster.csv"))
cat("\n=== r per row × step ===\n"); print(dcast(ann, rowlab~step, value.var="r")[, lapply(.SD,function(x) if(is.numeric(x)) round(x,2) else x)])

pal <- c(P1="#D7191C",P2="#FDAE61",P3="#74C476",P4="#1A9850")
shp <- c(P1=16,P2=17,P3=15,P4=18)
lim <- range(c(rows$Delta_sim, rows$Delta_obs), na.rm=TRUE); lim <- lim + c(-1,1)*diff(lim)*0.06
slab <- c("Baseline","Step 1","Step 2","Step 3","Step 4")
rows[, step := factor(step, levels=0:4, labels=slab)]
ann[,  step := factor(step, levels=0:4, labels=slab)]
ann[, lab := paste0(ifelse(is.finite(r), sprintf("r = %+.2f\n", r), ""), sprintf("RMSE = %.2f °C\nn = %d", RMSE, n))]
p <- ggplot(rows, aes(Delta_obs, Delta_sim)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55") +
  geom_smooth(method="lm", formula=y~x, se=TRUE, colour="#E69F00", fill="#E69F00", alpha=0.16, linewidth=0.6) +
  geom_point(aes(colour=pt_clu, shape=pt_clu), size=1.7, alpha=0.85) +
  geom_text(data=ann, aes(x=-Inf, y=Inf, label=lab), hjust=-0.07, vjust=1.2, size=2.9, lineheight=0.95, inherit.aes=FALSE) +
  geom_text(data=ann, aes(x=Inf, y=-Inf, label=add), hjust=1.04, vjust=-0.7, size=2.7, fontface="italic", colour="grey25", inherit.aes=FALSE) +
  facet_grid(rowlab ~ step) + coord_equal(xlim=lim, ylim=lim) +
  scale_colour_manual(values=pal, name="Cluster") + scale_shape_manual(values=shp, name="Cluster") +
  labs(x=expression("Field "*Delta*T[max]~" ("*degree*"C)"), y=expression("Simulated "*Delta*T[max]~" ("*degree*"C)")) +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold"), legend.position="bottom")
ggsave("out_files/Chapter1/figures/FigSh_forward_percluster.png", p, width=13, height=14, dpi=200, bg="white")
cat("\nDONE -> out_files/Chapter1/figures/FigSh_forward_percluster.png\n")
