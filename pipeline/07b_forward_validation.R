# ==============================================================================
# PIPELINE STAGE 07b — Forward validation in SHAPLEY-RANK order, per cluster + All.
# Adds traits cumulatively in the order of the pooled |median φ| ranking (05b/06b)
# and validates each step's simulated mean ΔTmax against the 53 HOBO at 1 m.
# Rows = P1..P4 + All ; cols = Baseline → +trait → … → REF. r / RMSE per panel.
# Output: fig_validation_forward_shapley.{png,pdf} + tab_validation_forward.csv
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
cli_h1("STAGE 07b — forward validation (Shapley order)")

coal <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))
clu  <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))
V    <- readRDS(file.path(PIPE$OUT_DATA, "ref_validation.rds"))
df_macro <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq))
mean_macro <- mean(df_macro$Tmax_macro, na.rm = TRUE)

# obs per-plot mean ΔTmax
obs_pp <- V$obs_daily[, .(Delta_obs = mean(Delta_obs, na.rm=TRUE)), by=id_plot]

# Shapley-rank order (pooled, ΔTmax, all period) from 06b ranking table
rk <- fread(file.path(PIPE$OUT_TAB, "tab_shapley_ranking.csv"))
ord <- rk[metric=="Tmax" & period=="All period" & Cluster=="All"][order(rank), trait]
cli_alert("Shapley order: {paste(ord, collapse=' > ')}")

pos <- c(LAI=1L, Hmax=2L, fCover=3L, LAD=4L)               # bit positions
mk  <- function(active){ b<-rep("0",4L); b[active]<-"1"; paste(b,collapse="") }
steps <- list(list(bit="0000", lab="Baseline")); active <- integer(0)
for (k in seq_along(ord)) {
  active <- c(active, pos[[ord[k]]])
  lab <- if (k==length(ord)) sprintf("+%s (REF)", ord[k]) else sprintf("+%s", ord[k])
  steps[[k+1L]] <- list(bit=mk(active), lab=lab)
}
step_lab <- vapply(steps, `[[`, "", "lab")

# build per-plot sim ΔTmax at each step
long <- rbindlist(lapply(seq_along(steps), function(i){
  b <- steps[[i]]$bit
  d <- coal[bit==b, .(id_plot, Delta_sim = Tmax_all - mean_macro)]
  d <- merge(merge(d, obs_pp, by="id_plot"), clu, by="id_plot")
  d[, `:=`(step=i, step_lab=factor(steps[[i]]$lab, levels=step_lab))]; d
}))
# add pooled "All" cluster
longA <- copy(long)[, Cluster := "All"]
long <- rbind(long, longA)
long[, Cluster := factor(Cluster, levels=c(paste0("P",1:4),"All"))]

# per (cluster, step) fit stats
ann <- long[, { ok <- is.finite(Delta_sim)&is.finite(Delta_obs)
  .(r=cor(Delta_sim[ok],Delta_obs[ok]),
    RMSE=sqrt(mean((Delta_sim[ok]-Delta_obs[ok])^2)), n=sum(ok)) },
  by=.(Cluster, step, step_lab)]
ann[, lab := sprintf("italic(r)==%.2f", r)]
ann[, lab2:= sprintf("RMSE==%.2f", RMSE)]
fwrite(ann, file.path(PIPE$OUT_TAB, "tab_validation_forward.csv"))
cli_h2("Forward validation r per cluster × step"); print(dcast(ann, Cluster~step_lab, value.var="r")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)])

pal <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850", All="grey30")
lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm=TRUE); lim <- lim + c(-1,1)*diff(lim)*0.06
p <- ggplot(long, aes(Delta_obs, Delta_sim, colour=Cluster)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55") +
  geom_hline(yintercept=0, linetype="dotted", colour="grey85") +
  geom_vline(xintercept=0, linetype="dotted", colour="grey85") +
  geom_point(size=1.5, alpha=0.8) +
  geom_text(data=ann[is.finite(r)], aes(x=-Inf,y=Inf,label=lab), parse=TRUE, hjust=-0.12, vjust=1.5, size=3.6, inherit.aes=FALSE) +
  geom_text(data=ann, aes(x=-Inf,y=Inf,label=lab2), parse=TRUE, hjust=-0.10, vjust=3.0, size=3.4, inherit.aes=FALSE) +
  facet_grid(Cluster ~ step_lab) +
  coord_equal(xlim=lim, ylim=lim) + scale_colour_manual(values=pal, guide="none") +
  labs(x=expression("Observed mean "*Delta*T[max]~"("*degree*"C)"),
       y=expression("Simulated mean "*Delta*T[max]~"("*degree*"C)"),
       title="Forward validation in Shapley-rank order — per cluster + All",
       caption="Traits added cumulatively by pooled |Shapley φ| rank. 1:1 dashed. r/RMSE vs 53 HOBO at 1 m.") +
  theme_bw(base_size=15) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=13), plot.caption=element_text(size=10.5,colour="grey40"))
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_forward_shapley.png"), p, width=17, height=18, dpi=200, bg="white")
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_forward_shapley.pdf"), p, width=17, height=18, device=cairo_pdf)
cli_alert_success("Saved fig_validation_forward_shapley")
