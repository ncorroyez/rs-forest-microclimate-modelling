# ==============================================================================
# PIPELINE STAGE 10 (ANNEX) — Conditional Shapley for LAD, bounding the
# OFF-MANIFOLD artefact. LAD & Hmax are co-defined (LAD lives on [0,Hmax]); the
# interventional 2^4 toggles produce chimeric coalitions where a real LAD profile
# is compressed onto a baseline (mean) Hmax. This decomposes φ_LAD into:
#   φ_LAD | Hmax=real     (coherent domain)   vs   φ_LAD | Hmax=baseline (chimera)
# φ_LAD_global = mean of the two. If the two agree → off-manifold effect is small
# → the attribution is robust. Metric = ΔTmax (all summer, 1 m). Reuses .LAD_PAIRS
# from R/h1_shapley_conditional.R. Sign: φ<0 = LAD buffers (lowers Tmax).
# Outputs (OUT_FIG/annex): fig_annex_conditional_shapley_LAD.{png,pdf}
#                          tab_annex_conditional_shapley_LAD.csv
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
source(here::here("R/h1_shapley_conditional.R"))      # provides .LAD_PAIRS
ANN <- file.path(PIPE$OUT_FIG, "annex"); dir.create(ANN, recursive=TRUE, showWarnings=FALSE)
cli_h1("STAGE 10 (annex) — conditional Shapley LAD | Hmax (off-manifold bound)")

coal <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))
clu  <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))
phi05 <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "shapley_perplot.rds")))
bits <- PIPE$BITS

grp <- vapply(.LAD_PAIRS, function(pp) pp$group, "")
wts <- vapply(.LAD_PAIRS, function(pp) pp$w, 0)
cond_one <- function(vec) {                            # vec named by bit (ΔTmax value)
  d <- vapply(.LAD_PAIRS, function(pp) as.numeric(vec[pp$S_LAD]) - as.numeric(vec[pp$S]), 0)
  c(global    = sum(wts * d),                          # = standard φ_LAD (my sign)
    Hmax_real = 2 * sum(wts[grp=="A"] * d[grp=="A"]),  # conditional, coherent
    Hmax_base = 2 * sum(wts[grp=="B"] * d[grp=="B"]))  # conditional, chimera
}

W <- dcast(coal, id_plot ~ bit, value.var = "Tmax_all")
res <- rbindlist(lapply(seq_len(nrow(W)), function(i) {
  v <- as.numeric(unlist(W[i, ..bits])); names(v) <- bits
  c3 <- cond_one(v)
  data.table(id_plot = W$id_plot[i], global=c3["global"],
             Hmax_real=c3["Hmax_real"], Hmax_base=c3["Hmax_base"])
}))
res <- merge(res, clu, by="id_plot")

# consistency: φ_LAD(global, conditional) == φ_LAD from stage 05 (interventional)
chk <- merge(res[, .(id_plot, global)],
             phi05[metric=="Tmax" & period=="All period" & trait=="LAD", .(id_plot, phi)], by="id_plot")
pipe_assert(all(pipe_near(chk$global, chk$phi, 1e-6)),
            "Conditional decomposition sums to stage-05 interventional φ_LAD")

# bootstrap median CI per (cluster or All) × quantity
resC <- rbind(res, copy(res)[, Cluster := "All"]); resC[, Cluster := factor(Cluster, levels=c(paste0("P",1:4),"All"))]
long <- melt(resC, id.vars=c("id_plot","Cluster"),
             measure.vars=c("global","Hmax_real","Hmax_base"),
             variable.name="estimand", value.name="phi")
boot_ci <- function(x){x<-x[is.finite(x)];n<-length(x); if(n<2) return(c(median=if(n)median(x) else NA,NA,NA))
  set.seed(42L); m<-vapply(1:2000,function(b) median(x[sample.int(n,n,TRUE)]),numeric(1))
  c(median=median(x), quantile(m,c(.025,.975),names=FALSE))}
ci <- long[, as.list(boot_ci(phi)), by=.(Cluster, estimand)]
setnames(ci, c("Cluster","estimand","median","ci_lo","ci_hi"))
ci[, estimand := factor(estimand, levels=c("global","Hmax_real","Hmax_base"),
       labels=c("phi[LAD]~(global)","phi[LAD]*'|'*H[max]==real","phi[LAD]*'|'*H[max]==baseline"))]
fwrite(ci, file.path(PIPE$OUT_TAB, "tab_annex_conditional_shapley_LAD.csv"))

cli_h2("LAD conditional Shapley (ΔTmax, all) — pooled All")
print(ci[Cluster=="All", .(estimand, median=round(median,3), CI=sprintf("[%.3f,%.3f]",ci_lo,ci_hi))])
gap <- ci[Cluster=="All"][grepl("real",estimand), median] - ci[Cluster=="All"][grepl("baseline",estimand), median]
cli_alert("Off-manifold gap (real − baseline, All): {round(gap,3)} °C")

# figure
ecol <- c("phi[LAD]~(global)"="grey25", "phi[LAD]*'|'*H[max]==real"="#1A9850", "phi[LAD]*'|'*H[max]==baseline"="#D7191C")
p <- ggplot(ci, aes(median, estimand, colour=estimand)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey55") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.25, linewidth=0.9) +
  geom_point(size=3.3) +
  facet_wrap(~ Cluster, nrow=1) +
  scale_y_discrete(labels=function(x) parse(text=x)) +
  scale_colour_manual(values=ecol, guide="none") +
  labs(x=expression(varphi[LAD]~"on "*Delta*T[max]~"("*degree*"C, median ± bootstrap 95% CI)"), y=NULL,
       title="ANNEX — Conditional Shapley for LAD: off-manifold bound",
       caption="Hmax=real = coherent domain; Hmax=baseline = real LAD compressed onto mean Hmax (chimera). global = mean of the two. Agreement ⇒ off-manifold artefact small.") +
  theme_bw(base_size=15) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold"), axis.text.y=element_text(size=13),
        plot.caption=element_text(size=10, colour="grey40"))
ggsave(file.path(ANN,"fig_annex_conditional_shapley_LAD.png"), p, width=16, height=5.5, dpi=300, bg="white")
ggsave(file.path(ANN,"fig_annex_conditional_shapley_LAD.pdf"), p, width=16, height=5.5, device=cairo_pdf)
cli_alert_success("Saved conditional Shapley LAD figure + table")
