# ==============================================================================
# PIPELINE STAGE 10b (ANNEX) — Conditional Shapley for LAI and fCover, bounding
# the OFF-MANIFOLD artefact that the LAI↔fCover collinearity (r ≈ 0.76) creates.
#
# Motivation: φ_LAI ≈ φ_fCover in the standard attribution, but LAI and fCover are
# collinear, so coalitions like "LAI=real high, fCover=baseline" (and vice-versa)
# are off the natural manifold. We CANNOT claim them independently co-dominant
# unless their attribution is stable whether or not the correlated partner is real.
#
# Same machinery as stage 10 (LAD|Hmax). We decompose each trait's φ by the level
# of its correlated partner in the coalition S:
#   φ_LAI    = mean( φ_LAI | fCover=real ,  φ_LAI | fCover=baseline )
#   φ_fCover = mean( φ_fCover | LAI=real ,  φ_fCover | LAI=baseline )
# "real" = coherent / on-manifold ; "baseline" = chimera / off-manifold.
# Agreement of the two conditionals ⇒ the off-manifold artefact is small ⇒ the
# trait's contribution is robust and the LAI/fCover split is defensible. Large gap
# ⇒ the split is labile and must be reported as a combined quantity/cover factor.
#
# Metric = ΔTmax (all summer, 1 m). Bit order: LAI | Hmax | fCover | LAD.
# Sign: φ<0 = trait buffers (lowers Tmax). d = value(S∪{trait}) − value(S).
# Outputs (OUT_FIG/annex): fig_annex_conditional_shapley_lai_fcover.{png,pdf}
#                          tab_annex_conditional_shapley_lai_fcover.csv
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
ANN <- file.path(PIPE$OUT_FIG, "annex"); dir.create(ANN, recursive=TRUE, showWarnings=FALSE)
cli_h1("STAGE 10b (annex) — conditional Shapley LAI|fCover & fCover|LAI (off-manifold bound)")

coal  <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))
clu   <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))
phi05 <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "shapley_perplot.rds")))
bits  <- PIPE$BITS

# Shapley weights w(|S|) for n=4 : s=0→1/4, s=1→1/12, s=2→1/12, s=3→1/4.
.w_of <- function(S) { s <- sum(utf8ToInt(S) == utf8ToInt("1")); c(1/4,1/12,1/12,1/4)[s+1] }

# Build the 8 pairs for toggling `tog_pos` (1=LAI .. 4=LAD), conditioned on the
# level of `cond_pos` in S : group "A" = partner real (bit=1), "B" = baseline.
make_pairs <- function(tog_pos, cond_pos) {
  Svec <- bits[substr(bits, tog_pos, tog_pos) == "0"]          # coalitions without the toggled trait
  lapply(Svec, function(S) {
    S_on <- S; substr(S_on, tog_pos, tog_pos) <- "1"
    list(S = S, S_on = S_on, w = .w_of(S),
         group = if (substr(S, cond_pos, cond_pos) == "1") "A" else "B")
  })
}
.LAI_PAIRS <- make_pairs(tog_pos = 1, cond_pos = 3)   # φ_LAI conditioned on fCover (pos 3)
.FCV_PAIRS <- make_pairs(tog_pos = 3, cond_pos = 1)   # φ_fCover conditioned on LAI  (pos 1)

# φ decomposition for one plot's 16-value vector (named by bit).
cond_one <- function(vec, PAIRS) {
  grp <- vapply(PAIRS, function(pp) pp$group, "")
  wts <- vapply(PAIRS, function(pp) pp$w, 0)
  d   <- vapply(PAIRS, function(pp) as.numeric(vec[pp$S_on]) - as.numeric(vec[pp$S]), 0)
  c(global       = sum(wts * d),                       # = standard interventional φ
    partner_real = 2 * sum(wts[grp=="A"] * d[grp=="A"]),
    partner_base = 2 * sum(wts[grp=="B"] * d[grp=="B"]))
}

W <- dcast(coal, id_plot ~ bit, value.var = "Tmax_all")
decomp <- function(PAIRS, trait_name) {
  out <- rbindlist(lapply(seq_len(nrow(W)), function(i) {
    v <- as.numeric(unlist(W[i, ..bits])); names(v) <- bits
    c3 <- cond_one(v, PAIRS)
    data.table(id_plot=W$id_plot[i], trait=trait_name,
               global=c3["global"], partner_real=c3["partner_real"], partner_base=c3["partner_base"])
  }))
  merge(out, clu, by="id_plot")
}
resL <- decomp(.LAI_PAIRS, "LAI")
resF <- decomp(.FCV_PAIRS, "fCover")

# ---- consistency: global conditional φ == stage-05 interventional φ ----------
# period stored as factor "All period"/"10% hottest"; coerce robustly.
phi05[, period_chr := as.character(period)]
all_lbl <- if ("All period" %in% phi05$period_chr) "All period" else sort(unique(phi05$period_chr))[1]
chk_one <- function(res, tr) {
  ref <- phi05[metric=="Tmax" & period_chr==all_lbl & trait==tr, .(id_plot, phi)]
  m <- merge(res[, .(id_plot, global)], ref, by="id_plot")
  pipe_assert(nrow(m) > 0 && all(pipe_near(m$global, m$phi, 1e-5)),
              sprintf("Conditional decomposition sums to stage-05 interventional φ_%s", tr))
}
chk_one(resL, "LAI"); chk_one(resF, "fCover")

# ---- bootstrap median CI per (cluster|All) × estimand ------------------------
boot_ci <- function(x){ x<-x[is.finite(x)]; n<-length(x)
  if(n<2) return(c(median=if(n)median(x) else NA, NA, NA))
  set.seed(42L); m<-vapply(1:2000, function(b) median(x[sample.int(n,n,TRUE)]), numeric(1))
  c(median=median(x), quantile(m, c(.025,.975), names=FALSE)) }

build_ci <- function(res) {
  resC <- rbind(res, copy(res)[, Cluster := "All"])
  resC[, Cluster := factor(Cluster, levels=c(paste0("P",1:4),"All"))]
  long <- melt(resC, id.vars=c("id_plot","trait","Cluster"),
               measure.vars=c("global","partner_real","partner_base"),
               variable.name="estimand", value.name="phi")
  ci <- long[, as.list(boot_ci(phi)), by=.(trait, Cluster, estimand)]
  setnames(ci, c("trait","Cluster","estimand","median","ci_lo","ci_hi")); ci
}
ci <- rbind(build_ci(resL), build_ci(resF))
fwrite(ci, file.path(PIPE$OUT_TAB, "tab_annex_conditional_shapley_lai_fcover.csv"))

# ---- console summary + off-manifold gaps -------------------------------------
report <- function(tr, partner) {
  z <- ci[trait==tr & Cluster=="All"]
  cli_h2(sprintf("%s conditional Shapley (ΔTmax, all) — pooled All [partner = %s]", tr, partner))
  print(z[, .(estimand, median=round(median,3), CI=sprintf("[%.3f,%.3f]", ci_lo, ci_hi))])
  gap <- z[estimand=="partner_real", median] - z[estimand=="partner_base", median]
  cli_alert("{tr}: off-manifold gap (partner real − baseline, All) = {round(gap,3)} °C")
  invisible(gap)
}
gL <- report("LAI", "fCover")
gF <- report("fCover", "LAI")
# on-manifold (partner real) comparison — the defensible LAI-vs-fCover contrast
rL <- ci[trait=="LAI"    & Cluster=="All" & estimand=="partner_real", median]
rF <- ci[trait=="fCover" & Cluster=="All" & estimand=="partner_real", median]
cli_alert_info("On-manifold (partner real): φ_LAI={round(rL,3)} vs φ_fCover={round(rF,3)} °C")

# ---- figure : two panels (LAI|fCover, fCover|LAI), pooled All + clusters -----
labmap <- c(global="varphi~(global)", partner_real="varphi*'|'*partner==real",
            partner_base="varphi*'|'*partner==baseline")
ci[, est_lab := factor(labmap[as.character(estimand)],
       levels=c("varphi~(global)","varphi*'|'*partner==real","varphi*'|'*partner==baseline"))]
ci[, trait_lab := factor(trait, levels=c("LAI","fCover"),
       labels=c("varphi[LAI]*'  | fCover'", "varphi[fCover]*'  | LAI'"))]
ecol <- c("varphi~(global)"="grey25",
          "varphi*'|'*partner==real"="#1A9850",
          "varphi*'|'*partner==baseline"="#D7191C")
p <- ggplot(ci[Cluster=="All"], aes(median, est_lab, colour=est_lab)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey55") +
  geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.22, linewidth=0.9) +
  geom_point(size=3.3) +
  facet_wrap(~ trait_lab, ncol=1, scales="free_y", labeller=label_parsed) +
  scale_y_discrete(labels=function(x) parse(text=x)) +
  scale_colour_manual(values=ecol, guide="none") +
  labs(x=expression(varphi~"on "*Delta*T[max]~"("*degree*"C, median ± bootstrap 95% CI)"), y=NULL,
       title="Conditional Shapley for LAI & fCover: off-manifold bound under collinearity",
       caption=paste(strwrap(paste(
         "partner = real: coherent / on-manifold (correlated trait present).",
         "partner = baseline: chimera (correlated trait stripped to baseline).",
         "global = mean of the two. Agreement of the two conditionals ⇒ collinearity (off-manifold)",
         "artefact small ⇒ each trait's contribution is robust and the LAI/fCover split is defensible."),
         width=110), collapse="\n")) +
  theme_bw(base_size=14) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=13), axis.text.y=element_text(size=12),
        plot.title=element_text(size=13, face="bold"),
        plot.caption=element_text(size=9, colour="grey40", hjust=0))
ggsave(file.path(ANN,"fig_annex_conditional_shapley_lai_fcover.png"), p, width=12, height=7, dpi=300, bg="white")
ggsave(file.path(ANN,"fig_annex_conditional_shapley_lai_fcover.pdf"), p, width=12, height=7, device=cairo_pdf)
cli_alert_success("Saved conditional Shapley LAI|fCover figure + table")
