# ==============================================================================
# PIPELINE STAGE 07 — HOBO validation of REF (1111) — the "glass ceiling".
#   • per-plot mean ΔTmax sim vs obs  → r  (assert ≈ 0.92)
#   • buffering/amplifying split from REF sim log slope (fixed 1 m) → legacy assert 45/8 (v10 branch)
# Figures: per-plot mean scatter (coloured by cluster) + per-cluster facets.
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
source(here::here("scripts/_article_style.R"))   # PAL_CLUSTER, theme_article (style partagé, anti-dérive)
cli_h1("STAGE 07 — HOBO validation (REF)")

V   <- readRDS(file.path(PIPE$OUT_DATA, "ref_validation.rds"))
clu <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))

# ---- per-plot mean ΔTmax sim vs obs ----------------------------------------
sim_pp <- V$ref_daily[, .(Delta_sim = mean(Delta_sim, na.rm=TRUE)), by=id_plot]
obs_pp <- V$obs_daily[, .(Delta_obs = mean(Delta_obs, na.rm=TRUE)), by=id_plot]
pp <- Reduce(function(a,b) merge(a,b,by="id_plot"), list(sim_pp, obs_pp, clu))
r_pp   <- cor(pp$Delta_sim, pp$Delta_obs)
rmse_pp<- sqrt(mean((pp$Delta_sim - pp$Delta_obs)^2)); bias_pp <- mean(pp$Delta_sim - pp$Delta_obs)
cli_alert("Per-plot mean ΔTmax: r={round(r_pp,3)} RMSE={round(rmse_pp,2)}°C bias={round(bias_pp,2)}°C n={nrow(pp)}")

# ---- buffering / amplifying split (REF sim log slope, fixed 1 m) ------------
sl <- V$ref_slope
n_buf <- sum(sl$log_slope_sim < 0, na.rm=TRUE); n_amp <- sum(sl$log_slope_sim > 0, na.rm=TRUE)
cli_alert("REF sim slope split: {n_buf} buffering / {n_amp} amplifying (n={nrow(sl)})")

# ---- figures (style partagé ; pas de titre dans la figure = convention maison) --
pal <- PAL_CLUSTER
mae_pp <- mean(abs(pp$Delta_sim - pp$Delta_obs))
lim <- range(c(pp$Delta_sim, pp$Delta_obs)); lim <- lim + c(-1,1)*diff(lim)*0.06
pA <- ggplot(pp, aes(Delta_obs, Delta_sim, colour=Cluster, shape=Cluster)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55") +
  geom_hline(yintercept=0, linetype="dotted", colour="grey80") +
  geom_vline(xintercept=0, linetype="dotted", colour="grey80") +
  geom_point(size=3.4, alpha=0.9) +
  annotate("text", x=-Inf, y=Inf, hjust=-0.08, vjust=1.5, size=5,
           label=sprintf("italic(r)==%.2f~~RMSE==%.2f~degree*C~~MAE==%.2f~~bias==%+.2f", r_pp, rmse_pp, mae_pp, bias_pp), parse=TRUE) +
  coord_equal(xlim=lim, ylim=lim) + scale_colour_manual(values=pal) +
  labs(x=expression("Observed mean "*Delta*T[max]~"("*degree*"C)"),
       y=expression("Simulated mean "*Delta*T[max]~"("*degree*"C)")) +
  theme_article(18)
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_perplot.png"), pA, width=8, height=8, dpi=300, bg="white")
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_perplot.pdf"), pA, width=8, height=8, device=cairo_pdf)

# daily sim vs obs faceted by cluster
dd <- merge(merge(V$ref_daily, V$obs_daily, by=c("id_plot","date")), clu, by="id_plot")
annC <- dd[, .(r=cor(Delta_sim, Delta_obs), RMSE=sqrt(mean((Delta_sim-Delta_obs)^2)), n=.N), by=Cluster]
annC[, lab := sprintf("italic(r)==%.2f~~RMSE==%.2f", r, RMSE)]
limd <- range(c(dd$Delta_sim, dd$Delta_obs))
pB <- ggplot(dd, aes(Delta_obs, Delta_sim, colour=Cluster)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55") +
  geom_point(size=0.7, alpha=0.25) + facet_wrap(~Cluster, nrow=1) +
  geom_text(data=annC, aes(x=-Inf,y=Inf,label=lab), parse=TRUE, hjust=-0.06, vjust=1.5, size=4.2, inherit.aes=FALSE) +
  coord_equal(xlim=limd, ylim=limd) + scale_colour_manual(values=pal, guide="none") +
  labs(x=expression("Observed daily "*Delta*T[max]~"("*degree*"C)"),
       y=expression("Simulated daily "*Delta*T[max]~"("*degree*"C)"),
       title="REF daily ΔTmax sim vs obs, by cluster") +
  theme_bw(base_size=16) + theme(panel.grid.minor=element_blank(),
       strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"))
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_by_cluster.png"), pB, width=16, height=5, dpi=300, bg="white")
ggsave(file.path(PIPE$OUT_FIG,"fig_validation_by_cluster.pdf"), pB, width=16, height=5, device=cairo_pdf)

fwrite(pp,  file.path(PIPE$OUT_TAB, "tab_validation_perplot.csv"))
fwrite(annC, file.path(PIPE$OUT_TAB, "tab_validation_by_cluster.csv"))

# ---- regression asserts (branch-aware; enforced on legacy AND published z05) --
.vkey <- if (PIPE$BRANCH == "") "legacy" else PIPE$BRANCH
.tg   <- PIPE$ASSERT$VAL[[.vkey]]
if (!is.null(.tg)) {
  pipe_assert(pipe_near(r_pp, .tg$r_pp, 0.03), sprintf("[%s] Per-plot ΔTmax r = %.3f (target %.2f)", .vkey, r_pp, .tg$r_pp))
  pipe_assert(n_buf == .tg$buf && n_amp == .tg$amp,
              sprintf("[%s] Buffering/amplifying split = %d/%d (target %d/%d)", .vkey, n_buf, n_amp, .tg$buf, .tg$amp))
} else cli_alert_info("[branch {PIPE$BRANCH}] per-plot r = {round(r_pp,3)} ; split {n_buf}/{n_amp} (no target defined — not enforced)")
cli_alert_success("Stage 07 done.")
