# ==============================================================================
# Le plafond de verre est-il de la topographie / hétérogénéité que MuSICA ignore ?
#
# MuSICA (1-D placette) ne voit ni relief, ni TWI, ni hétérogénéité latérale, ni
# lisières. On teste si le RÉSIDU de validation (ΔTmax_obs HOBO − ΔTmax_sim REF)
# par placette est expliqué par des covariables topo/structure (NC_Full, Oak_Only) :
#   élévation (dtm), pente, TWI, northness, SD hauteur (std), rumple, gap_fraction.
# + comparaison des 8 placettes amplificatrices vs 45 tamponnantes sur ces covariables.
#
# Reads :
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv  (stage A1)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv      (stage A2: fCover)
#         ~/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked/{dtm,slope_res_10_m,twi_res_10_m,northness_res_10_m,std_res_10_m,rumple_res_10_m}.tif
#         the six topographic covariates. OUTSIDE the repository (171 GB; Not_Masked, NOT
#         Oak_Only, because the oak mask drops the amplifying loggers). Absent -> the script stops.
# Writes: outputs/figures_pipeline/annex/fig_residual_vs_topo.{png,pdf}
#          + tab_residual_vs_topo.csv
# ==============================================================================

suppressMessages({ library(here); library(data.table); library(terra); library(tidyverse) })
source(here::here("scripts/_article_style.R"))
OUT <- here::here("outputs/figures_pipeline/annex"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NCF <- path.expand("~/Documents/NC_Full")
# Covariates from Not_Masked (NOT Oak_Only): the oak mask drops loggers on
# non-oak / gap / edge pixels — i.e. exactly the amplifying plots we must explain.
# Using Not_Masked keeps all 53 loggers (no NA) and characterizes the actual
# topographic / heterogeneity neighbourhood each logger sits in.
NM  <- file.path(NCF, "03_RESULTS/Blois/Metrics/Not_Masked")

# ---- 1. résidu de validation par placette + coords ---------------------------
# wind-corrected baseline: per-plot ΔTmax from the corrected full-model validation
# LINEAGE FIX 2026-07-29: this read tab_hobo_frblo_validation.csv, i.e. the superseded
# `frblo` lineage, so the validation residual decomposed here was NOT the residual of the
# chapter's reported model. The native20 validation is the canonical one (and is on the
# aligned clock); the two residuals correlate at r = 0.988 but their means differ,
# -0.724 against -0.980, with a maximum per-plot difference of 1.66 degC.
V <- fread(here::here("out_files/Chapter1/tables/tab_hobo_native20_validation.csv"))
res <- V[is.finite(sim_dt) & is.finite(obs_dt),
  .(id_plot, dTmax_sim = sim_dt, dTmax_obs = obs_dt)][
  , residual := dTmax_obs - dTmax_sim]                       # >0 : HOBO + chaud que le modèle
lad <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))[, .(id_plot, x, y)]
D <- merge(res, lad, by = "id_plot")

# ---- 2. covariables topo / hétérogénéité (NC_Full) ---------------------------
cov_files <- c(
  elevation   = file.path(NM, "dtm.tif"),
  slope       = file.path(NM, "slope_res_10_m.tif"),
  twi         = file.path(NM, "twi_res_10_m.tif"),
  northness   = file.path(NM, "northness_res_10_m.tif"),
  sd_height   = file.path(NM, "std_res_10_m.tif"),
  rumple      = file.path(NM, "rumple_res_10_m.tif")
)
# FAIL LOUDLY 2026-07-31. This used to drop missing covariates silently, so a run with
# ~/Documents/NC_Full unmounted produced a decomposition over whatever happened to be
# present and reported it as if complete. NC_Full is 171 GB, lives outside the repo and
# is documented nowhere; the failure mode was a quietly incomplete Appendix C.
.missing <- names(cov_files)[!file.exists(cov_files)]
if (length(.missing))
  stop("make_residual_vs_topo: missing covariate raster(s): ", paste(.missing, collapse = ", "),
       "\n  expected under: ", NM,
       "\n  This directory is outside the repository. Mount it, or the residual decomposition",
       "\n  of Appendix C would be computed over an incomplete covariate set.")
cat("Covariables trouvées :", paste(names(cov_files), collapse = ", "), "\n")

# points placettes -> CRS de chaque raster (coords supposées UTM31N / EPSG:32631)
pts <- vect(as.data.frame(D[, .(x, y, id_plot)]), geom = c("x","y"), crs = "EPSG:32631")
for (nm in names(cov_files)) {
  r <- rast(cov_files[[nm]])
  p <- if (crs(r) != crs(pts)) project(pts, crs(r)) else pts
  D[[nm]] <- terra::extract(r, p, ID = FALSE)[[1]]
}
# Fractional cover actually used to drive MuSICA (floored at 0.5) — the model's own
# openness variable, in place of a separate raster gap fraction.
D <- merge(D, fread(here::here("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv"))[, .(id_plot, fCover)], by = "id_plot")
D[, fCover := pmax(fCover, 0.5)]
fwrite(D, file.path(OUT, "tab_residual_vs_topo.csv"))
# Guard: every logger must have complete covariates (no silent lm() row drop).
cc <- complete.cases(as.data.frame(D)[, names(cov_files)])
cat(sprintf("Plots with complete covariates: %d / %d\n", sum(cc), nrow(D)))
if (any(!cc)) cat("  WARNING — plots still NA:", paste(D$id_plot[!cc], collapse=", "), "\n")

# ---- 3. corrélations + régression multiple -----------------------------------
covs <- c(names(cov_files), "fCover")
cors <- sapply(covs, function(c) {
  ok <- is.finite(D[[c]]) & is.finite(D$residual)
  if (sum(ok) < 5) return(c(r = NA, p = NA))
  ct <- cor.test(D[[c]][ok], D$residual[ok]); c(r = unname(ct$estimate), p = ct$p.value)
})
cat("\n=== corrélation covariable ↔ résidu (HOBO − modèle) ===\n")
print(round(t(cors), 3))
Dz <- copy(D); for (c in covs) Dz[[c]] <- scale(Dz[[c]])[,1]   # standardisé
fit <- lm(reformulate(covs, "residual"), data = Dz)
cat(sprintf("\nRégression multiple résidu ~ topo+hétéro : R²=%.2f (ajusté %.2f), p=%.3g\n",
            summary(fit)$r.squared, summary(fit)$adj.r.squared,
            pf(summary(fit)$fstatistic[1], summary(fit)$fstatistic[2], summary(fit)$fstatistic[3], lower.tail = FALSE)))
print(round(summary(fit)$coefficients, 3))

# Variance-inflation factors, quoted in the Fig. C1 caption (fCover ~ 1.3; sd_height
# and rumple mutually collinear). Added 2026-07-31: the caption cited VIFs that no
# script computed, so they could not be checked against the covariate set actually used.
# Computed directly as 1 / (1 - R2_j) of each covariate on the others, no extra package.
vifs <- sapply(covs, function(c)
  1 / (1 - summary(lm(reformulate(setdiff(covs, c), c), data = Dz))$r.squared))
cat("\n=== variance-inflation factors (Fig. C1 caption) ===\n")
print(round(vifs, 2))
fwrite(data.table(covariate = names(vifs), vif = round(unname(vifs), 3)),
       file.path(OUT, "tab_residual_vif.csv"))

# ---- 4. amplificatrices (slope micro-macro > 1) vs tamponnantes ---------------
slp <- V[, .(id_plot, slope_obs = obs_sl)]
D <- merge(D, slp, by = "id_plot"); D[, grp := ifelse(slope_obs > 1, "amplifying", "buffering")]
cat(sprintf("\nGroupes terrain : %d amplifying / %d buffering\n", sum(D$grp=="amplifying"), sum(D$grp=="buffering")))
cmp <- rbindlist(lapply(covs, function(c) {
  a <- D[grp=="amplifying"][[c]]; b <- D[grp=="buffering"][[c]]
  data.table(cov=c, amplifying=mean(a,na.rm=T), buffering=mean(b,na.rm=T),
             p_ttest = tryCatch(t.test(a,b)$p.value, error=function(e) NA)) }))
cat("\n=== moyennes amplifying vs buffering ===\n"); print(round_df <- cmp[, lapply(.SD, function(v) if(is.numeric(v)) round(v,3) else v)])

# ---- 5. figure : résidu vs fractional cover (variable du modèle) + groupes ----
best <- "fCover"
labf <- sprintf("residual ~ fractional cover: r=%.2f (p=%.3f) | full model adj. R²=%.2f (n=%d, %d covariates)",
                cors["r", best], cors["p", best],
                summary(fit)$adj.r.squared, nobs(fit), length(covs))
p1 <- ggplot(D, aes(.data[[best]], residual, colour = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = PAL_GRP, name = NULL) +
  labs(x = "Fractional cover (model input, floored at 0.5)",
       y = "Validation residual  (observed − model ΔTmax, °C)", subtitle = labf) +
  theme_article() + legend_corner()
ggsave_article(file.path(OUT, "fig_residual_vs_topo"), p1, 7.5, 5.5)
saveRDS(p1, file.path(OUT, "_rds_A8a_topo.rds"))      # pour le composite A8
cli::cli_alert_success("Saved fig_residual_vs_topo (png+pdf) + table")
