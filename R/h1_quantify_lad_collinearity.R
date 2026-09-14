# ==============================================================================
# Quantify LAD scalar collinearity with LAI / Hmax / fCover.
# LAD scalars : FPC1 (from joint cLHS+HOBO FPCA) and LAD_L2 (deviation from uniform).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

OUT_TAB <- here::here("outputs/lovb_floor05/tables")
OUT_FIG <- here::here("outputs/lovb_floor05/figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ---- cLHS sample (n=400) ----------------------------------------------------
df_c <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds")))
cli_h1("LAD collinearity in cLHS sample (n = 400)")

vars_to_test <- c("LAI", "Hmax", "fCover")
res_clhs <- list()
for (v in vars_to_test) {
  for (lad_scalar in c("FPC1")) {
    res_clhs[[length(res_clhs)+1L]] <- data.table(
      scope = "cLHS (n=400)", LAD_scalar = lad_scalar, trait = v,
      pearson = round(cor(df_c[[lad_scalar]], df_c[[v]],
                            use = "complete.obs"), 3),
      spearman = round(cor(df_c[[lad_scalar]], df_c[[v]],
                              method = "spearman", use = "complete.obs"), 3)
    )
  }
}

# Multiple regression : LAD scalar = b1*LAI + b2*Hmax + b3*fCover + b0
fit_clhs <- lm(FPC1 ~ LAI + Hmax + fCover, data = df_c)
sum_clhs <- summary(fit_clhs)
cli_alert("Multiple regression  FPC1 ~ LAI + Hmax + fCover  (cLHS) :")
print(round(coef(sum_clhs), 4))
cli_alert("Adjusted R² : {round(sum_clhs$adj.r.squared, 3)}")
cli_alert("F-stat p : {format.pval(pf(sum_clhs$fstatistic[1], sum_clhs$fstatistic[2], sum_clhs$fstatistic[3], lower.tail=FALSE))}")

# ---- HOBO sample (n=53) -----------------------------------------------------
DT_h <- readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds"))
cli_h1("LAD collinearity in HOBO sample (n = 53)")
for (v in vars_to_test) {
  for (lad_scalar in c("FPC1", "LAD_L2")) {
    res_clhs[[length(res_clhs)+1L]] <- data.table(
      scope = "HOBO (n=53)", LAD_scalar = lad_scalar, trait = v,
      pearson = round(cor(DT_h[[lad_scalar]], DT_h[[v]],
                            use = "complete.obs"), 3),
      spearman = round(cor(DT_h[[lad_scalar]], DT_h[[v]],
                              method = "spearman", use = "complete.obs"), 3)
    )
  }
}
tab_corr <- rbindlist(res_clhs)
cli_alert("Pairwise correlations LAD scalar vs trait :")
print(tab_corr)

fit_hobo_fpc1 <- lm(FPC1 ~ LAI + Hmax + fCover, data = DT_h)
sum_hobo_fpc1 <- summary(fit_hobo_fpc1)
cli_alert("Multiple regression  FPC1 ~ LAI + Hmax + fCover  (HOBO) :")
print(round(coef(sum_hobo_fpc1), 4))
cli_alert("Adjusted R² : {round(sum_hobo_fpc1$adj.r.squared, 3)}  |  F-stat p : {format.pval(pf(sum_hobo_fpc1$fstatistic[1], sum_hobo_fpc1$fstatistic[2], sum_hobo_fpc1$fstatistic[3], lower.tail=FALSE))}")

fit_hobo_l2 <- lm(LAD_L2 ~ LAI + Hmax + fCover, data = DT_h)
sum_hobo_l2 <- summary(fit_hobo_l2)
cli_alert("Multiple regression  LAD_L2 ~ LAI + Hmax + fCover  (HOBO) :")
print(round(coef(sum_hobo_l2), 4))
cli_alert("Adjusted R² : {round(sum_hobo_l2$adj.r.squared, 3)}  |  F-stat p : {format.pval(pf(sum_hobo_l2$fstatistic[1], sum_hobo_l2$fstatistic[2], sum_hobo_l2$fstatistic[3], lower.tail=FALSE))}")

# ---- Quantify leakage : how much of LVA(LAD) signal at HOBO is explained by
#       the implicit LAI/fCover information leaking through LAD ?
cli_h1("Quantify LVA(LAD) signal leakage")
DT_h_full <- merge(DT_h,
                     readRDS(here::here("outputs/lovb_floor05/data/DT_contrib_HOBO_floor05.rds")),
                     by = "id_plot")
DT_h_full[, LVA_LAD := Tmax_mean_LVA_LAD - Tmax_mean_NULL]
cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
DT_h_full <- merge(DT_h_full, hobo_t, by = "id_plot")
DT_h_full[, dT_signed := -dTmax_obs]

# Raw correlation : LVA(LAD) vs observed
rho_raw <- cor(DT_h_full$dT_signed, DT_h_full$LVA_LAD, method = "spearman",
                 use = "complete.obs")
cli_alert("Raw rho(LVA(LAD), dT_signed) = {round(rho_raw, 3)}")

# Partial correlation : LVA(LAD) vs observed, controlling for LAI + Hmax + fCover
# Use complete-case subset to ensure aligned residuals
ok <- complete.cases(DT_h_full[, .(LVA_LAD, dT_signed, LAI, Hmax, fCover)])
sub <- DT_h_full[ok]
res_lva <- residuals(lm(LVA_LAD ~ LAI + Hmax + fCover, data = sub))
res_dT  <- residuals(lm(dT_signed ~ LAI + Hmax + fCover, data = sub))
rho_partial <- cor(res_dT, res_lva, method = "spearman")
cli_alert("Partial rho(LVA(LAD), dT_signed | LAI+Hmax+fCover) = {round(rho_partial, 3)}")
cli_alert("Leakage : {round((rho_raw - rho_partial)/rho_raw * 100, 1)}% of LVA(LAD) signal removed by controlling for LAI+Hmax+fCover")

# Save summary
tab_summary <- data.table(
  metric = c("rho_raw_LVA_LAD_vs_dT",
              "rho_partial_LVA_LAD_vs_dT_given_LAI_Hmax_fCover",
              "leakage_pct"),
  value = c(round(rho_raw, 3), round(rho_partial, 3),
              round((rho_raw - rho_partial) / rho_raw * 100, 1))
)
fwrite(tab_corr, file.path(OUT_TAB, "tab_LAD_collinearity.csv"))
fwrite(tab_summary, file.path(OUT_TAB, "tab_LVA_LAD_leakage.csv"))
cli_alert_success("Saved tab_LAD_collinearity.csv and tab_LVA_LAD_leakage.csv")

# ---- Figure : pairwise scatter LAD vs traits --------------------------------
cli_h2("Figure : LAD scalar vs traits")
df_c[, scope := "cLHS (n=400)"]
DT_h[, scope := "HOBO (n=53)"]
long <- rbindlist(list(
  df_c[, .(scope, LAI, Hmax, fCover, FPC1)],
  DT_h[, .(scope, LAI, Hmax, fCover, FPC1)]
), use.names = TRUE)
long_melted <- melt(long, id.vars = c("scope", "FPC1"),
                      measure.vars = c("LAI","Hmax","fCover"),
                      variable.name = "trait", value.name = "trait_value")

# Compute correlations for annotation
ann <- long_melted[, .(rho = round(cor(FPC1, trait_value, method = "spearman",
                                            use = "complete.obs"), 2)),
                      by = .(scope, trait)]
ann[, label := sprintf("rho == %+.2f", rho)]
ann[, scope := factor(scope, levels = c("cLHS (n=400)", "HOBO (n=53)"))]
ann[, trait := factor(trait, levels = c("LAI","Hmax","fCover"))]

p <- ggplot(long_melted, aes(x = trait_value, y = FPC1)) +
  geom_point(alpha = 0.45, colour = "#2C5F2D", size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400", alpha = 0.2) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
             parse = TRUE, hjust = -0.1, vjust = 1.3, size = 4.5,
             fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
  facet_grid(scope ~ trait, scales = "free", switch = "y") +
  labs(
    title    = "LAD scalar (FPC1) vs structural traits — collinearity",
    subtitle = "Spearman ρ annotated per panel    |    Linear fit + 95% CI shown for reference",
    x = "Trait value",
    y = "LAD shape (FPC1 score)"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text = element_text(face = "bold", size = 12),
         plot.title = element_text(face = "bold", size = 14),
         plot.subtitle = element_text(colour = "grey25", size = 11),
         panel.grid.minor = element_blank())
ggsave(file.path(OUT_FIG, "fig_LAD_collinearity.png"),
        p, width = 13, height = 7.5, dpi = 300)
cli_alert_success("Saved fig_LAD_collinearity.png")
