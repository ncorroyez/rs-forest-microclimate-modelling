# ==============================================================================
# Comparison: MuSICA before_Klara (pre-3.2.3) vs modifs_Klara (3.2.3)
# Replays the main.R post-processing on both output sets and overlays figures.
# ==============================================================================

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

suppressPackageStartupMessages({
  library(ncdf4)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
})

invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# ---- Config ------------------------------------------------------------------

VERSIONS <- list(
  before_Klara = "out_files/before_Klara",
  modifs_Klara = "out_files/modifs_Klara"
)
REF_SC       <- "REF_all_real"
FORWARD_ORDER <- c("H1f_0_Null_baseline", "H1f_1_LAI_only",
                    "H1f_2_LAI_Hmax",      "H1f_3_LAI_Hmax_fCover",
                    "H1f_4_Full_real")
OUT_DIR <- "outputs/version_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load shared data --------------------------------------------------------

df_sample <- readRDS(file.path(CFG$out_dir, "clhs_sample.rds"))
df_macro  <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---- Helper: build named path vector for a version root ----------------------

build_paths <- function(root, scenarios_rel) {
  paths <- file.path(root, scenarios_rel)
  names(paths) <- basename(scenarios_rel)
  paths[dir.exists(paths)]
}

# Relative paths for each analysis block
REL_H1_FORWARD <- file.path("H1_forward",
  c("REF_all_real", "H1f_0_Null_baseline", "H1f_1_LAI_only",
    "H1f_2_LAI_Hmax", "H1f_3_LAI_Hmax_fCover", "H1f_4_Full_real"))

REL_H1_LOO <- file.path("H1_loo",
  c("H1l_dropLAI_mean", "H1l_dropHmax_mean",
    "H1l_dropfCover_mean", "H1l_dropLAD_uniform", "H1l_dropLAD_meanShape"))

REL_H2 <- file.path("H2_uniform_vs_real",
  c("H2_real_LAD", "H2_uniform_LAD", "H2_cluster_type_LAD"))

# ---- Extract ΔTmax for all versions ------------------------------------------

cat("Extracting scenarios...\n")

results <- imap(VERSIONS, function(root, vname) {
  cat(sprintf("  [%s]\n", vname))

  paths_h1 <- build_paths(root, c(REL_H1_FORWARD, REL_H1_LOO))
  paths_h2 <- build_paths(root, REL_H2)

  df_h1 <- extract_all_scenarios(paths_h1, df_macro, CFG$date_seq) %>%
    join_with_sample(df_sample)
  df_h2 <- extract_all_scenarios(paths_h2, df_macro, CFG$date_seq)

  list(
    h1_scores = score_scenarios_vs_reference(df_h1, REF_SC),
    df_h1     = df_h1,
    h2_wide   = summarise_h2(df_h2),
    df_h2     = df_h2
  )
})

saveRDS(results, file.path(OUT_DIR, "results_both_versions.rds"))

# ---- Figure 1 : Forward curve overlaid ---------------------------------------

df_fwd <- map_df(names(results), function(v) {
  results[[v]]$h1_scores %>%
    filter(scenario %in% FORWARD_ORDER) %>%
    mutate(
      version = v,
      step    = match(scenario, FORWARD_ORDER)
    )
})

p_fwd <- ggplot(df_fwd, aes(x = step, y = rmse, colour = version, group = version)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = seq_along(FORWARD_ORDER), labels = FORWARD_ORDER) +
  scale_colour_manual(values = c(before_Klara = "#d73027", modifs_Klara = "#31688e"),
                      labels = c(before_Klara = "before_Klara (v3.2.0)",
                                 modifs_Klara = "modifs_Klara (v3.2.3)")) +
  labs(title    = "H1 — Forward inclusion curve",
       subtitle = "RMSE vs reference REF_all_real",
       x = NULL, y = "RMSE (°C)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "h1_forward_curve_comparison.png"),
       p_fwd, width = 10, height = 5, dpi = 200)

# ---- Figure 2 : LOO scores side by side -------------------------------------

df_loo <- map_df(names(results), function(v) {
  results[[v]]$h1_scores %>%
    filter(str_starts(scenario, "H1l_")) %>%
    mutate(version = v)
})

p_loo <- ggplot(df_loo, aes(x = reorder(scenario, rmse), y = rmse,
                              fill = version, alpha = version)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", rmse)),
            position = position_dodge(0.7), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c(before_Klara = "#d73027", modifs_Klara = "#31688e")) +
  scale_alpha_manual(values = c(before_Klara = 0.75, modifs_Klara = 0.9)) +
  labs(title = "H1 — Leave-one-out: RMSE vs reference",
       x = NULL, y = "RMSE (°C)", fill = NULL, alpha = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(OUT_DIR, "h1_loo_comparison.png"),
       p_loo, width = 10, height = 5, dpi = 200)

# ---- Figure 3 : H2 mean diff (real - uniform) per version -------------------

df_h2_sum <- map_df(names(results), function(v) {
  results[[v]]$h2_wide %>%
    summarise(
      mean_diff = mean(diff, na.rm = TRUE),
      sd_diff   = sd(diff,   na.rm = TRUE),
      n         = n()
    ) %>%
    mutate(version = v, se = sd_diff / sqrt(n))
})

p_h2_bar <- ggplot(df_h2_sum, aes(x = version, y = mean_diff, fill = version)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_diff - se, ymax = mean_diff + se), width = 0.15) +
  geom_text(aes(label = sprintf("%.3f °C", mean_diff)),
            vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c(before_Klara = "#d73027", modifs_Klara = "#31688e")) +
  labs(title    = "H2 — Mean per-plot, per-day diff (Real − Uniform LAD)",
       subtitle = "Mesure de l'effet de la forme du profil LAD",
       x = NULL, y = "Δ ΔTmax moyen (°C)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# H2 diff distribution overlay
df_h2_long <- map_df(names(results), function(v) {
  results[[v]]$h2_wide %>% mutate(version = v)
})

p_h2_dens <- ggplot(df_h2_long, aes(x = diff, colour = version, fill = version)) +
  geom_density(alpha = 0.2, linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(before_Klara = "#d73027", modifs_Klara = "#31688e")) +
  scale_fill_manual(values  = c(before_Klara = "#d73027", modifs_Klara = "#31688e")) +
  labs(title  = "H2 — Distribution des diff (Real − Uniform)",
       x      = "Δ ΔTmax (°C)", y = "Densité",
       colour = NULL, fill = NULL) +
  theme_bw(base_size = 12)

p_h2 <- p_h2_bar + p_h2_dens +
  plot_annotation(title = "H2 — Impact de la version MuSICA sur l'effet profil LAD")

ggsave(file.path(OUT_DIR, "h2_comparison.png"),
       p_h2, width = 12, height = 5, dpi = 200)

# ---- Figure 4 : ΔTmax reference scenario — before vs after scatter ----------

df_ref_both <- map_df(names(results), function(v) {
  results[[v]]$df_h1 %>%
    filter(scenario == REF_SC) %>%
    mutate(version = v)
}) %>%
  pivot_wider(id_cols = c(x, y, date), names_from = version, values_from = Delta_Tmax)

r2_ref  <- cor(df_ref_both$before_Klara, df_ref_both$modifs_Klara, use = "complete.obs")^2
bias_ref <- mean(df_ref_both$modifs_Klara - df_ref_both$before_Klara, na.rm = TRUE)
rmse_ref <- sqrt(mean((df_ref_both$modifs_Klara - df_ref_both$before_Klara)^2, na.rm = TRUE))
lims_ref <- range(c(df_ref_both$before_Klara, df_ref_both$modifs_Klara), na.rm = TRUE)

p_ref_scatter <- ggplot(df_ref_both, aes(x = before_Klara, y = modifs_Klara)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(option = "C", trans = "log10", name = "Count\n(log)") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey20") +
  coord_fixed(xlim = lims_ref, ylim = lims_ref) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R² = %.4f\nBias = %.3f°C\nRMSE = %.3f°C",
                           r2_ref, bias_ref, rmse_ref),
           hjust = -0.1, vjust = 1.2, fontface = "bold", size = 4.5) +
  labs(title    = sprintf("ΔTmax scénario référence (%s)", REF_SC),
       subtitle = "Chaque point = 1 plot × 1 jour (été 2021)",
       x        = "ΔTmax before_Klara (°C)",
       y        = "ΔTmax modifs_Klara (°C)") +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_DIR, "ref_scenario_scatter.png"),
       p_ref_scatter, width = 7, height = 6, dpi = 200)

# ---- Summary table -----------------------------------------------------------

df_scores_all <- map_df(names(results), function(v) {
  results[[v]]$h1_scores %>% mutate(version = v)
})

cat("\n== H1 scores comparison ==\n")
df_scores_wide <- df_scores_all %>%
  pivot_wider(id_cols = scenario, names_from = version,
              values_from = c(rmse, mean_diff),
              names_glue = "{.value}_{version}") %>%
  mutate(delta_rmse = rmse_modifs_Klara - rmse_before_Klara) %>%
  arrange(abs(delta_rmse))

print(df_scores_wide, n = Inf)
write.csv(df_scores_wide, file.path(OUT_DIR, "h1_scores_comparison.csv"), row.names = FALSE)

cat(sprintf("\nH2 mean diff — before: %.3f°C  |  after: %.3f°C  |  delta: %.3f°C\n",
            df_h2_sum$mean_diff[df_h2_sum$version == "before_Klara"],
            df_h2_sum$mean_diff[df_h2_sum$version == "modifs_Klara"],
            diff(df_h2_sum$mean_diff)))

cat(sprintf("REF scatter — R²=%.4f  bias=%.3f°C  RMSE=%.3f°C\n",
            r2_ref, bias_ref, rmse_ref))

cat(sprintf("\nFigures saved to %s/\n", OUT_DIR))
