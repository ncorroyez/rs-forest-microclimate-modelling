# ==============================================================================
# compare_musica_versions_hobo.R
#
# Compare MuSICA < 3.2.3 (Blois binary) vs 3.2.3 (model-3.2.3 binary)
# on HOBO sensor locations using the REF_all_real scenario.
#
# Outputs → out_files/version_comparison_hobo/
#   v_old/REF_all_real/*.nc        — simulations with Blois binary
#   v3.2.3/REF_all_real/*.nc       — simulations with 3.2.3 binary
#   figures/comp_01_scatter.png    — ΔTmax obs vs sim scatter (hex density)
#   figures/comp_02_metrics.png    — R² / RMSE / bias bar chart
#   figures/comp_03_timeseries.png — Timeseries Tmax sur 3 HOBOs représentatifs
#   figures/comp_04_scatter_tmax.png — Scatter Tmax absolu obs vs sim
# ==============================================================================

setwd("/home/corroyez/Documents/z_Example_rmusica_31012025")

suppressPackageStartupMessages({
  library(musica.tools)
  library(rmusica)
  library(ncdf4)
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(viridis)
  library(sf)
  library(terra)
})

source("R/config.R")
source("R/io.R")
source("R/forest.R")
source("R/lad.R")
source("R/scenarios.R")
source("R/musica.R")
source("R/validation.R")

save_plot <- function(p, path, width = 12, height = 7, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
}

# ---- Config ------------------------------------------------------------------

# Blois site constants applied via musica.param to override version-specific NML defaults.
# Both NMLs (Blois & 3.2.3) have wrong FR-Bil coordinates (lat=44.49°N) and N_LEAF_AGE=3.
# Species/soil always point to Blois in_files. ABL_flag='none' disables yoyo in v3.2.3
# (its NML has ABL_flag='up'; the Blois binary predates ABL_flag so we don't pass it there).
.BLOIS_SETUP <- list(
  "site_latitude"    = 47.57,
  "site_longitude"   = 1.26,
  "site_altitude"    = 39.18,
  "forcing_timestep" = 3600,
  "n_leaf_age"       = 1,
  "param_files_path" = '"./in_files/Blois/in_files/"',
  "species_names"    = '"musica_veg1"',
  "soil_file_name"   = '"./in_files/Blois/in_files/musica_soil.nml"',
  "history_variables" = c("layer_thickness", "z_soil", "dz_soil", "t_soil",
                          "w_soil", "t_air", "w_air", "wind")
)

VERSIONS <- list(
  list(
    label       = "v_old",
    cmd         = normalizePath("in_files/Blois/musica"),
    musica_nml  = normalizePath("in_files/Blois/musica.nml"),
    musica_vars = normalizePath("in_files/Blois/variables.csv"),
    extra_setup = .BLOIS_SETUP,
    colour      = "#d8576b",
    linetype    = "dashed"
  ),
  list(
    label       = "v3.2.3",
    cmd         = normalizePath("in_files/model-3.2.3/musica"),
    musica_nml  = normalizePath("in_files/model-3.2.3/musica.nml"),
    musica_vars = normalizePath("in_files/model-3.2.3/variables.csv"),
    extra_setup = c(.BLOIS_SETUP, list("abl_flag" = '"none"')),
    colour      = "#31688e",
    linetype    = "solid"
  )
)

cat("--- Config check ---\n")
for (v in VERSIONS) {
  for (f in c(v$cmd, v$musica_nml, v$musica_vars)) {
    if (!file.exists(f)) stop("File not found: ", f)
  }
  if (file.access(v$cmd, mode = 1) != 0) {
    Sys.chmod(v$cmd, mode = "0755")
    cat(sprintf("  [chmod +x] %s\n", basename(v$cmd)))
  }
  cat(sprintf("[%s] binary=%s  nml=%s\n",
              v$label, basename(v$cmd), basename(v$musica_nml)))
}
cat("--------------------\n\n")

OUT_ROOT <- "out_files/version_comparison_hobo"
FIG_DIR  <- file.path(OUT_ROOT, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs -------------------------------------------------------------

cat("[1/4] Loading inputs...\n")
rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_macro       <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                   df_macro, CFG$ids_to_remove)
df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                     CFG$ids_to_remove)
df_sample      <- readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds")
ref_sc         <- scenario_reference(df_sample)

cat(sprintf("    %d HOBO plots | %d obs days | REF scenario: %s\n",
            nrow(df_hobo_inputs), nrow(df_macro), ref_sc$name))

# ---- Run simulations for each version ----------------------------------------

cat("[2/4] Running simulations...\n")

run_version <- function(v) {
  out_dir <- file.path(OUT_ROOT, v$label, ref_sc$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  existing <- list.files(out_dir, pattern = "\\.nc$")
  n_needed <- nrow(df_hobo_inputs)

  if (length(existing) >= n_needed) {
    cat(sprintf("  [%s] %d NC found — skipping\n", v$label, length(existing)))
    return(out_dir)
  }

  cat(sprintf("  [%s] %d/%d NC missing — running...\n",
              v$label, n_needed - length(existing), n_needed))

  for (i in seq_len(nrow(df_hobo_inputs))) {
    plot_row <- df_hobo_inputs[i, ]
    out_file <- file.path(out_dir,
                          sprintf("musica_out_HOBO_%s.nc", plot_row$id_plot))
    if (file.exists(out_file)) next
    cat(sprintf("    [%s] HOBO %s\n", v$label, plot_row$id_plot))
    run_musica_one(plot_row, ref_sc, out_file, CFG$forcing_file, v$cmd,
                   musica_nml       = v$musica_nml,
                   musica_variables = v$musica_vars,
                   extra_setup      = v$extra_setup)
  }
  out_dir
}

sim_dirs <- setNames(
  lapply(VERSIONS, run_version),
  sapply(VERSIONS, `[[`, "label")
)

# ---- Extract ΔTmax -----------------------------------------------------------

cat("[3/4] Extracting ΔTmax...\n")

extract_version <- function(v, out_dir) {
  nc_files <- list.files(out_dir, pattern = "\\.nc$", full.names = TRUE)
  map_dfr(nc_files, function(f) {
    id_str <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)")
    res    <- extract_deltatmax_one(f, df_macro, CFG$date_seq)
    if (is.null(res)) return(NULL)
    res %>% mutate(id_plot = id_str, version = v$label)
  })
}

df_sims <- map2_dfr(VERSIONS, sim_dirs, extract_version)

df_paired <- df_sims %>%
  inner_join(df_hobo_daily %>% dplyr::select(id_plot, date, Delta_obs, Tmax_obs),
             by = c("id_plot", "date")) %>%
  rename(Delta_sim = Delta_Tmax,
         Tmax_sim  = Tmax_micro) %>%
  mutate(version = factor(version, levels = sapply(VERSIONS, `[[`, "label")))

# ---- Metrics -----------------------------------------------------------------

df_metrics <- df_paired %>%
  group_by(version) %>%
  summarise(
    n    = n(),
    r2   = cor(Delta_obs, Delta_sim, use = "complete.obs")^2,
    rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)),
    bias = mean(Delta_sim - Delta_obs, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Metrics ---\n")
print(as.data.frame(df_metrics))
write.csv(df_metrics, file.path(OUT_ROOT, "version_metrics.csv"), row.names = FALSE)

# ---- Figure 1: Scatter obs vs sim per version --------------------------------

cat("[4/4] Generating figures...\n")

pal     <- setNames(sapply(VERSIONS, `[[`, "colour"), sapply(VERSIONS, `[[`, "label"))
lims    <- range(c(df_paired$Delta_obs, df_paired$Delta_sim), na.rm = TRUE)

df_stats <- df_metrics %>%
  mutate(label = sprintf("R²=%.2f  RMSE=%.2f°C\nbias=%+.2f°C", r2, rmse, bias),
         x_pos = lims[1], y_pos = lims[2])

p_scatter <- ggplot(df_paired, aes(x = Delta_obs, y = Delta_sim, colour = version)) +
  geom_hex(aes(fill = after_stat(log10(count))), colour = NA, bins = 40) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, aes(colour = version)) +
  geom_label(data = df_stats,
             aes(x = x_pos, y = y_pos, label = label, colour = version),
             fill = "white", label.size = 0.3, size = 3.2,
             hjust = -0.05, vjust = 1.1, fontface = "bold",
             inherit.aes = FALSE) +
  facet_wrap(~ version, nrow = 1) +
  scale_fill_viridis_c(name = "log₁₀(n)") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_fixed(xlim = lims, ylim = lims) +
  labs(
    title    = "MuSICA version comparison — ΔTmax at HOBO locations",
    subtitle = sprintf("Scénario : %s | n HOBO = %d | %d jours",
                       ref_sc$name, n_distinct(df_paired$id_plot),
                       n_distinct(df_paired$date)),
    x = expression("HOBO obs" ~ Delta * T[max] ~ (degree*C)),
    y = expression("MuSICA sim" ~ Delta * T[max] ~ (degree*C))
  ) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12))

save_plot(p_scatter, file.path(FIG_DIR, "comp_01_scatter.png"), width = 12, height = 6)

# ---- Figure 2: Metrics comparison -------------------------------------------

df_m_long <- df_metrics %>%
  dplyr::select(version, r2, rmse, bias) %>%
  pivot_longer(-version, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("rmse", "bias", "r2"),
                          labels = c("RMSE (°C)", "Bias (°C)", "R²")))

p_metrics <- ggplot(df_m_long, aes(x = version, y = value, fill = version)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = round(value, 3)),
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(title = "Metrics per model version", x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

save_plot(p_metrics, file.path(FIG_DIR, "comp_02_metrics.png"), width = 10, height = 5)

# ---- Figure 5: Per-plot slope Tmax_micro ~ Tmax_macro ------------------------
# Sense : slope=1 ⇒ pas de tampon | slope<1 ⇒ canopée dampe le signal macro.
# C'est la métrique standard de buffering (De Frenne et al.). Calculée par
# régression linéaire intra-plot sur la saison estivale.

df_slope <- df_paired %>%
  group_by(id_plot, version) %>%
  summarise(
    slope_sim = tryCatch(coef(lm(Tmax_sim ~ Tmax_macro))[2], error = function(e) NA_real_),
    slope_obs = tryCatch(coef(lm(Tmax_obs ~ Tmax_macro))[2], error = function(e) NA_real_),
    n_days    = n(),
    .groups   = "drop"
  ) %>% dplyr::filter(!is.na(slope_sim), !is.na(slope_obs), n_days >= 30)

df_slope_long <- bind_rows(
  df_slope %>% transmute(id_plot, source = "HOBO obs",            slope = slope_obs) %>% distinct(),
  df_slope %>% transmute(id_plot, source = as.character(version), slope = slope_sim)
) %>% mutate(source = factor(source, levels = c("HOBO obs", levels(df_paired$version))))

# Stats console
cat("\n--- Slope summary (Tmax_micro ~ Tmax_macro per plot) ---\n")
print(df_slope_long %>%
        group_by(source) %>%
        summarise(n = n(),
                  mean_slope   = round(mean(slope, na.rm = TRUE),   3),
                  median_slope = round(median(slope, na.rm = TRUE), 3),
                  sd_slope     = round(sd(slope,   na.rm = TRUE),   3),
                  .groups = "drop"))
write.csv(df_slope, file.path(OUT_ROOT, "version_slopes.csv"), row.names = FALSE)

pal_slope <- c("HOBO obs" = "#111111", pal)

p_slope_box <- ggplot(df_slope_long,
                      aes(x = source, y = slope, fill = source)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_boxplot(width = 0.55, alpha = 0.82, outlier.size = 1.1,
               outlier.alpha = 0.5) +
  geom_jitter(width = 0.13, size = 1.2, alpha = 0.32, colour = "grey25") +
  scale_fill_manual(values = pal_slope, guide = "none") +
  labs(title    = "A. Slope Tmax_micro ~ Tmax_macro — distribution per plot",
       subtitle = "Slope = 1 ⇒ no buffering | Slope < 1 ⇒ canopy dampens macro signal",
       x = NULL, y = "Regression slope (°C / °C)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

lims_sl <- range(c(df_slope$slope_obs, df_slope$slope_sim), na.rm = TRUE)
p_slope_sc <- ggplot(df_slope, aes(x = slope_obs, y = slope_sim, colour = version)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 2.2) +
  facet_wrap(~ version, nrow = 1) +
  scale_colour_manual(values = pal, guide = "none") +
  coord_fixed(xlim = lims_sl, ylim = lims_sl) +
  labs(title = "B. Per-plot slope: simulated vs HOBO observed",
       x = "HOBO observed slope", y = "MuSICA simulated slope") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

p_slope <- p_slope_box / p_slope_sc + patchwork::plot_layout(heights = c(1, 1.15))
save_plot(p_slope, file.path(FIG_DIR, "comp_05_slope.png"), width = 12, height = 11)

# ---- Figure 3: Timeseries for 3 representative HOBOs -------------------------
# Colours chosen for max discrimination:
#   HOBO obs  → near-black (#111111), solid thick   — ground truth
#   ERA5 macro → medium grey (#999999), longdash thin — reference climate
#   v_old     → vivid crimson (#E63946), solid medium
#   v3.2.3    → bright teal  (#2a9d8f), solid medium
#   → crimson vs teal = maximum perceptual distance; both distinct from black/grey

rep_hobos  <- pick_representative_hobos(df_hobo_daily)
rep_ids    <- rep_hobos$id_plot

df_ts_obs <- df_hobo_daily %>%
  filter(id_plot %in% rep_ids) %>%
  transmute(id_plot, date, value = Tmax_obs, source = "HOBO obs")

df_ts_era <- tidyr::expand_grid(id_plot = rep_ids, df_macro) %>%
  transmute(id_plot, date, value = Tmax_macro, source = "ERA5 macro")

df_ts_sims <- df_paired %>%
  filter(id_plot %in% rep_ids) %>%
  transmute(id_plot, date, value = Tmax_sim,
            source = paste0("MuSICA ", as.character(version)))

df_ts_all <- bind_rows(df_ts_obs, df_ts_era, df_ts_sims) %>%
  inner_join(rep_hobos %>% dplyr::select(id_plot, facet_label), by = "id_plot") %>%
  mutate(
    facet_label = factor(facet_label, levels = rep_hobos$facet_label),
    source      = factor(source,
                         levels = c("HOBO obs", "ERA5 macro",
                                    "MuSICA v_old", "MuSICA v3.2.3"))
  )

src_pal  <- c("HOBO obs"      = "#111111",
              "ERA5 macro"    = "#999999",
              "MuSICA v_old"  = "#E63946",
              "MuSICA v3.2.3" = "#2a9d8f")
src_lty  <- c("HOBO obs"      = "solid",
              "ERA5 macro"    = "longdash",
              "MuSICA v_old"  = "solid",
              "MuSICA v3.2.3" = "solid")
src_lwd  <- c("HOBO obs"      = 1.2,
              "ERA5 macro"    = 0.6,
              "MuSICA v_old"  = 0.9,
              "MuSICA v3.2.3" = 0.9)

p_ts <- ggplot(df_ts_all,
               aes(x = date, y = value,
                   colour   = source,
                   linetype = source,
                   linewidth = source)) +
  geom_line(alpha = 0.95) +
  facet_wrap(~ facet_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = src_pal, name = NULL) +
  scale_linetype_manual(values = src_lty, name = NULL) +
  scale_linewidth_manual(values = src_lwd, name = NULL) +
  guides(colour   = guide_legend(override.aes = list(linewidth = 1.2)),
         linetype = guide_legend(),
         linewidth = "none") +
  labs(
    title    = "Séries temporelles Tmax — 3 HOBOs représentatifs",
    subtitle = sprintf("Scénario : %s | Noir = HOBO obs | Gris long-tiret = ERA5 | Cramoisi = v_old | Vert-bleu = v3.2.3",
                       ref_sc$name),
    x = NULL,
    y = expression(T[max] ~ (degree*C))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 11),
        legend.key.width = unit(1.8, "cm"),
        strip.text       = element_text(size = 10, lineheight = 1.2))

save_plot(p_ts, file.path(FIG_DIR, "comp_03_timeseries.png"),
          width = 14, height = 11)

# ---- Figure 4: Scatter Tmax absolu obs vs sim --------------------------------
# Complémente la fig 1 (ΔTmax) en montrant les températures absolues.
# Révèle si le biais (+2°C) est constant sur toute la gamme de Tmax ou
# s'il amplifie/se réduit lors des extrêmes (canicule vs journées fraîches).

df_tmax_abs <- df_paired %>%
  mutate(Tmax_obs_abs = Tmax_obs)

lims_abs <- range(c(df_tmax_abs$Tmax_obs_abs, df_tmax_abs$Tmax_sim), na.rm = TRUE)

df_stats_abs <- df_tmax_abs %>%
  group_by(version) %>%
  summarise(
    r2   = cor(Tmax_obs_abs, Tmax_sim, use = "complete.obs")^2,
    rmse = sqrt(mean((Tmax_obs_abs - Tmax_sim)^2, na.rm = TRUE)),
    bias = mean(Tmax_sim - Tmax_obs_abs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("R²=%.2f  RMSE=%.2f°C\nbias=%+.2f°C", r2, rmse, bias))

p_scatter_abs <- ggplot(df_tmax_abs,
                        aes(x = Tmax_obs_abs, y = Tmax_sim)) +
  geom_hex(aes(fill = after_stat(log10(count))), colour = NA, bins = 45) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "black", linewidth = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9,
              colour = "white", fill = "white", alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9,
              aes(colour = version)) +
  geom_label(data = df_stats_abs,
             aes(label = label, colour = version),
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.1,
             fill = "white", linewidth = 0.3,
             size = 3.2, fontface = "bold",
             inherit.aes = FALSE) +
  facet_wrap(~ version, nrow = 1) +
  scale_fill_viridis_c(name = "log₁₀(n)", option = "magma") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_fixed(xlim = lims_abs, ylim = lims_abs) +
  labs(
    title    = "MuSICA version comparison — Tmax absolu aux capteurs HOBO",
    subtitle = sprintf("Scénario : %s | n HOBO = %d | %d jours | bissectrice = accord parfait",
                       ref_sc$name,
                       n_distinct(df_tmax_abs$id_plot),
                       n_distinct(df_tmax_abs$date)),
    x = expression("HOBO obs" ~ T[max] ~ (degree*C)),
    y = expression("MuSICA sim" ~ T[max] ~ (degree*C))
  ) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 12))

save_plot(p_scatter_abs, file.path(FIG_DIR, "comp_04_scatter_tmax.png"),
          width = 12, height = 6)

cat(sprintf("\nDone. Figures:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(FIG_DIR, "comp_01_scatter.png"),
            file.path(FIG_DIR, "comp_02_metrics.png"),
            file.path(FIG_DIR, "comp_03_timeseries.png"),
            file.path(FIG_DIR, "comp_04_scatter_tmax.png")))
