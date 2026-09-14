# ==============================================================================
# Generate V10 figures (fCover_b = mean = 0.869) :
#   1. Heatmap LOO Δ in °C (replaces V9 — answers Eva's question)
#   2. Inverse forward HOBO (slope) — V10 baseline applied where fCover=0
# Plus a side-by-side comparison V9 vs V10 heatmap (so we can discuss with the
# committee why the conclusion changes).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools); library(patchwork)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

# Helper : Tmax mean from V10 (if NC exists) else V9
archetype_tmax <- function(cl, bit, version = c("V9", "V10")) {
  version <- match.arg(version)
  arch_lbl <- sprintf("Arch_C%s", cl)
  sc_nm    <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))[bit]
  v10_path <- file.path("out_files/H1_archetypes_v10_fcovmean",
                        arch_lbl, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
  v9_path  <- file.path("out_files/H1_archetypes_v9",
                        arch_lbl, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
  path <- if (version == "V10" && file.exists(v10_path)) v10_path else v9_path
  res <- tryCatch(extract_deltatmax_one(path, df_macro, CFG$date_seq,
                                        z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

build_wide <- function(version) {
  rows <- list()
  for (cl in 1:4) for (bit in unname(.ARCH_COAL_MAP)) {
    rows[[length(rows) + 1L]] <- data.table(
      Cluster = cl, bit = bit,
      Tmax = archetype_tmax(cl, bit, version = version))
  }
  dcast(rbindlist(rows), Cluster ~ bit, value.var = "Tmax")
}

VARS <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

loo_delta_long <- function(W, version_lbl) {
  delta_long <- list()
  for (v in VARS) for (i in seq_len(nrow(W))) {
    delta_long[[length(delta_long) + 1L]] <- data.table(
      version = version_lbl,
      Cluster = W$Cluster[i], variable = v,
      delta   = W[["1111"]][i] - W[[loo_bits[v]]][i])
  }
  rbindlist(delta_long)
}

W_V9  <- build_wide("V9")
W_V10 <- build_wide("V10")
DL    <- rbind(loo_delta_long(W_V9, "V9 (fcov_b=1)"),
                loo_delta_long(W_V10, "V10 (fcov_b=mean)"))
# Cluster = RAW k-means code; map to display label via cluster_relabel
# (raw 1->P4 densest, 2->P2, 3->P1 sparsest, 4->P3). Direct P%d was WRONG (scrambled).
.RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
DL[, Profile := unname(.RELAB[as.character(Cluster)])]
avg_rows <- DL[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                  by = .(version, variable)]
DL_full <- rbind(DL[, .(version, variable, Profile, delta)], avg_rows)

# Use the slope-based forward order for row ordering (consistent with V9 narrative)
forward_order <- c("LAI", "Hmax", "LAD", "fCover")
DL_full[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]
var_lab_order <- ifelse(forward_order == "Hmax", "H[max]", forward_order)
DL_full[, var_lab := factor(var_lab, levels = var_lab_order)]
DL_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
DL_full[, version := factor(version, levels = c("V9 (fcov_b=1)", "V10 (fcov_b=mean)"))]

max_abs <- max(abs(DL_full$delta), na.rm = TRUE) * 1.05

make_panel <- function(d) {
  ggplot(d, aes(x = Profile, y = var_lab, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f", delta)),
              fontface = "bold", size = 4.5,
              colour = ifelse(abs(d$delta) > max_abs * 0.55, "white", "grey15")) +
    geom_vline(xintercept = 4.5, colour = "white", linewidth = 4) +
    geom_vline(xintercept = 4.5, colour = "grey40", linewidth = 1.0) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = c(-max_abs, max_abs),
                          name = expression(Delta[v]^"LOO" ~ "(°C)")) +
    scale_y_discrete(limits = rev,
                     labels = function(x) parse(text = as.character(x))) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(face = "bold", size = 13),
          axis.text.y = element_text(face = "bold", size = 13),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 14))
}

p_v9  <- make_panel(DL_full[version == "V9 (fcov_b=1)"]) +
  labs(title = "V9 : fCover_baseline = 1 (closed canopy)")
p_v10 <- make_panel(DL_full[version == "V10 (fcov_b=mean)"]) +
  labs(title = "V10 : fCover_baseline = mean(fCover) = 0.869")

# Single V10 heatmap
ggsave(file.path(OUT, "fig_heatmap_LOO_v10_fcovmean.png"),
        p_v10 + theme(legend.position = "right"),
        width = 12, height = 5, dpi = 300)
cli_alert_success("Saved fig_heatmap_LOO_v10_fcovmean.png")

# V9 vs V10 side-by-side
combo <- (p_v9 + theme(legend.position = "none")) | p_v10
ggsave(file.path(OUT, "fig_heatmap_LOO_v9_vs_v10_fcover_baseline.png"),
        combo, width = 18, height = 5, dpi = 300)
cli_alert_success("Saved fig_heatmap_LOO_v9_vs_v10_fcover_baseline.png")

# ----------------------------------------------------------------------------
# Inverse forward HOBO with V10 baseline
# ----------------------------------------------------------------------------
cli_h1("V10 inverse forward HOBO (slope)")

collect_step_v10 <- function(bit, metric = c("slope", "Tmax")) {
  metric <- match.arg(metric)
  v10_dir <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9_dir  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10_dir) && length(list.files(v10_dir, "\\.nc$")) > 0) v10_dir else v9_dir
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    if (metric == "slope") {
      res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                                z_target = CFG$tair_target_height),
                       error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id, sim_metric = res$slope)
    } else {
      res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                             z_target = CFG$tair_target_height),
                       error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id, sim_metric = mean(res$Delta_Tmax, na.rm = TRUE))
    }
  }
  rbindlist(rows)
}

hobo_raw <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a",
         as.Date(datetime) %in% CFG$date_seq,
         !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(id_plot, time) %>%
  summarise(t_hobo = mean(t_hobo, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()
hobo_with_era5 <- merge(hobo_raw, era5_hourly, by = "time")
hobo_slope_obs <- hobo_with_era5[, {
  fit <- lm(t_hobo ~ Tair_era5)
  .(slope_obs = as.numeric(coef(fit)[2]))
}, by = id_plot]

inv_slope <- list(
  bits   = c("0000","0001","0101","1101","1111"),
  labels = c("Baseline","+LAD","+LAD+H_max","+LAD+H_max+LAI","REF (+fCover)"))

long <- list(); fits <- list()
for (k in seq_along(inv_slope$bits)) {
  bit <- inv_slope$bits[k]
  d <- collect_step_v10(bit, "slope")
  if (nrow(d) == 0) next
  d <- merge(d, hobo_slope_obs[, .(id_plot, obs_metric = slope_obs)], by = "id_plot")
  d[, step := k]; d[, label := inv_slope$labels[k]]; d[, bit := bit]
  long[[k]] <- d
  ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
  x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
  fits[[k]] <- data.table(step = k, label = inv_slope$labels[k],
                           r = suppressWarnings(cor(x, y)),
                           RMSE = sqrt(mean((y - x)^2)),
                           MAE  = mean(abs(y - x)), n = length(x))
}
long <- rbindlist(long, fill = TRUE)
fits <- rbindlist(fits)
long$label <- factor(as.character(long$label), levels = inv_slope$labels)

ann <- fits[, .(label,
                 lab_r    = sprintf("italic(r)==%+.2f", r),
                 lab_rmse = sprintf("RMSE==%.3f", RMSE),
                 lab_mae  = sprintf("MAE==%.3f", MAE))]
ann$label <- factor(ann$label, levels = inv_slope$labels)
xy_lim <- range(c(long$sim_metric, long$obs_metric), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.06; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

p_inv <- ggplot(long, aes(x = obs_metric, y = sim_metric)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey80") +
  geom_point(colour = "#2C5F2D", alpha = 0.75, size = 2.0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
            parse = TRUE, hjust = -0.08, vjust = 4.2, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ label, nrow = 1) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  labs(title = "Inverse forward (V10, fCover_baseline = 0.869) — weakest → strongest",
        x = "Field slope", y = "Simulated slope") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_forward_selection_HOBO_v10_inverse_slope.png"),
        p_inv, width = 18, height = 5, dpi = 300)
cli_alert_success("Saved fig_forward_selection_HOBO_v10_inverse_slope.png")
cli_alert("V10 inverse slope fits :"); print(fits)
fwrite(fits, file.path(OUT, "tab_v10_HOBO_inverse_slope.csv"))
