# ==============================================================================
# V10 (fCover_baseline = mean = 0.869) — 4 figures slope vs Tmax :
#   1. fig_heatmap_LOO_v10_slope.png       LOO Δ on hourly slope (per archetype)
#   2. fig_heatmap_LOO_v10_Tmax.png        LOO Δ on daily Tmax mean (°C)
#   3. fig_forward_HOBO_v10_slope.png      Inverse forward, slope-based
#   4. fig_forward_HOBO_v10_Tmax.png       Inverse forward, ΔTmax-based
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

# ============================================================================
# Helpers — archetype Tmax and slope, V10 fallback to V9 for fCover=1 bits
# ============================================================================
sc_name_of <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

arch_nc_path <- function(cl, bit) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  v10_path <- file.path("out_files/H1_archetypes_v10_fcovmean",
                        arch_lbl, sprintf("musica_out_%s_%s.nc",
                                          arch_lbl, sc_name_of[bit]))
  v9_path  <- file.path("out_files/H1_archetypes_v9",
                        arch_lbl, sprintf("musica_out_%s_%s.nc",
                                          arch_lbl, sc_name_of[bit]))
  if (file.exists(v10_path)) v10_path else v9_path
}

extract_metric <- function(path, metric = c("Tmax","slope")) {
  metric <- match.arg(metric)
  if (!file.exists(path)) return(NA_real_)
  if (metric == "Tmax") {
    res <- tryCatch(extract_deltatmax_one(path, df_macro, CFG$date_seq,
                                          z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) return(NA_real_)
    mean(res$Tmax_micro, na.rm = TRUE)
  } else {
    res <- tryCatch(extract_hourly_slope_one(path, era5_hourly, CFG$date_seq,
                                             z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) return(NA_real_)
    res$slope
  }
}

# Build wide table : Cluster × bit, for one metric
build_wide_v10 <- function(metric) {
  rows <- list()
  for (cl in 1:4) for (bit in unname(.ARCH_COAL_MAP)) {
    rows[[length(rows) + 1L]] <- data.table(
      Cluster = cl, bit = bit,
      val = extract_metric(arch_nc_path(cl, bit), metric))
  }
  dcast(rbindlist(rows), Cluster ~ bit, value.var = "val")
}

VARS <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

build_heatmap <- function(W, title_str, fill_lab, value_fmt, png_path,
                           row_order = c("LAD", "Hmax", "fCover", "LAI")) {
  # Row labels with italic H_max + subscripts, plotmath-friendly, AND units.
  # variable -> plotmath expression string
  var_lab_pm <- c(
    LAI    = "atop(LAI, \"(m\"^\"2\"~\"m\"^\"-2\"*\")\")",
    Hmax   = "atop(italic(H)[max], \"(m)\")",
    fCover = "atop(fCover, \"(\\u2013)\")",
    LAD    = "atop(LAD, \"(m\"^\"2\"~\"m\"^\"-3\"*\")\")"
  )
  delta_long <- list()
  for (v in VARS) for (i in seq_len(nrow(W))) {
    delta_long[[length(delta_long) + 1L]] <- data.table(
      Cluster = W$Cluster[i], variable = v,
      delta   = W[["1111"]][i] - W[[loo_bits[v]]][i])
  }
  long <- rbindlist(delta_long)
  # Relabel cluster IDs by ASCENDING LAI (CLUSTER_RELABEL helper from
  # R/cluster_relabel.R) : old C3 → P1, C2 → P2, C4 → P3, C1 → P4.
  long[, Profile := as.character(relabel_cluster(Cluster))]
  avg_rows <- long[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                       by = variable]
  long_full <- rbind(long[, .(variable, Profile, delta)], avg_rows)
  long_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
  # Row order : weakest → strongest (LAD top → LAI bottom).
  # Coherent avec le forward inverse Baseline → +LAD → +Hmax → +fCover → REF +LAI.
  long_full[, var_lab := var_lab_pm[as.character(variable)]]
  var_lab_levels <- var_lab_pm[row_order]
  long_full[, var_lab := factor(var_lab, levels = var_lab_levels)]

  max_abs <- max(abs(long_full$delta), na.rm = TRUE) * 1.05

  # Axes inverses : variables en X (LAD → LAI gauche-droite),
  # profils P1..P4 + Avg en Y (P1 sparse en haut, P4 dense en bas, Avg séparée).
  p <- ggplot(long_full, aes(x = var_lab, y = Profile, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf(value_fmt, delta)),
              fontface = "bold", size = 8,
              colour = ifelse(abs(long_full$delta) > max_abs * 0.55, "white", "grey15")) +
    # Séparateur horizontal entre les 4 profils et la ligne Avg
    geom_hline(yintercept = 1.5, colour = "white", linewidth = 4) +
    geom_hline(yintercept = 1.5, colour = "grey40", linewidth = 1.0) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = c(-max_abs, max_abs),
                          name = fill_lab) +
    scale_x_discrete(labels = function(x) parse(text = as.character(x))) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 24) +
    theme(axis.text.x = element_text(face = "bold", size = 24),
          axis.text.y = element_text(face = "bold", size = 24),
          axis.ticks  = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22, face = "bold"),
          legend.key.height = unit(1.6, "cm"),
          legend.position = "right")
  # Slide 16:9 — 16x9 in
  ggsave(png_path, p, width = 16, height = 9, dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", png_path), p, width = 16, height = 9,
          device = cairo_pdf)
  cli_alert_success("Saved {basename(png_path)} (+ .pdf)")
}

# ============================================================================
# 1 & 2 — Heatmaps V10 (slope and Tmax)
# ============================================================================
cli_h1("Heatmaps V10")

W_tmax  <- build_wide_v10("Tmax")
W_slope <- build_wide_v10("slope")

build_heatmap(W_tmax,
              "V10 — LOO Δ on daily Tmax (°C)",
              expression(Delta[v]^"LOO" ~ "(°C)"),
              "%+.2f",
              file.path(OUT, "fig_heatmap_LOO_v10_Tmax.png"))

build_heatmap(W_slope,
              "V10 — LOO Δ on hourly slope",
              expression(Delta[v]^"LOO" ~ "(slope)"),
              "%+.3f",
              file.path(OUT, "fig_heatmap_LOO_v10_slope.png"))

# ============================================================================
# 3 & 4 — Forward HOBO V10 (slope and Tmax)
# ============================================================================
cli_h1("Forward HOBO V10 (inverse path : weakest → strongest)")

collect_step_v10 <- function(bit, metric = c("slope", "Tmax")) {
  metric <- match.arg(metric)
  v10_dir <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9_dir  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10_dir) &&
                 length(list.files(v10_dir, "\\.nc$")) > 0) v10_dir else v9_dir
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    if (metric == "Tmax") {
      res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                      error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id,
                                sim_metric = mean(res$Delta_Tmax, na.rm = TRUE))
    } else {
      res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                                z_target = CFG$tair_target_height),
                      error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id, sim_metric = res$slope)
    }
  }
  rbindlist(rows)
}

# Observed
hobo_obs_Tmax <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                                df_macro, CFG$ids_to_remove))[
  , .(obs_metric = mean(Delta_obs, na.rm = TRUE)), by = id_plot]

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
hobo_obs_slope <- hobo_with_era5[, .(obs_metric = as.numeric(coef(lm(t_hobo ~ Tair_era5))[2])),
                                  by = id_plot]

build_forward <- function(metric, png_path, title_str, x_lab, y_lab, fmt_axis) {
  inv <- if (metric == "Tmax") list(
    bits   = c("0000","0001","0101","0111","1111"),
    labels = c("Baseline","+LAD","+LAD+H_max","+LAD+H_max+fCover","REF (+LAI)")
  ) else list(
    bits   = c("0000","0001","0101","1101","1111"),
    labels = c("Baseline","+LAD","+LAD+H_max","+LAD+H_max+LAI","REF (+fCover)")
  )
  obs <- if (metric == "Tmax") hobo_obs_Tmax else hobo_obs_slope

  long <- list(); fits <- list()
  for (k in seq_along(inv$bits)) {
    bit <- inv$bits[k]
    d <- collect_step_v10(bit, metric)
    if (nrow(d) == 0) next
    d <- merge(d, obs, by = "id_plot")
    d[, step := k]; d[, label := inv$labels[k]]
    long[[k]] <- d
    ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
    x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
    fits[[k]] <- data.table(step = k, label = inv$labels[k],
                             r = suppressWarnings(cor(x, y)),
                             RMSE = sqrt(mean((y - x)^2)),
                             MAE  = mean(abs(y - x)),
                             bias = mean(y - x), n = length(x))
  }
  long <- rbindlist(long, fill = TRUE)
  fits <- rbindlist(fits)
  long$label <- factor(as.character(long$label), levels = inv$labels)

  ann <- fits[, .(label,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf(fmt_axis$lab_rmse, RMSE),
                   lab_mae  = sprintf(fmt_axis$lab_mae,  MAE))]
  ann$label <- factor(ann$label, levels = inv$labels)
  xy_lim <- range(c(long$sim_metric, long$obs_metric), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.06; xy_lim <- xy_lim + c(-xy_pad, xy_pad)
  ref_line <- if (metric == "slope") 1 else 0

  p <- ggplot(long, aes(x = obs_metric, y = sim_metric)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = ref_line, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = ref_line, linetype = "dotted", colour = "grey80") +
    geom_point(colour = "#2C5F2D", alpha = 0.75, size = 2.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400",
                linewidth = 0.9, alpha = 0.22) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
              parse = TRUE, hjust = -0.08, vjust = 1.4, size = 5.5,
              inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
              parse = TRUE, hjust = -0.08, vjust = 2.8, size = 5.5,
              inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
              parse = TRUE, hjust = -0.08, vjust = 4.2, size = 5.5,
              inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ label, nrow = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    labs(title = title_str, x = x_lab, y = y_lab) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 14),
          panel.grid.minor = element_blank())
  ggsave(png_path, p, width = 18, height = 5, dpi = 300)
  cli_alert_success("Saved {basename(png_path)}")
  cli_alert("Fits :"); print(fits)
}

build_forward("Tmax",
               file.path(OUT, "fig_forward_HOBO_v10_Tmax.png"),
               "V10 — Inverse forward selection on ΔTmax (°C)",
               expression("Field " ~ Delta * T[max] ~ "(°C)"),
               expression("Simulated " ~ Delta * T[max] ~ "(°C)"),
               list(lab_rmse = "RMSE==%.2f * ' °C'",
                    lab_mae  = "MAE==%.2f * ' °C'"))

build_forward("slope",
               file.path(OUT, "fig_forward_HOBO_v10_slope.png"),
               "V10 — Inverse forward selection on hourly slope",
               "Field slope (Tair_HOBO ~ Tair_ERA5)",
               "Simulated slope (Tair_sim ~ Tair_ERA5)",
               list(lab_rmse = "RMSE==%.3f",
                    lab_mae  = "MAE==%.3f"))

cli_h1("Done — 4 V10 standalone figures (slope vs Tmax).")
