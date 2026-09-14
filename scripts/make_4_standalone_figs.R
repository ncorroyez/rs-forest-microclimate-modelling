# ==============================================================================
# 4 figures standalone (Tmax-based, slide-ready) :
#   1. fig_heatmap_LOO_Tmax_v9.png         — V9 (fCover_baseline = 1)
#   2. fig_heatmap_LOO_Tmax_v10.png        — V10 (fCover_baseline = mean = 0.869)
#   3. fig_forward_HOBO_Tmax_v9.png        — V9 inverse path
#   4. fig_forward_HOBO_Tmax_v10.png       — V10 inverse path (recommended)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ============================================================================
# Helpers
# ============================================================================
arch_tmax <- function(cl, bit, version = c("V9", "V10")) {
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
      Cluster = cl, bit = bit, Tmax = arch_tmax(cl, bit, version = version))
  }
  dcast(rbindlist(rows), Cluster ~ bit, value.var = "Tmax")
}

VARS <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

build_heatmap <- function(W, version_lbl, png_path) {
  delta_long <- list()
  for (v in VARS) for (i in seq_len(nrow(W))) {
    delta_long[[length(delta_long) + 1L]] <- data.table(
      Cluster = W$Cluster[i], variable = v,
      delta   = W[["1111"]][i] - W[[loo_bits[v]]][i])
  }
  long <- rbindlist(delta_long)
  # Cluster = RAW k-means code; map to display label via cluster_relabel
  # (raw 1->P4 densest, 2->P2, 3->P1 sparsest, 4->P3). Direct P%d was WRONG (scrambled).
  .RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
  long[, Profile := unname(.RELAB[as.character(Cluster)])]
  avg_rows <- long[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                       by = variable]
  long_full <- rbind(long[, .(variable, Profile, delta)], avg_rows)
  long_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]
  forward_order <- c("LAI", "Hmax", "LAD", "fCover")
  long_full[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]
  var_lab_order <- ifelse(forward_order == "Hmax", "H[max]", forward_order)
  long_full[, var_lab := factor(var_lab, levels = var_lab_order)]

  max_abs <- max(abs(long_full$delta), na.rm = TRUE) * 1.05

  p <- ggplot(long_full, aes(x = Profile, y = var_lab, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f", delta)),
              fontface = "bold", size = 5.5,
              colour = ifelse(abs(long_full$delta) > max_abs * 0.55, "white", "grey15")) +
    geom_vline(xintercept = 4.5, colour = "white", linewidth = 4) +
    geom_vline(xintercept = 4.5, colour = "grey40", linewidth = 1.0) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, limits = c(-max_abs, max_abs),
                          name = expression(Delta[v]^"LOO" ~ "(°C)")) +
    scale_y_discrete(limits = rev,
                     labels = function(x) parse(text = as.character(x))) +
    labs(x = NULL, y = NULL, title = version_lbl) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(face = "bold", size = 15),
          axis.text.y = element_text(face = "bold", size = 15),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 15))
  ggsave(png_path, p, width = 12, height = 5, dpi = 300)
  cli_alert_success("Saved {basename(png_path)}")
}

# ============================================================================
# Heatmaps (LOO Δ Tmax, °C) — Fig 1 & 2
# ============================================================================
cli_h1("Heatmaps")
W_V9  <- build_wide("V9")
W_V10 <- build_wide("V10")

build_heatmap(W_V9,  "Baseline fCover = 1 (closed canopy)",
              file.path(OUT, "fig_heatmap_LOO_Tmax_v9.png"))
build_heatmap(W_V10, "Baseline fCover = mean(fCover) = 0.869",
              file.path(OUT, "fig_heatmap_LOO_Tmax_v10.png"))

# ============================================================================
# Forward HOBO (Tmax) — Fig 3 & 4
# ============================================================================
cli_h1("Forward HOBO (Tmax-based, inverse path : weakest → strongest)")

collect_step <- function(bit, version = c("V9","V10")) {
  version <- match.arg(version)
  v10_dir <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9_dir  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (version == "V10" && dir.exists(v10_dir) &&
                 length(list.files(v10_dir, "\\.nc$")) > 0) v10_dir else v9_dir
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                          z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- data.table(id_plot = id,
                              sim_metric = mean(res$Delta_Tmax, na.rm = TRUE))
  }
  rbindlist(rows)
}

hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                            df_macro, CFG$ids_to_remove))[
  , .(obs_metric = mean(Delta_obs, na.rm = TRUE)), by = id_plot]

build_forward_fig <- function(version, png_path, title_str) {
  inv_tmax <- list(
    bits   = c("0000","0001","0101","0111","1111"),
    labels = c("Baseline","+LAD","+LAD+H_max","+LAD+H_max+fCover","REF (+LAI)"))

  long <- list(); fits <- list()
  for (k in seq_along(inv_tmax$bits)) {
    bit <- inv_tmax$bits[k]
    d <- collect_step(bit, version)
    if (nrow(d) == 0) next
    d <- merge(d, hobo_obs, by = "id_plot")
    d[, step := k]; d[, label := inv_tmax$labels[k]]; d[, bit := bit]
    long[[k]] <- d
    ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
    x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
    fits[[k]] <- data.table(step = k, label = inv_tmax$labels[k],
                             r = suppressWarnings(cor(x, y)),
                             RMSE = sqrt(mean((y - x)^2)),
                             MAE  = mean(abs(y - x)), bias = mean(y - x), n = length(x))
  }
  long <- rbindlist(long, fill = TRUE)
  fits <- rbindlist(fits)
  long$label <- factor(as.character(long$label), levels = inv_tmax$labels)

  ann <- fits[, .(label,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                   lab_mae  = sprintf("MAE==%.2f * ' °C'", MAE))]
  ann$label <- factor(ann$label, levels = inv_tmax$labels)
  xy_lim <- range(c(long$sim_metric, long$obs_metric), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.06; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

  p <- ggplot(long, aes(x = obs_metric, y = sim_metric)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
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
    labs(title = title_str,
         x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
         y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 14),
          panel.grid.minor = element_blank())
  ggsave(png_path, p, width = 18, height = 5, dpi = 300)
  cli_alert_success("Saved {basename(png_path)}")
  print(fits)
}

build_forward_fig("V9",
                   file.path(OUT, "fig_forward_HOBO_Tmax_v9.png"),
                   "Baseline fCover = 1 (closed canopy)")
build_forward_fig("V10",
                   file.path(OUT, "fig_forward_HOBO_Tmax_v10.png"),
                   "Baseline fCover = mean(fCover) = 0.869")

cli_h1("Done — 4 standalone figures generated in outputs/figs_MEB2026_final/")
