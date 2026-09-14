# ==============================================================================
# 4 forward selection figures (Tmax, V10) :
#   A. 5 × 5 grid — per-cluster GREEDY (each row recomputes order) + 1 row "All"
#   B. 5 × 5 grid — per-cluster NON-GREEDY fixed order (LAD→Hmax→fCover→LAI) + "All"
#   C. 1 × 5     — per-plot mean, GREEDY (LAI→fCover→Hmax→LAD), points coloured by cluster
#   D. 1 × 5     — per-plot mean, NON-GREEDY (LAD→Hmax→fCover→LAI), points coloured by cluster
#
# Daily × cluster aggregation for A/B ; per-plot mean for C/D.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools); library(sf)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---- HOBO cluster assignment -----------------------------------------------
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[
  , .(x, y, Cluster = as.character(Cluster))]
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
hobo_xy <- sf::st_coordinates(hobo_pts)
df_hobo_xy <- data.table(id_plot = hobo_pts$id_plot,
                          x = hobo_xy[, "X"], y = hobo_xy[, "Y"])
df_hobo_clu <- df_hobo_xy[, .(id_plot, Cluster = sapply(seq_len(.N), function(i) {
  d2 <- (df_forest$x - x[i])^2 + (df_forest$y - y[i])^2
  df_forest$Cluster[which.min(d2)]
}))]
df_hobo_clu[, Cluster := relabel_cluster(Cluster)]

# ---- Daily obs and sim -----------------------------------------------------
hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                           df_macro, CFG$ids_to_remove))[
  , .(id_plot, date, Delta_obs)]

collect_step_daily <- function(bit) {
  v10 <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                          z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- as.data.table(res)[, .(id_plot = id, date,
                                          Delta_sim = Delta_Tmax)]
  }
  rbindlist(rows)
}

# Pre-compute 16 coalitions (daily×cluster + per-plot mean variants)
cli_h1("Pre-computing all 16 coalitions")
COAL_dc <- list()   # daily × cluster aggregation
COAL_pp <- list()   # per-plot mean
for (bit in unname(.ARCH_COAL_MAP)) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) next
  m <- merge(s, hobo_obs, by = c("id_plot", "date"))
  mc <- merge(m, df_hobo_clu, by = "id_plot")
  COAL_dc[[bit]] <- mc[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
                            Delta_obs = mean(Delta_obs, na.rm = TRUE)),
                       by = .(Cluster, date)]
  COAL_pp[[bit]] <- mc[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
                            Delta_obs = mean(Delta_obs, na.rm = TRUE)),
                       by = .(id_plot, Cluster)]
}
cli_alert("16 coalitions pre-computed")

# ---- Helpers : fit + greedy path ------------------------------------------
fit_overall <- function(d) {
  ok <- !is.na(d$Delta_sim) & !is.na(d$Delta_obs)
  x <- d$Delta_obs[ok]; y <- d$Delta_sim[ok]
  if (length(x) < 5) return(list(r = NA, RMSE = NA, MAE = NA, n = length(x)))
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       n = length(x))
}
fit_subset <- function(d, cluster) {
  if (cluster == "All") return(fit_overall(d))
  fit_overall(d[Cluster == cluster])
}

VARS <- c("LAI", "Hmax", "fCover", "LAD")
bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")

greedy_path_subset <- function(coal_list, cluster) {
  active <- character(0); path <- list()
  b0 <- bit_from_set(active); f0 <- fit_subset(coal_list[[b0]], cluster)
  path[[1]] <- data.table(step = 1, Cluster = cluster, added = "Baseline",
                          bit = b0, r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE, n = f0$n)
  for (k in 2:5) {
    remaining <- setdiff(VARS, active)
    cand_bits <- sapply(remaining, function(v) bit_from_set(c(active, v)))
    cand_r <- sapply(cand_bits, function(b) fit_subset(coal_list[[b]], cluster)$r)
    if (all(is.na(cand_r))) break
    best_i <- which.max(cand_r)
    best_v <- remaining[best_i]
    active <- c(active, best_v); bit <- cand_bits[best_i]
    f <- fit_subset(coal_list[[bit]], cluster)
    path[[k]] <- data.table(step = k, Cluster = cluster, added = best_v,
                            bit = bit, r = f$r, RMSE = f$RMSE, MAE = f$MAE, n = f$n)
    if (length(active) == length(VARS)) break
  }
  rbindlist(path)
}

# Fixed (non-greedy) order : LAD → Hmax → fCover → LAI (inverse forward)
FIXED_ORDER <- c("LAD", "Hmax", "fCover", "LAI")
fixed_path_subset <- function(coal_list, cluster) {
  active <- character(0); path <- list()
  b0 <- bit_from_set(active); f0 <- fit_subset(coal_list[[b0]], cluster)
  path[[1]] <- data.table(step = 1, Cluster = cluster, added = "Baseline",
                          bit = b0, r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE, n = f0$n)
  for (k in 2:5) {
    active <- c(active, FIXED_ORDER[k-1]); bit <- bit_from_set(active)
    f <- fit_subset(coal_list[[bit]], cluster)
    path[[k]] <- data.table(step = k, Cluster = cluster,
                            added = FIXED_ORDER[k-1],
                            bit = bit, r = f$r, RMSE = f$RMSE, MAE = f$MAE, n = f$n)
  }
  rbindlist(path)
}

clusters <- c("P1", "P2", "P3", "P4", "All")

# ============================================================================
# Build common figure helper
# ============================================================================
make_grid_fig <- function(paths_list, coal_list, png_path, title_str) {
  # paths_list : list of data.tables per cluster (incl. "All")
  long <- list()
  for (cl in names(paths_list)) {
    p <- paths_list[[cl]]
    for (k in seq_len(nrow(p))) {
      bit <- p$bit[k]
      d <- coal_list[[bit]]
      if (is.null(d)) next
      if (cl == "All") {
        d_panel <- copy(d)
      } else {
        d_panel <- d[Cluster == cl]
      }
      d_panel[, Cluster_row := cl]
      d_panel[, step := k]
      d_panel[, added := p$added[k]]
      long[[length(long) + 1L]] <- d_panel
    }
  }
  long <- rbindlist(long, fill = TRUE)
  long[, Cluster_row := factor(Cluster_row, levels = clusters)]

  # Per-row labels (each row has its own variable ordering)
  ann_long <- list()
  for (cl in names(paths_list)) {
    p <- paths_list[[cl]]
    for (k in seq_len(nrow(p))) {
      lab <- if (k == 1) "Baseline"
              else if (k == nrow(p)) sprintf("REF (+%s)",
                gsub("Hmax", "H_max", p$added[k]))
              else paste("+", paste(gsub("Hmax", "H_max", p$added[2:k]),
                                      collapse = " + "))
      ann_long[[length(ann_long) + 1L]] <- data.table(
        Cluster_row = cl, step = k, label = lab,
        lab_r = sprintf("italic(r)==%+.2f", p$r[k]),
        lab_rmse = sprintf("RMSE==%.2f * ' °C'", p$RMSE[k]))
    }
  }
  ann <- rbindlist(ann_long)
  ann[, Cluster_row := factor(Cluster_row, levels = clusters)]
  long <- merge(long, ann[, .(Cluster_row, step, label)],
                  by = c("Cluster_row", "step"))

  xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

  # Add column label combining step + variable
  ann[, col_lab := sprintf("Step %d", step)]
  long[, col_lab := sprintf("Step %d", step)]

  # Slide-quality grid : large fonts, large points, 600 dpi PNG + PDF
  p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim, colour = Cluster)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_point(alpha = 0.75, size = 2.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey25", linewidth = 1.0, aes(group = 1)) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
              parse = TRUE, hjust = -0.08, vjust = 1.4,
              size = 5, fontface = "bold",
              inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
              parse = TRUE, hjust = -0.08, vjust = 2.8,
              size = 4.5, fontface = "bold",
              inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = Inf, y = -Inf, label = label),
              hjust = 1.05, vjust = -0.5,
              size = 4.2, fontface = "italic",
              inherit.aes = FALSE, colour = "grey25") +
    facet_grid(Cluster_row ~ step) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                    P3 = "#7570B3", P4 = "#E7298A"),
                          na.translate = FALSE) +
    labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
         y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
    theme_bw(base_size = 22) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text  = element_text(face = "bold", size = 20),
          axis.title  = element_text(face = "bold", size = 22),
          axis.text   = element_text(size = 14),
          legend.text  = element_text(size = 18),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          panel.grid.minor = element_blank(),
          # 16:9 figure : 5 cols × 5 rows. Total plot area ~ 18 x 9 in (legend
          # + axes ~ 2 in vertical). Panel aspect = 18/5 : 9/5 = 3.6/1.8 = 2.
          # We keep aspect.ratio slightly < 1 so each panel is wider than tall,
          # consistent with the 16:9 outer ratio.
          aspect.ratio = 0.6)
  # Slide format 16:9 : 24 x 13.5 in @ 600 dpi PNG + PDF vector
  ggsave(png_path, p, width = 24, height = 13.5, dpi = 600, bg = "white")
  ggsave(sub("\\.png$", ".pdf", png_path), p, width = 24, height = 13.5,
          device = cairo_pdf)
  cli_alert_success("Saved {basename(png_path)} (+ .pdf, slide 16:9)")
  invisible(list(long = long, ann = ann))
}

# ============================================================================
# A. Greedy per cluster + "All" row
# ============================================================================
cli_h1("A. Greedy per cluster + All (daily × cluster)")
paths_g <- lapply(clusters, function(cl) greedy_path_subset(COAL_dc, cl))
names(paths_g) <- clusters
make_grid_fig(paths_g, COAL_dc,
              file.path(OUT, "fig_forward_HOBO_v10_greedy_per_cluster_plus_all.png"),
              "V10 — Greedy forward selection PER CLUSTER + All (daily × cluster)")

# ============================================================================
# B. Non-greedy fixed order (LAD→Hmax→fCover→LAI) per cluster + All
# ============================================================================
cli_h1("B. Non-greedy fixed (LAD→Hmax→fCover→LAI) per cluster + All")
paths_f <- lapply(clusters, function(cl) fixed_path_subset(COAL_dc, cl))
names(paths_f) <- clusters
make_grid_fig(paths_f, COAL_dc,
              file.path(OUT, "fig_forward_HOBO_v10_fixed_per_cluster_plus_all.png"),
              "V10 — Fixed forward order LAD→Hmax→fCover→LAI, per cluster + All (daily × cluster)")

# ============================================================================
# C. Per-plot mean — greedy
# ============================================================================
cli_h1("C. Per-plot mean, greedy")
# Run greedy on per-plot mean using only one dataset (53 plots)
fit_pp <- function(d) {
  ok <- !is.na(d$Delta_sim) & !is.na(d$Delta_obs)
  x <- d$Delta_obs[ok]; y <- d$Delta_sim[ok]
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       n = length(x))
}
active <- character(0); path_pp_g <- list()
b0 <- bit_from_set(active); f0 <- fit_pp(COAL_pp[[b0]])
path_pp_g[[1]] <- data.table(step = 1, added = "Baseline", bit = b0,
                              r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE, n = f0$n)
for (k in 2:5) {
  remaining <- setdiff(VARS, active)
  cand_bits <- sapply(remaining, function(v) bit_from_set(c(active, v)))
  cand_r <- sapply(cand_bits, function(b) fit_pp(COAL_pp[[b]])$r)
  best_i <- which.max(cand_r); best_v <- remaining[best_i]
  active <- c(active, best_v); bit <- cand_bits[best_i]
  f <- fit_pp(COAL_pp[[bit]])
  path_pp_g[[k]] <- data.table(step = k, added = best_v, bit = bit,
                                r = f$r, RMSE = f$RMSE, MAE = f$MAE, n = f$n)
}
greedy_pp <- rbindlist(path_pp_g)
print(greedy_pp)

make_perplot_fig <- function(path, png_path, title_str) {
  labels_path <- sapply(seq_len(nrow(path)), function(k) {
    if (k == 1) "Baseline"
    else if (k == nrow(path)) sprintf("REF (+%s)",
       gsub("Hmax", "H_max", path$added[k]))
    else paste("+", paste(gsub("Hmax", "H_max", path$added[2:k]),
                              collapse = " + "))
  })
  long <- list()
  for (k in seq_len(nrow(path))) {
    d <- copy(COAL_pp[[ path$bit[k] ]])
    d[, step := k]; d[, label := labels_path[k]]
    long[[k]] <- d
  }
  long <- rbindlist(long, fill = TRUE)
  long$label <- factor(as.character(long$label), levels = labels_path)
  ann <- path[, .(label = labels_path,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                   lab_mae  = sprintf("MAE==%.2f * ' °C'", MAE))]
  ann$label <- factor(ann$label, levels = labels_path)
  xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.06; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

  p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim, colour = Cluster, shape = Cluster)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_point(alpha = 0.85, size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.9, alpha = 0.22, aes(group = 1)) +
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
    coord_fixed(xlim = xy_lim, ylim = xy_lim) +
    scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                    P3 = "#7570B3", P4 = "#E7298A")) +
    scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)) +
    labs(title = title_str,
         x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
         y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
  ggsave(png_path, p, width = 18, height = 6, dpi = 300)
  cli_alert_success("Saved {basename(png_path)}")
}

make_perplot_fig(greedy_pp,
                  file.path(OUT, "fig_forward_HOBO_v10_perplot_greedy_clustercolour.png"),
                  "V10 — Per-plot mean, GREEDY forward (coloured by cluster)")
fwrite(greedy_pp, file.path(OUT, "tab_v10_HOBO_perplot_greedy_path.csv"))

# ============================================================================
# D. Per-plot mean — non-greedy fixed (LAD→Hmax→fCover→LAI)
# ============================================================================
cli_h1("D. Per-plot mean, non-greedy fixed (LAD→Hmax→fCover→LAI)")
active <- character(0); path_pp_f <- list()
b0 <- bit_from_set(active); f0 <- fit_pp(COAL_pp[[b0]])
path_pp_f[[1]] <- data.table(step = 1, added = "Baseline", bit = b0,
                              r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE, n = f0$n)
for (k in 2:5) {
  active <- c(active, FIXED_ORDER[k-1]); bit <- bit_from_set(active)
  f <- fit_pp(COAL_pp[[bit]])
  path_pp_f[[k]] <- data.table(step = k, added = FIXED_ORDER[k-1], bit = bit,
                                r = f$r, RMSE = f$RMSE, MAE = f$MAE, n = f$n)
}
fixed_pp <- rbindlist(path_pp_f)
print(fixed_pp)
make_perplot_fig(fixed_pp,
                  file.path(OUT, "fig_forward_HOBO_v10_perplot_fixed_clustercolour.png"),
                  "V10 — Per-plot mean, FIXED order LAD→Hmax→fCover→LAI (coloured by cluster)")
fwrite(fixed_pp, file.path(OUT, "tab_v10_HOBO_perplot_fixed_path.csv"))

cli_h1("Done — 4 figures generated")
