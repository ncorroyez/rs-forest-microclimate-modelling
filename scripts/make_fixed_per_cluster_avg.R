# ==============================================================================
# Fixed-order forward (LAD→Hmax→fCover→LAI) per cluster + All, BUT each panel
# shows per-PLOT mean ΔTmax (one point per plot, instead of one point per day).
#
# Companion to fig_forward_HOBO_v10_fixed_per_cluster_plus_all.png (the daily
# version). 5 rows (P1..P4 + All) × 5 cols (forward steps).
#
# Row 1..4 : n = nb plots in cluster (14, 16, 8, 15)
# Row 5    : n = 53 (all plots)
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

# Cluster assignment per HOBO plot
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

# Pre-compute per-plot mean for each of the 5 coalitions of the fixed path
FIXED_BITS   <- c("0000", "0001", "0101", "0111", "1111")
FIXED_LABELS <- c("Baseline","+LAD","+LAD+H_max","+LAD+H_max+fCover","REF (+LAI)")

per_plot_mean <- function(bit) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) return(NULL)
  m <- merge(s, hobo_obs, by = c("id_plot", "date"))
  m <- merge(m, df_hobo_clu, by = "id_plot")
  m[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
        Delta_obs = mean(Delta_obs, na.rm = TRUE)),
    by = .(id_plot, Cluster)]
}

cli_h1("Per-plot mean for 5 fixed-path coalitions")
PP <- list()
for (b in FIXED_BITS) PP[[b]] <- per_plot_mean(b)

clusters <- c("P1", "P2", "P3", "P4", "All")

fit_subset <- function(d, cluster) {
  if (cluster == "All") dc <- d else dc <- d[Cluster == cluster]
  ok <- !is.na(dc$Delta_sim) & !is.na(dc$Delta_obs)
  x <- dc$Delta_obs[ok]; y <- dc$Delta_sim[ok]
  list(r = if (length(x) < 3) NA else suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       n = length(x))
}

cli_h1("Build long data + fits")
long <- list(); ann <- list()
for (cl in clusters) {
  for (k in seq_along(FIXED_BITS)) {
    bit <- FIXED_BITS[k]
    d <- PP[[bit]]
    if (is.null(d)) next
    panel_d <- if (cl == "All") copy(d) else d[Cluster == cl]
    panel_d[, Cluster_row := cl]; panel_d[, step := k]
    panel_d[, label := FIXED_LABELS[k]]
    long[[length(long) + 1L]] <- panel_d
    fit <- fit_subset(d, cl)
    # Harmonised single annotation block (one entry per stat, NA-safe)
    ann[[length(ann) + 1L]] <- data.table(
      Cluster_row = cl, step = k, label = FIXED_LABELS[k],
      r_val    = fit$r,
      rmse_val = fit$RMSE,
      n_val    = fit$n)
  }
}
long <- rbindlist(long, fill = TRUE)
long[, Cluster_row := factor(Cluster_row, levels = clusters)]
ann <- rbindlist(ann)
ann[, Cluster_row := factor(Cluster_row, levels = clusters)]
long$label <- factor(as.character(long$label), levels = FIXED_LABELS)
ann$label <- factor(as.character(ann$label), levels = FIXED_LABELS)

# Build single-line annotation strings ; skip r if NA
ann[, lab_r    := ifelse(is.na(r_val), NA_character_,
                          sprintf("italic(r)==%+.2f", r_val))]
ann[, lab_rmse := sprintf("RMSE==%.2f * ' °C'", rmse_val)]
ann[, lab_n    := sprintf("n==%d", n_val)]

xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

# Harmonised annotation font size / colour / face — single style for all stats
STAT_SIZE   <- 6.0      # all three stats same size
STAT_FACE   <- "bold"
STAT_COLOUR <- "grey15"

p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim, colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(alpha = 0.85, size = 3.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22, aes(group = 1)) +
  # r — only rendered if not NA
  geom_text(data = ann[!is.na(lab_r)], aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.10, vjust = 1.4,
            size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
            inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.10, vjust = 2.8,
            size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
            inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_n),
            parse = TRUE, hjust = -0.10, vjust = 4.2,
            size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
            inherit.aes = FALSE) +
  facet_grid(Cluster_row ~ label) +
  coord_fixed(xlim = xy_lim, ylim = xy_lim) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A"),
                       na.translate = FALSE) +
  scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18),
                      na.translate = FALSE) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 22),
        axis.title  = element_text(size = 22, face = "bold"),
        axis.text   = element_text(size = 16),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.png"),
        p, width = 24, height = 26, dpi = 600, bg = "white")
cli_alert_success("Saved fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.png (24x26 in @ 600 dpi)")

# Vector backup for max-quality slide insertion
ggsave(file.path(OUT, "fig_forward_HOBO_v10_fixed_per_cluster_plus_all_AVG.pdf"),
        p, width = 24, height = 26, device = cairo_pdf)
cli_alert_success("Saved PDF vector version")

# Stats CSV
cli_h1("Stats per (cluster, step)")
stats <- ann[, .(Cluster_row, step, label, lab_r, lab_rmse, lab_n)]
print(stats)
fwrite(stats, file.path(OUT, "tab_v10_HOBO_fixed_per_cluster_plus_all_AVG.csv"))
