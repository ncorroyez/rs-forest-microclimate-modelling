# ==============================================================================
# Daily × cluster aggregation — both HOBO greedy forward and the by-cluster
# validation scatter.
#
# Aggregation rule :
#   - Underlying data : daily ΔTmax per (plot, day)
#   - For each cluster × day, average across plots in that cluster
#   - 4 clusters × ~120 days = ~480 points per panel
#
# Outputs :
#   1. fig_forward_HOBO_v10_true_greedy_Tmax_daily_cluster.png
#   2. fig_HOBO_validation_v10_daily_cluster_panels.png  (4 panneaux P1..P4)
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

# ---- Cluster assignment per HOBO plot (nearest forest cell) ----------------
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
cli_alert("HOBO clusters : "); print(table(df_hobo_clu$Cluster))

# ---- HOBO obs daily (per plot × date) ---------------------------------------
hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                           df_macro, CFG$ids_to_remove))[
  , .(id_plot, date, Delta_obs)]

# ---- Sim daily for one coalition (V10 fallback V9) -------------------------
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
    d <- as.data.table(res)[, .(id_plot = id, date, Delta_sim = Delta_Tmax,
                                 Tmax_macro)]
    rows[[id]] <- d
  }
  rbindlist(rows)
}

# ---- Aggregate per (cluster, date) -----------------------------------------
agg_cluster_daily <- function(sim_long) {
  m <- merge(sim_long, hobo_obs, by = c("id_plot", "date"))
  m <- merge(m, df_hobo_clu, by = "id_plot")
  m[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
        Delta_obs = mean(Delta_obs, na.rm = TRUE),
        Tmax_macro = mean(Tmax_macro, na.rm = TRUE),
        n_plots = .N),
    by = .(Cluster, date)]
}

# ============================================================================
# 1. Greedy forward HOBO Tmax (daily × cluster)
# ============================================================================
cli_h1("Greedy forward HOBO Tmax (daily × cluster)")

# Fit all 16 coalitions at the daily×cluster level
fit_bit_cluster <- function(bit) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) return(list(r = NA, RMSE = NA, MAE = NA, bias = NA, n = 0))
  m <- agg_cluster_daily(s)
  ok <- !is.na(m$Delta_sim) & !is.na(m$Delta_obs)
  x <- m$Delta_obs[ok]; y <- m$Delta_sim[ok]
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       bias = mean(y - x),
       n = length(x))
}

bit_table <- data.table(bit = unname(.ARCH_COAL_MAP))
fits <- list()
for (b in bit_table$bit) fits[[b]] <- fit_bit_cluster(b)
all_fits <- data.table(bit = names(fits),
                       r    = sapply(fits, `[[`, "r"),
                       RMSE = sapply(fits, `[[`, "RMSE"),
                       MAE  = sapply(fits, `[[`, "MAE"),
                       bias = sapply(fits, `[[`, "bias"),
                       n    = sapply(fits, `[[`, "n"))
print(all_fits)
fwrite(all_fits, file.path(OUT, "tab_v10_HOBO_all16_fits_Tmax_daily_cluster.csv"))

# Greedy path
VARS <- c("LAI", "Hmax", "fCover", "LAD")
bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")
active <- character(0)
path <- list()
b0 <- bit_from_set(active); f0 <- fits[[b0]]
path[[1]] <- data.table(step = 1, added = "Baseline", active_set = "",
                        bit = b0, r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE,
                        bias = f0$bias, n = f0$n)
for (k in 2:5) {
  remaining <- setdiff(VARS, active)
  cand_bits <- sapply(remaining, function(v) bit_from_set(c(active, v)))
  cand_r <- sapply(cand_bits, function(b) fits[[b]]$r)
  best_i <- which.max(cand_r)
  best_v <- remaining[best_i]
  active <- c(active, best_v); bit <- cand_bits[best_i]
  path[[k]] <- data.table(step = k, added = best_v,
                          active_set = paste(active, collapse = "+"),
                          bit = bit, r = cand_r[best_i],
                          RMSE = fits[[bit]]$RMSE, MAE = fits[[bit]]$MAE,
                          bias = fits[[bit]]$bias, n = fits[[bit]]$n)
  if (length(active) == length(VARS)) break
}
greedy <- rbindlist(path)
print(greedy)
fwrite(greedy, file.path(OUT, "tab_v10_HOBO_true_greedy_Tmax_daily_cluster.csv"))

labels_path <- sapply(seq_len(nrow(greedy)), function(k) {
  if (k == 1) "Baseline"
  else if (k == nrow(greedy)) sprintf("REF (+%s)", greedy$added[k])
  else paste("+", paste(strsplit(greedy$active_set[k], "\\+")[[1]], collapse = " + "))
})
labels_path <- gsub("Hmax", "H_max", labels_path)

long <- list()
for (k in seq_len(nrow(greedy))) {
  bit <- greedy$bit[k]
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) next
  m <- agg_cluster_daily(s)
  m[, step := k]; m[, label := labels_path[k]]
  long[[k]] <- m
}
long <- rbindlist(long, fill = TRUE)
long$label <- factor(as.character(long$label), levels = labels_path)

ann <- greedy[, .(label = labels_path,
                   lab_r    = sprintf("italic(r)==%+.2f", r),
                   lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                   lab_mae  = sprintf("MAE==%.2f * ' °C'", MAE))]
ann$label <- factor(ann$label, levels = labels_path)
xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

p_g <- ggplot(long, aes(x = Delta_obs, y = Delta_sim, colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(alpha = 0.55, size = 1.8) +
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
  labs(title = "V10 — TRUE greedy forward HOBO ΔTmax (daily × cluster aggregation)",
       x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave(file.path(OUT, "fig_forward_HOBO_v10_true_greedy_Tmax_daily_cluster.png"),
        p_g, width = 18, height = 6, dpi = 300)
cli_alert_success("Saved fig_forward_HOBO_v10_true_greedy_Tmax_daily_cluster.png")

# ============================================================================
# 2. Validation by-cluster (REF only, 4-panel, daily × cluster)
# ============================================================================
cli_h1("Validation by-cluster — daily × cluster (REF=1111)")

s_ref <- collect_step_daily("1111")
m_ref <- agg_cluster_daily(s_ref)

stats_C <- m_ref[, {
  ok <- !is.na(Delta_obs) & !is.na(Delta_sim)
  x <- Delta_obs[ok]; y <- Delta_sim[ok]
  .(r = suppressWarnings(cor(x, y)),
    RMSE = sqrt(mean((y - x)^2)),
    MAE = mean(abs(y - x)),
    bias = mean(y - x), n = length(x))
}, by = Cluster]
print(stats_C)
fwrite(stats_C, file.path(OUT, "tab_v10_HOBO_validation_REF_daily_cluster_stats.csv"))

ann_C <- stats_C[, .(Cluster,
                      lab_r    = sprintf("italic(r)==%+.2f", r),
                      lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE),
                      lab_n    = sprintf("n==%d", n))]
xy_lim2 <- range(c(m_ref$Delta_sim, m_ref$Delta_obs), na.rm = TRUE)
xy_pad2 <- diff(xy_lim2) * 0.05; xy_lim2 <- xy_lim2 + c(-xy_pad2, xy_pad2)

p_cl <- ggplot(m_ref, aes(x = Delta_obs, y = Delta_sim, colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey25", linewidth = 0.9) +
  geom_text(data = ann_C, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann_C, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8, size = 5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann_C, aes(x = -Inf, y = Inf, label = lab_n),
            parse = TRUE, hjust = -0.08, vjust = 4.2, size = 4.2,
            inherit.aes = FALSE, colour = "grey35") +
  facet_wrap(~ Cluster, nrow = 1) +
  coord_fixed(xlim = xy_lim2, ylim = xy_lim2) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A")) +
  scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)) +
  labs(title = "V10 REF — HOBO validation by cluster (daily × cluster aggregation)",
       x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  guides(colour = "none", shape = "none") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_HOBO_validation_v10_daily_cluster_panels.png"),
        p_cl, width = 14, height = 4.5, dpi = 300)
cli_alert_success("Saved fig_HOBO_validation_v10_daily_cluster_panels.png")
