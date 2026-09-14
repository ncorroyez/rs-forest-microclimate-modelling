# ==============================================================================
# Forward greedy PER CLUSTER — each row = one cluster (P1..P4), each col = one
# forward step. Within each cluster, the greedy order is recomputed from that
# cluster's daily HOBO data only.
#
# Aggregation : daily × cluster (mean across plots within cluster per day).
# Output : fig_forward_HOBO_v10_greedy_per_cluster.png (4 × 5 grid)
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

# Cluster assignment
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

# Per-cluster daily aggregation
agg_cluster_daily <- function(sim_long) {
  m <- merge(sim_long, hobo_obs, by = c("id_plot", "date"))
  m <- merge(m, df_hobo_clu, by = "id_plot")
  m[, .(Delta_sim = mean(Delta_sim, na.rm = TRUE),
        Delta_obs = mean(Delta_obs, na.rm = TRUE)),
    by = .(Cluster, date)]
}

# Pre-compute daily-cluster series for all 16 coalitions
cli_h1("Pre-computing 16 coalitions at daily × cluster")
COAL_DT <- list()
for (bit in unname(.ARCH_COAL_MAP)) {
  s <- collect_step_daily(bit)
  if (nrow(s) == 0) { COAL_DT[[bit]] <- NULL; next }
  COAL_DT[[bit]] <- agg_cluster_daily(s)
}
cli_alert("Coalitions pre-computed : {length(COAL_DT)}")

# Fit per (bit, cluster)
fit_bc <- function(bit, cluster) {
  d <- COAL_DT[[bit]]
  if (is.null(d)) return(list(r = NA, RMSE = NA, MAE = NA, bias = NA, n = 0))
  dc <- d[Cluster == cluster]
  ok <- !is.na(dc$Delta_sim) & !is.na(dc$Delta_obs)
  x <- dc$Delta_obs[ok]; y <- dc$Delta_sim[ok]
  if (length(x) < 5) return(list(r = NA, RMSE = NA, MAE = NA, bias = NA, n = length(x)))
  list(r = suppressWarnings(cor(x, y)),
       RMSE = sqrt(mean((y - x)^2)),
       MAE = mean(abs(y - x)),
       bias = mean(y - x), n = length(x))
}

# Greedy path per cluster
VARS <- c("LAI", "Hmax", "fCover", "LAD")
bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")

greedy_for_cluster <- function(cluster) {
  active <- character(0); path <- list()
  b0 <- bit_from_set(active); f0 <- fit_bc(b0, cluster)
  path[[1]] <- data.table(step = 1, Cluster = cluster, added = "Baseline",
                          active_set = "", bit = b0,
                          r = f0$r, RMSE = f0$RMSE, MAE = f0$MAE,
                          bias = f0$bias, n = f0$n)
  for (k in 2:5) {
    remaining <- setdiff(VARS, active)
    cand_bits <- sapply(remaining, function(v) bit_from_set(c(active, v)))
    cand_r <- sapply(cand_bits, function(b) fit_bc(b, cluster)$r)
    # If all NA fall back to any
    if (all(is.na(cand_r))) break
    best_i <- which.max(cand_r)
    best_v <- remaining[best_i]
    active <- c(active, best_v); bit <- cand_bits[best_i]
    fit <- fit_bc(bit, cluster)
    path[[k]] <- data.table(step = k, Cluster = cluster, added = best_v,
                            active_set = paste(active, collapse = "+"),
                            bit = bit, r = fit$r, RMSE = fit$RMSE,
                            MAE = fit$MAE, bias = fit$bias, n = fit$n)
    if (length(active) == length(VARS)) break
  }
  rbindlist(path)
}

clusters <- c("P1", "P2", "P3", "P4")
greedy_all <- rbindlist(lapply(clusters, greedy_for_cluster))
cli_alert("Greedy paths per cluster :")
print(greedy_all[, .(Cluster, step, added, active_set, r = round(r, 3),
                     RMSE = round(RMSE, 3))])
fwrite(greedy_all, file.path(OUT, "tab_v10_HOBO_greedy_per_cluster.csv"))

# Build long data for plot : need to recall the daily×cluster points for each
# (cluster, bit) in the greedy path
long <- list()
for (i in seq_len(nrow(greedy_all))) {
  r <- greedy_all[i]
  d <- COAL_DT[[r$bit]]
  if (is.null(d)) next
  dc <- d[Cluster == r$Cluster][, step := r$step][, label := r$added]
  long[[i]] <- dc
}
long <- rbindlist(long, fill = TRUE)
long[, Cluster := factor(Cluster, levels = clusters)]

# Build step-specific labels per cluster
greedy_all[, panel_label := {
  if (step == 1) "Baseline"
  else if (active_set == paste(VARS, collapse = "+") |
           length(strsplit(active_set, "\\+")[[1]]) == 4) {
    sprintf("REF (+%s)", added)
  } else paste("+", gsub("Hmax", "H_max",
                          paste(strsplit(active_set, "\\+")[[1]], collapse = " + ")))
}, by = .(Cluster, step)]

long <- merge(long,
              greedy_all[, .(Cluster, step, panel_label)],
              by = c("Cluster", "step"))

# Annotation per (Cluster, step)
ann <- greedy_all[, .(Cluster, step, panel_label,
                       lab_r    = sprintf("italic(r)==%+.2f", r),
                       lab_rmse = sprintf("RMSE==%.2f * ' °C'", RMSE))]

# Plot 4 × 5 grid
xy_lim <- range(c(long$Delta_sim, long$Delta_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

# Make a panel label combining step + variable for each (cluster, step)
long[, panel := sprintf("Step %d : %s", step, panel_label)]
ann[, panel := sprintf("Step %d : %s", step, panel_label)]

p <- ggplot(long, aes(x = Delta_obs, y = Delta_sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(aes(colour = Cluster), alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey25", linewidth = 0.7) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 3.5,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.6, size = 3.2,
            inherit.aes = FALSE, colour = "grey15") +
  facet_grid(Cluster ~ step, scales = "fixed",
              labeller = labeller(step = function(x) {
                sapply(x, function(s) {
                  sub <- ann[step == as.integer(s)]
                  if (nrow(sub) == 0) return(paste("Step", s))
                  # Different label per cluster — use the cluster's own
                  paste("Step", s)
                })
              })) +
  coord_fixed(xlim = xy_lim, ylim = xy_lim) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A")) +
  labs(title = "V10 — Greedy forward selection PER CLUSTER (daily × cluster, V10 Tmax)",
       subtitle = "Each row = cluster, columns = forward step. Greedy order recomputed within each cluster.",
       x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "none",
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_forward_HOBO_v10_greedy_per_cluster.png"),
        p, width = 16, height = 13, dpi = 300)
cli_alert_success("Saved fig_forward_HOBO_v10_greedy_per_cluster.png")

# Also a compact summary table : greedy order per cluster
cli_h1("Greedy order per cluster (summary)")
summary_table <- greedy_all[, .(order = paste(added[-1], collapse = " → "),
                                  r_step2 = r[2], r_REF = r[5],
                                  RMSE_REF = RMSE[5]),
                              by = Cluster]
print(summary_table)
fwrite(summary_table, file.path(OUT, "tab_v10_HOBO_greedy_order_per_cluster.csv"))
