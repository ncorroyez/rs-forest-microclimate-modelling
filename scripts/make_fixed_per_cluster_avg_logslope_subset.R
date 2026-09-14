# ==============================================================================
# Per-plot AVG fig — log(slope), SUBSET version (lightweight for slides).
# Same structure as make_fixed_per_cluster_avg_logslope.R but only displays
# selected rows (default : P1, P4, All) instead of P1..P4 + All.
#
# Configure CLUSTERS_SUBSET at the top to pick the rows you want.
#
# Output : fig_forward_HOBO_v10_fixed_subset_AVG_logslope.{png,pdf}
# ==============================================================================

# ---- USER CONFIG : which rows to display (subset of P1..P4, All) -----------
CLUSTERS_SUBSET <- c("P1", "P4", "All")

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
df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

stopifnot(all(CLUSTERS_SUBSET %in% c("P1", "P2", "P3", "P4", "All")))

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

# ---- HOBO observed slope per plot ------------------------------------------
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
hobo_obs <- hobo_with_era5[, {
  fit <- lm(t_hobo ~ Tair_era5)
  sl  <- as.numeric(coef(fit)[2])
  .(slope_obs = sl, log_slope_obs = log(sl))
}, by = id_plot]

# ---- Simulated slope per plot per coalition --------------------------------
collect_step_slope <- function(bit) {
  v10 <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                              z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- data.table(id_plot = id,
                              slope_sim     = res$slope,
                              log_slope_sim = log(res$slope))
  }
  rbindlist(rows)
}

cli_h1("Pre-computing all 16 HOBO coalitions (per-plot log slope)")
PP <- list()
for (bit in unname(.ARCH_COAL_MAP)) {
  s <- collect_step_slope(bit)
  if (nrow(s) == 0) next
  m <- merge(s, hobo_obs[, .(id_plot, log_slope_obs)], by = "id_plot")
  m <- merge(m, df_hobo_clu, by = "id_plot")
  PP[[bit]] <- m
}

# ---- Δv per cluster (from heatmap Tmax) ------------------------------------
sc_name_of <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
arch_nc_path <- function(cl, bit) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  v10 <- file.path("out_files/H1_archetypes_v10_fcovmean", arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_name_of[bit]))
  v9  <- file.path("out_files/H1_archetypes_v9", arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_name_of[bit]))
  if (file.exists(v10)) v10 else v9
}
arch_tmax <- function(cl, bit) {
  res <- tryCatch(extract_deltatmax_one(arch_nc_path(cl, bit), df_macro,
                                          CFG$date_seq,
                                          z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}
VARS     <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
arch_clusters <- c("3", "2", "4", "1")
new_labels    <- c("P1", "P2", "P3", "P4")
dv_per_cluster <- list()
for (i in seq_along(arch_clusters)) {
  cl_old <- arch_clusters[i]; profile <- new_labels[i]
  tmax_ref <- arch_tmax(cl_old, "1111")
  dv <- sapply(VARS, function(v) tmax_ref - arch_tmax(cl_old, loo_bits[v]))
  dv_per_cluster[[profile]] <- dv
}
dv_avg <- colMeans(do.call(rbind, dv_per_cluster), na.rm = TRUE)
dv_per_cluster[["All"]] <- dv_avg
order_per_cluster <- lapply(dv_per_cluster, function(dv) names(sort(dv, decreasing = TRUE)))

# ---- Subset of rows ---------------------------------------------------------
clusters <- CLUSTERS_SUBSET
cli_alert("Displaying subset : {paste(clusters, collapse=' + ')}")

bit_from_set <- function(set) paste(ifelse(VARS %in% set, "1", "0"), collapse = "")
to_pm <- function(v) ifelse(v == "Hmax", "italic(H)[max]", v)
paths <- list()
for (cl in clusters) {
  ord <- order_per_cluster[[cl]]
  bits <- c(bit_from_set(character(0)))
  active <- character(0)
  for (v in ord) { active <- c(active, v); bits <- c(bits, bit_from_set(active)) }
  labels_step <- c("'Baseline'",
                    sapply(seq_along(ord), function(k) {
                      cum_vars <- sapply(ord[1:k], to_pm)
                      joined <- paste(cum_vars, collapse = "*' + '*")
                      if (k == length(ord)) paste0("'REF: '*", joined)
                      else paste0("'+ '*", joined)
                    }))
  paths[[cl]] <- data.table(step = seq_along(bits), Cluster = cl,
                             bit = bits, label = labels_step)
}

# ---- Build long data + per-(cluster, step) fits ----------------------------
cli_h1("Build long data + fits")
long <- list(); ann <- list()
for (cl in clusters) {
  p <- paths[[cl]]
  for (k in seq_len(nrow(p))) {
    bit <- p$bit[k]
    d <- PP[[bit]]
    if (is.null(d)) next
    panel_d <- if (cl == "All") copy(d) else d[Cluster == cl]
    panel_d[, Cluster_row := cl]; panel_d[, step := k]
    panel_d[, label := p$label[k]]
    long[[length(long) + 1L]] <- panel_d
    ok <- !is.na(panel_d$log_slope_sim) & !is.na(panel_d$log_slope_obs)
    x <- panel_d$log_slope_obs[ok]; y <- panel_d$log_slope_sim[ok]
    n <- length(x)
    r_val <- if (n < 3) NA_real_ else suppressWarnings(cor(x, y))
    rmse_val <- if (n == 0) NA_real_ else sqrt(mean((y - x)^2))
    mae_val  <- if (n == 0) NA_real_ else mean(abs(y - x))
    ann[[length(ann) + 1L]] <- data.table(
      Cluster_row = cl, step = k, label = p$label[k],
      r_val = r_val, rmse_val = rmse_val, mae_val = mae_val, n_val = n)
  }
}
long <- rbindlist(long, fill = TRUE)
ann  <- rbindlist(ann)
long[, Cluster_row := factor(Cluster_row, levels = clusters)]
ann[ , Cluster_row := factor(Cluster_row, levels = clusters)]

ann[, lab_r    := ifelse(is.na(r_val), NA_character_,
                          sprintf("italic(r)==%+.2f", r_val))]
ann[, lab_rmse := sprintf("RMSE==%.3f", rmse_val)]
ann[, lab_mae  := sprintf("MAE==%.3f", mae_val)]
ann[, lab_n    := sprintf("n==%d", n_val)]
ann[, panel_lab := label]

xy_lim <- range(c(long$log_slope_sim, long$log_slope_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.05; xy_lim <- xy_lim + c(-xy_pad, xy_pad)

STAT_SIZE   <- 5.5
STAT_FACE   <- "bold"
STAT_COLOUR <- "grey15"

# Layout : 16:9 slide ; height scales with #rows so panels stay square-ish.
n_rows <- length(clusters)
fig_h  <- max(6.5, 4.5 * n_rows)        # 1 row ≈ 4.5 in tall, min 6.5 in
fig_w  <- 24

p <- ggplot(long, aes(x = log_slope_obs, y = log_slope_sim,
                       colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(alpha = 0.85, size = 3.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22, aes(group = 1)) +
  geom_text(data = ann[!is.na(lab_r)],
             aes(x = -Inf, y = Inf, label = lab_r),
             parse = TRUE, hjust = -0.10, vjust = 1.4,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
             parse = TRUE, hjust = -0.10, vjust = 2.8,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
             parse = TRUE, hjust = -0.10, vjust = 4.2,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_n),
             parse = TRUE, hjust = -0.10, vjust = 5.6,
             size = STAT_SIZE, fontface = STAT_FACE, colour = STAT_COLOUR,
             inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = Inf, y = -Inf, label = panel_lab),
             hjust = 1.05, vjust = -0.5,
             size = 4.5, fontface = "italic", parse = TRUE,
             inherit.aes = FALSE, colour = "grey25") +
  facet_grid(Cluster_row ~ step,
              labeller = labeller(step = c("1" = "Baseline",
                                              "2" = "Step 1",
                                              "3" = "Step 2",
                                              "4" = "Step 3",
                                              "5" = "Step 4"))) +
  coord_fixed(xlim = xy_lim, ylim = xy_lim) +
  scale_colour_manual(values = c(PAL_CLUSTER), na.translate = FALSE) +
  scale_shape_manual(values = c(PAL_SHAPES), na.translate = FALSE) +
  labs(x = expression("Field " ~ log(slope[micro/macro])),
       y = expression("Simulated " ~ log(slope[micro/macro]))) +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text  = element_text(face = "bold", size = 22),
        axis.title  = element_text(face = "bold", size = 24),
        axis.text   = element_text(size = 16),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 24, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        panel.grid = element_blank())

suffix <- paste(clusters, collapse = "_")
ggsave(file.path(OUT, sprintf("fig_forward_HOBO_v10_subset_%s_AVG_logslope.png", suffix)),
        p, width = fig_w, height = fig_h, dpi = 600, bg = "white", limitsize = FALSE)
ggsave(file.path(OUT, sprintf("fig_forward_HOBO_v10_subset_%s_AVG_logslope.pdf", suffix)),
        p, width = fig_w, height = fig_h, device = cairo_pdf, limitsize = FALSE)
cli_alert_success("Saved fig_forward_HOBO_v10_subset_{suffix}_AVG_logslope (png + pdf)")

fwrite(ann[, .(Cluster_row, step, label,
                r = round(r_val, 3), RMSE = round(rmse_val, 4),
                MAE = round(mae_val, 4), n = n_val)],
        file.path(OUT, sprintf("tab_v10_HOBO_subset_%s_AVG_logslope.csv", suffix)))
