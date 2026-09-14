# ==============================================================================
# Figure : true GREEDY forward feature selection on HOBO sims (V9).
# Two panel rows (or two figures) — one per criterion :
#   - Greedy(Tmax)  : Baseline → +LAI → +fCover → +Hmax → REF
#   - Greedy(slope) : Baseline → +fCover → +LAI → +Hmax → REF
#
# Each panel shows simulated metric vs observed (HOBO) metric across the
# 53 sensors, with r / RMSE / MAE annotated.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))

OUT <- here::here("outputs/figs_MEB2026_final")

cli_h1("Greedy forward selection — HOBO validation (V9)")

df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)
hobo_root   <- here::here("out_files/musica_hobo_v9")

# Helper : per-plot slope or Tmax mean for a given coalition
collect_step <- function(bit, metric = c("slope", "Tmax")) {
  metric <- match.arg(metric)
  subdir <- file.path(hobo_root, bit)
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

# HOBO obs : slope hourly + mean Delta_Tmax
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
hobo_dT_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                              df_macro, CFG$ids_to_remove))[
  , .(dT_obs = mean(Delta_obs, na.rm = TRUE)), by = id_plot]
hobo_obs <- merge(hobo_slope_obs, hobo_dT_obs, by = "id_plot")
cli_alert("HOBO obs rows : {nrow(hobo_obs)}")

# Path definitions
paths <- list(
  Greedy_Tmax  = list(
    name   = "Greedy (Tmax-based)",
    bits   = c("0000","1000","1010","1110","1111"),
    labels = c("Baseline","+LAI","+LAI+fCover","+LAI+fCover+H_max","REF"),
    metric = "Tmax",
    obs_col = "dT_obs",
    sim_lab = "Simulated mean ΔTmax (°C)",
    obs_lab = "Field mean ΔTmax (°C)"
  ),
  Greedy_slope = list(
    name   = "Greedy (slope-based)",
    bits   = c("0000","0010","1010","1110","1111"),
    labels = c("Baseline","+fCover","+fCover+LAI","+fCover+LAI+H_max","REF"),
    metric = "slope",
    obs_col = "slope_obs",
    sim_lab = "Simulated slope",
    obs_lab = "Field slope"
  )
)

make_path_fig <- function(p_cfg, png_name) {
  long <- list(); fits <- list()
  for (k in seq_along(p_cfg$bits)) {
    bit <- p_cfg$bits[k]
    d <- collect_step(bit, p_cfg$metric)
    if (nrow(d) == 0) {
      cli_alert_warning("bit {bit} : no NCs found in {file.path(hobo_root, bit)}")
      next
    }
    d <- merge(d, hobo_obs[, .(id_plot, obs_metric = get(p_cfg$obs_col))],
                by = "id_plot")
    d[, step := k]; d[, label := p_cfg$labels[k]]; d[, bit := bit]
    long[[k]] <- d
    ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
    x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
    fits[[k]] <- data.table(step = k, label = p_cfg$labels[k],
                             r = suppressWarnings(cor(x, y)),
                             RMSE = sqrt(mean((y - x)^2)),
                             MAE  = mean(abs(y - x)),
                             bias = mean(y - x),
                             n = length(x))
  }
  long <- rbindlist(long, fill = TRUE)
  fits <- rbindlist(fits)
  long$label <- factor(as.character(long$label), levels = p_cfg$labels)
  cli_alert("Per-step fits ({p_cfg$name}) :"); print(fits)
  fwrite(fits, file.path(OUT, sprintf("tab_v9_HOBO_%s.csv",
                                       gsub(" ", "_", tolower(p_cfg$name)))))
  xy_lim <- range(c(long$sim_metric, long$obs_metric), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.06
  xy_lim <- xy_lim + c(-xy_pad, xy_pad)
  ann <- fits[, .(label,
                    lab_r    = sprintf("italic(r)==%+.2f", r),
                    lab_rmse = sprintf("RMSE==%.3f", RMSE),
                    lab_mae  = sprintf("MAE==%.3f", MAE))]
  ann$label <- factor(ann$label, levels = p_cfg$labels)

  p <- ggplot(long, aes(x = obs_metric, y = sim_metric)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = if (p_cfg$metric == "slope") 1 else 0,
                linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = if (p_cfg$metric == "slope") 1 else 0,
                linetype = "dotted", colour = "grey80") +
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
    labs(title = p_cfg$name, x = p_cfg$obs_lab, y = p_cfg$sim_lab) +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 12),
           plot.title = element_text(face = "bold", size = 14),
           panel.grid.minor = element_blank())
  ggsave(file.path(OUT, png_name), p, width = 18, height = 5, dpi = 300)
  cli_alert_success("Saved {png_name}")
  invisible(list(long = long, fits = fits))
}

cli_h2("Greedy (Tmax)")
res_t <- make_path_fig(paths$Greedy_Tmax,
                        "fig_forward_selection_HOBO_v9_greedy_Tmax.png")

cli_h2("Greedy (slope)")
res_s <- make_path_fig(paths$Greedy_slope,
                        "fig_forward_selection_HOBO_v9_greedy_slope.png")

saveRDS(list(Tmax = res_t, slope = res_s),
         here::here("outputs/v9_greedy_forward_HOBO.rds"))
cli_h1("Greedy HOBO figures done.")
