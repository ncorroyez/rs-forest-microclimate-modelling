# ==============================================================================
# Hourly MuSICA (y) vs HOBO (x), REF coalition (1111) — V11 (model-3.2.3).
# Same convention as V9 fig : nair == 1, no time shift.
# Two figs : all 53 plots + buffering subset.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))

OUT <- here::here("outputs/figs_MEB2026_final")

cli_h1("Hourly MuSICA vs HOBO REF — V11 (model-3.2.3)")

hobo_obs <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a",
         as.Date(datetime) %in% CFG$date_seq,
         !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(id_plot, time) %>%
  summarise(Tair_obs = mean(t_hobo, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()
cli_alert("HOBO hourly rows : {nrow(hobo_obs)} | {length(unique(hobo_obs$id_plot))} plots")

# MuSICA V11 outputs (model-3.2.3)
sim_dir <- here::here("out_files/musica_hobo_v11_model323/1111")
nc_files <- list.files(sim_dir, pattern = "musica_out_HOBO_.+[.]nc$", full.names = TRUE)
cli_alert("REF NCs : {length(nc_files)}")

sim_rows <- list()
for (f in nc_files) {
  id <- sub("musica_out_HOBO_(.+)[.]nc$", "\\1", basename(f))
  nc <- tryCatch(nc_open(f), error = function(e) NULL); if (is.null(nc)) next
  raw <- tryCatch(musica.tools::get_variable(nc, "Tair_z"), error = function(e) NULL)
  nc_close(nc)
  if (is.null(raw) || nrow(raw) == 0) next
  d <- as.data.table(raw)[nair == 1,
                            .(time = lubridate::floor_date(time, "hour"),
                              Tair_sim = Tair_z - 273.15)]
  d[, id_plot := id]
  d <- d[as.Date(time) %in% CFG$date_seq]
  sim_rows[[id]] <- d
}
DT_sim <- rbindlist(sim_rows)

DT <- merge(DT_sim, hobo_obs, by = c("id_plot", "time"))
DT <- DT[!is.na(Tair_sim) & !is.na(Tair_obs)]
cli_alert("Paired rows : {nrow(DT)} | {length(unique(DT$id_plot))} plots")

era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)
DT_era <- merge(DT, as.data.table(era5_hourly), by = "time")
slope_per_plot <- DT_era[, .(slope = as.numeric(coef(lm(Tair_sim ~ Tair_era5))[2])),
                          by = id_plot]
buf_ids <- slope_per_plot[slope < 1, id_plot]
amp_ids <- slope_per_plot[slope > 1, id_plot]
cli_alert("V11 buffering : {length(buf_ids)} | amplifying : {length(amp_ids)}")

stats_block <- function(dt) {
  x <- dt$Tair_obs; y <- dt$Tair_sim
  list(r = cor(x, y), r2 = cor(x, y)^2,
       rmse = sqrt(mean((y - x)^2)),
       mae = mean(abs(y - x)),
       bias = mean(y - x), n = nrow(dt))
}

make_fig <- function(dt, stats, png_path, subtitle) {
  lim <- range(c(dt$Tair_obs, dt$Tair_sim), na.rm = TRUE)
  pad <- diff(lim) * 0.04; lim <- lim + c(-pad, pad)
  p <- ggplot(dt, aes(x = Tair_obs, y = Tair_sim)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(option = "viridis", trans = "log10",
                          name = "Count", labels = scales::label_log()) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "white", linewidth = 1.0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#FFB400", linewidth = 1.0) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 1.5,
              label = sprintf("italic(r)==%+.3f", stats$r),
              parse = TRUE, size = 8, fontface = "bold", colour = "grey15") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 3.0,
              label = sprintf("italic(R)^2==%.3f", stats$r2),
              parse = TRUE, size = 8, fontface = "bold", colour = "grey15") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 4.5,
              label = sprintf("RMSE==%.2f * ' °C'", stats$rmse),
              parse = TRUE, size = 8, fontface = "bold", colour = "grey15") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 6.0,
              label = sprintf("MAE==%.2f * ' °C'", stats$mae),
              parse = TRUE, size = 8, fontface = "bold", colour = "grey15") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 7.5,
              label = sprintf("bias==%+.2f * ' °C'", stats$bias),
              parse = TRUE, size = 7, fontface = "italic", colour = "grey25") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.10, vjust = 9.0,
              label = sprintf("n = %s (hourly)", format(stats$n, big.mark = " ")),
              size = 7, fontface = "italic", colour = "grey35") +
    coord_fixed(xlim = lim, ylim = lim) +
    labs(x = expression("Field hourly " * italic(T)[air] ~ "(HOBO, °C)"),
         y = expression("Simulated hourly " * italic(T)[air] ~ "(MuSICA v11, °C)"),
         subtitle = subtitle) +
    theme_bw(base_size = 24) +
    theme(axis.title  = element_text(face = "bold", size = 26),
          axis.text   = element_text(size = 20),
          plot.subtitle = element_text(size = 16, face = "italic", colour = "grey35"),
          legend.text  = element_text(size = 20),
          legend.title = element_text(size = 22, face = "bold"),
          legend.key.height = unit(2, "cm"),
          legend.position = "right",
          panel.grid = element_blank())
  ggsave(png_path, p, width = 12, height = 10, dpi = 600, bg = "white")
  ggsave(sub("[.]png$", ".pdf", png_path), p, width = 12, height = 10,
          device = cairo_pdf)
  cli_alert_success("Saved {basename(png_path)} (+ pdf)")
}

st_all <- stats_block(DT)
cli_alert("ALL : r={round(st_all$r,4)} R²={round(st_all$r2,4)} RMSE={round(st_all$rmse,3)}°C MAE={round(st_all$mae,3)}°C bias={round(st_all$bias,3)}°C")
make_fig(DT, st_all,
          file.path(OUT, "fig_hourly_sim_vs_obs_REF_v11.png"),
          sprintf("V11 (model-3.2.3) — all %d HOBO plots", length(unique(DT$id_plot))))

DT_buf <- DT[id_plot %in% buf_ids]
if (nrow(DT_buf) > 0) {
  st_buf <- stats_block(DT_buf)
  cli_alert("BUF : r={round(st_buf$r,4)} R²={round(st_buf$r2,4)} RMSE={round(st_buf$rmse,3)}°C bias={round(st_buf$bias,3)}°C")
  make_fig(DT_buf, st_buf,
            file.path(OUT, "fig_hourly_sim_vs_obs_REF_v11_buf.png"),
            sprintf("V11 — buffering plots only (slope < 1, n=%d)", length(buf_ids)))
}

if (length(amp_ids) > 0) {
  DT_amp <- DT[id_plot %in% amp_ids]
  st_amp <- stats_block(DT_amp)
  make_fig(DT_amp, st_amp,
            file.path(OUT, "fig_hourly_sim_vs_obs_REF_v11_amp.png"),
            sprintf("V11 — amplifying plots only (slope > 1, n=%d)", length(amp_ids)))
}

fwrite(rbind(
  data.table(version = "v11", subset = "all", n_plots = length(unique(DT$id_plot)),
              r = st_all$r, R2 = st_all$r2, RMSE = st_all$rmse,
              MAE = st_all$mae, bias = st_all$bias, n_hourly = st_all$n),
  if (length(buf_ids) > 0) {
    data.table(version = "v11", subset = "buffering", n_plots = length(buf_ids),
                r = st_buf$r, R2 = st_buf$r2, RMSE = st_buf$rmse,
                MAE = st_buf$mae, bias = st_buf$bias, n_hourly = st_buf$n)
  } else NULL
), file.path(OUT, "tab_hourly_sim_vs_obs_REF_v11_stats.csv"))
