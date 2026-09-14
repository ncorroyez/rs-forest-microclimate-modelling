# ==============================================================================
# Inverse forward HOBO — Tmax-based — with V10 baseline (fCover_b = 0.869)
#
# Path : Baseline → +LAD → +LAD+Hmax → +LAD+Hmax+fCover → REF (+LAI)
# Coalitions :
#   0000  fCover=0  → V10 NCs (baseline-dependent)
#   0001  fCover=0  → V10 NCs
#   0101  fCover=0  → V10 NCs
#   0111  fCover=1  → V9 NCs (uses real fCover, unaffected by baseline)
#   1111  fCover=1  → V9 NCs (REF)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# V10 path : use V10 NC if exists, else V9
collect_step_v10_tmax <- function(bit) {
  v10_dir <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  v9_dir  <- file.path("out_files/musica_hobo_v9", bit)
  subdir <- if (dir.exists(v10_dir) && length(list.files(v10_dir, "\\.nc$")) > 0) v10_dir else v9_dir
  cli_alert("bit {bit} -> {subdir}")
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

hobo_dT_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                              df_macro, CFG$ids_to_remove))[
  , .(obs_metric = mean(Delta_obs, na.rm = TRUE)), by = id_plot]

inv_tmax <- list(
  bits   = c("0000","0001","0101","0111","1111"),
  # Labels are plotmath expressions so H_max renders italic H + max subscript.
  labels = c("Baseline",
              "'+LAD'",
              "'+LAD +'~italic(H)[max]",
              "'+LAD +'~italic(H)[max]~'+ fCover'",
              "'REF: +LAD +'~italic(H)[max]~'+ fCover + LAI'"))

long <- list(); fits <- list()
for (k in seq_along(inv_tmax$bits)) {
  bit <- inv_tmax$bits[k]
  d <- collect_step_v10_tmax(bit)
  if (nrow(d) == 0) next
  d <- merge(d, hobo_dT_obs, by = "id_plot")
  d[, step := k]; d[, label := inv_tmax$labels[k]]; d[, bit := bit]
  long[[k]] <- d
  ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
  x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
  fits[[k]] <- data.table(step = k, label = inv_tmax$labels[k],
                           r = suppressWarnings(cor(x, y)),
                           RMSE = sqrt(mean((y - x)^2)),
                           MAE  = mean(abs(y - x)),
                           bias = mean(y - x), n = length(x))
}
long <- rbindlist(long, fill = TRUE)
fits <- rbindlist(fits)
long$label <- factor(as.character(long$label), levels = inv_tmax$labels)

cli_alert("V10 inverse Tmax fits :"); print(fits)
fwrite(fits, file.path(OUT, "tab_v10_HOBO_inverse_Tmax.csv"))

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
  geom_point(colour = "#2C5F2D", alpha = 0.75, size = 3.0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 1.0, alpha = 0.22) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4,
            size = 7.5, fontface = "bold",
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8,
            size = 7, fontface = "bold",
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
            parse = TRUE, hjust = -0.08, vjust = 4.2,
            size = 7, fontface = "bold",
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ label, nrow = 1, labeller = label_parsed) +
  coord_fixed(xlim = xy_lim, ylim = xy_lim) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  theme_bw(base_size = 22) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 20),
        axis.title = element_text(face = "bold", size = 22),
        axis.text  = element_text(size = 16),
        # Pas de quadrillage — fond propre
        panel.grid       = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# 5 panneaux carrés alignés ; format 24 x 7 in @ 600 dpi PNG + PDF vector
ggsave(file.path(OUT, "fig_forward_selection_HOBO_v10_inverse_Tmax.png"),
        p, width = 24, height = 7, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_forward_selection_HOBO_v10_inverse_Tmax.pdf"),
        p, width = 24, height = 7, device = cairo_pdf)
cli_alert_success("Saved fig_forward_selection_HOBO_v10_inverse_Tmax.png (+ .pdf, square panels)")
