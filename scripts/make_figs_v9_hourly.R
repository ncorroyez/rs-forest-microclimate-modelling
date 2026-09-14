# ==============================================================================
# V9 figures with HOURLY-SLOPE metric :
#   - Fig 1 : LOO heatmap (archetypes × LOO_v) — Δ slope = REF - LOO_v
#             (negative Δ ⇒ removing v moves slope toward 1 ⇒ v was buffering)
#   - Fig 2 : HOBO forward-selection (simulated slope vs observed slope per plot)
#   - Tab  : buf/amp regime per coalition (slope < 1 = buffer)
#
# Reads the V9 sims produced by scripts/run_v9_legacy_binary.R.
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
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

state <- readRDS(here::here("outputs/v9_pipeline_state.rds"))
forward_order <- state$forward_order
bits_seq      <- state$bits_seq
DT_arch_w     <- state$arch_slopes

era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

# ---------------------------------------------------------------------------
# Fig 1 : Heatmap (LOO Δ in °C — daily Tmax_micro, like V2)
# ---------------------------------------------------------------------------
cli_h1("Fig 1 — LOO heatmap (ΔTmax in °C, daily)")

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

extract_tmax_one <- function(nc_path) {
  if (!file.exists(nc_path) || file.size(nc_path) < 1e6) return(NA_real_)
  res <- tryCatch(extract_deltatmax_one(nc_path, df_macro, CFG$date_seq,
                                        z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

# Archetype clusters from saved state (preserve order)
archetypes_clusters <- DT_arch_w$Cluster

arch_root <- here::here("out_files/H1_archetypes_v9")
arch_rows <- list()
for (cl in archetypes_clusters) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  for (sc_nm in names(.ARCH_COAL_MAP)) {
    bit <- .ARCH_COAL_MAP[sc_nm]
    nc <- file.path(arch_root, arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    arch_rows[[length(arch_rows) + 1L]] <- data.table(
      Cluster = cl, bit_code = bit, Tmax_micro = extract_tmax_one(nc))
  }
}
DT_arch_tmax <- rbindlist(arch_rows)
DT_arch_tmax_w <- dcast(DT_arch_tmax, Cluster ~ bit_code, value.var = "Tmax_micro")

VARS <- c("LAI","Hmax","fCover","LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

delta_long <- list()
for (v in VARS) {
  for (i in seq_len(nrow(DT_arch_tmax_w))) {
    delta_long[[length(delta_long) + 1L]] <- data.table(
      Cluster  = DT_arch_tmax_w$Cluster[i],
      variable = v,
      delta    = DT_arch_tmax_w[[ "1111" ]][i] - DT_arch_tmax_w[[ loo_bits[v] ]][i])
  }
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
long_full[, variable := factor(variable, levels = VARS)]
long_full[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]

# Order rows by the forward_order (strongest buffer at top)
var_lab_order <- ifelse(forward_order == "Hmax", "H[max]", forward_order)
long_full[, var_lab := factor(var_lab, levels = var_lab_order)]

max_abs <- max(abs(long_full$delta), na.rm = TRUE) * 1.05
sep_x <- 4.5

p_hm <- ggplot(long_full, aes(x = Profile, y = var_lab, fill = delta)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.2f", delta)),
            fontface = "bold", size = 5,
            colour = ifelse(abs(long_full$delta) > max_abs * 0.55,
                            "white", "grey15")) +
  geom_vline(xintercept = sep_x, colour = "white", linewidth = 4) +
  geom_vline(xintercept = sep_x, colour = "grey40", linewidth = 1.0) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-max_abs, max_abs),
                       name = expression(Delta[v]^"LOO" ~ "(°C)")) +
  scale_y_discrete(limits = rev,
                   labels = function(x) parse(text = as.character(x))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        legend.position = "right")
fig_hm <- file.path(OUT, "fig_heatmap_LOO_v9_hourly.png")
ggsave(fig_hm, p_hm, width = 12, height = 5, dpi = 300)
cli_alert_success("Saved {.path {fig_hm}}")

# ---------------------------------------------------------------------------
# Helper : extract per-plot slope from HOBO sims of a given coalition
# ---------------------------------------------------------------------------
collect_step_slope <- function(bit) {
  subdir <- here::here(file.path("out_files/musica_hobo_v9", bit))
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(
      extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                               z_target = CFG$tair_target_height),
      error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- data.table(id_plot = id, slope_sim = res$slope,
                             dT_sim = res$dT_mean)
  }
  rbindlist(rows)
}

# ---------------------------------------------------------------------------
# HOBO observed slope (Tair_hobo ~ Tair_era5 hourly)
# ---------------------------------------------------------------------------
cli_h1("Computing HOBO observed slope per plot")

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
hobo_slopes <- hobo_with_era5[, {
  fit <- lm(t_hobo ~ Tair_era5)
  sl <- coef(fit)[2]
  .(slope_obs = as.numeric(sl), dT_obs = mean(t_hobo - Tair_era5, na.rm = TRUE),
    n = .N)
}, by = id_plot]
cli_alert("HOBO observed slopes : n={nrow(hobo_slopes)}, "
          ~ "buffer(slope<1)={sum(hobo_slopes$slope_obs < 1)} / amp={sum(hobo_slopes$slope_obs > 1)}")

# ---------------------------------------------------------------------------
# Fig 2 : HOBO forward selection (slope_sim ~ slope_obs across coalitions)
# ---------------------------------------------------------------------------
cli_h1("Fig 2 — HOBO forward selection (slope)")

build_label <- function(k) {
  if (k == 1) return("Baseline")
  cum_vars <- forward_order[1:(k-1)]
  cum_vars_disp <- ifelse(cum_vars == "Hmax", "H_max", cum_vars)
  paste("+", paste(cum_vars_disp, collapse = " + "))
}
labels <- vapply(seq_len(5), build_label, character(1))

long_fwd <- list()
buf_amp_summary <- list()
for (k in seq_along(bits_seq)) {
  df_k <- collect_step_slope(bits_seq[k])
  if (nrow(df_k) == 0) next
  df_k <- merge(df_k, hobo_slopes[, .(id_plot, slope_obs, dT_obs)], by = "id_plot")
  df_k[, step := k]; df_k[, label := labels[k]]; df_k[, bit := bits_seq[k]]
  long_fwd[[k]] <- df_k
  buf_amp_summary[[k]] <- data.table(
    step = k, label = labels[k], bit = bits_seq[k],
    n = nrow(df_k),
    n_buffer = sum(df_k$slope_sim < 1),
    n_amp    = sum(df_k$slope_sim > 1))
}
long_fwd <- rbindlist(long_fwd, fill = TRUE)
long_fwd$label <- factor(as.character(long_fwd$label), levels = labels)
buf_amp_summary <- rbindlist(buf_amp_summary)
cli_alert("buf/amp summary :")
print(buf_amp_summary)
fwrite(buf_amp_summary, file.path(OUT, "tab_v9_buf_amp_per_coalition.csv"))

fits <- long_fwd[, {
  ok <- !is.na(slope_sim) & !is.na(slope_obs)
  x <- slope_obs[ok]; y <- slope_sim[ok]
  r <- suppressWarnings(cor(x, y))
  rmse <- sqrt(mean((y - x)^2))
  mae  <- mean(abs(y - x))
  .(n = length(x), r = r, RMSE = rmse, MAE = mae)
}, by = .(step, label)]
cli_alert("Per-step fits (slope_sim vs slope_obs) :"); print(fits)
fwrite(fits, file.path(OUT, "tab_v9_forward_selection_fits_slope.csv"))

xy_lim <- range(c(long_fwd$slope_sim, long_fwd$slope_obs), na.rm = TRUE)
xy_pad <- diff(xy_lim) * 0.06
xy_lim <- c(xy_lim[1] - xy_pad, xy_lim[2] + xy_pad)

ann <- fits[, .(label,
                lab_r    = sprintf("italic(r) == %+.2f", r),
                lab_rmse = sprintf("RMSE == %.3f", RMSE),
                lab_mae  = sprintf("MAE == %.3f", MAE))]

p_fwd <- ggplot(long_fwd, aes(x = slope_obs, y = slope_sim)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey80") +
  geom_point(colour = "#2C5F2D", alpha = 0.7, size = 2.0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400",
              linewidth = 0.9, alpha = 0.22) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r),
            parse = TRUE, hjust = -0.08, vjust = 1.4, size = 6,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
            parse = TRUE, hjust = -0.08, vjust = 2.8, size = 6,
            inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
            parse = TRUE, hjust = -0.08, vjust = 4.2, size = 6,
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ label, nrow = 1) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  labs(x = "Field slope (Tair_HOBO ~ Tair_ERA5, hourly)",
       y = "Simulated slope (Tair_sim ~ Tair_ERA5, hourly)") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())
fig_fwd <- file.path(OUT, "fig_forward_selection_HOBO_v9_hourly.png")
ggsave(fig_fwd, p_fwd, width = 18, height = 5, dpi = 300)
cli_alert_success("Saved {.path {fig_fwd}}")

cli_h1("V9 figures complete.")
