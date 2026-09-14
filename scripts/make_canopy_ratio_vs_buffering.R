# ==============================================================================
# Canopy Ratio (CR) vs microclimate buffering — answers Bouwen (2025) open
# question (Ch V.4.3, l.4607): "explore how slope values relate to CR and
# whether it offers additional explanatory power".
#
# CR definition (Starck et al. 2025, after Schneider et al. 2020 / Maeda 2022):
#   CR = (RH98 - RH25) / RH98
# where RHq = height at which fraction q of cumulative plant material (from the
# ground up) is reached. CR is the fraction of canopy height lying ABOVE RH25.
#   CR -> 1  => bottom-heavy (understory-rich, RH25 low)
#   CR -> 0  => top-heavy (overstory-concentrated, clear trunk, RH25 high)
#
# CR is NOT added to the Shapley/LOO lattice (it is derived from the LAD profile
# -> collinear with the LAD trait -> would reintroduce concurvity). Instead we
# test it correlationally and as an *added* predictor beyond LAI / VCI.
#
# Output : outputs/figs_MEB2026_final/fig_canopy_ratio_vs_buffering.{png,pdf}
#          + tab_canopy_ratio_metrics.csv  + tab_canopy_ratio_added_power.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/validation.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT   <- here::here("outputs/figs_MEB2026_final")
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

# ---- HOBO LiDAR inputs (LAD profiles + LAI/VCI/fCover/Hmax) ------------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                            CFG$ids_to_remove,
                                            hobo_buffer_mode = "buffer25"))
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
df_hobo$id_plot <- hobo_pts$id_plot

# ---- Canopy Ratio per plot --------------------------------------------------
lad_cols <- grep("^LAD_Layer_", names(df_hobo), value = TRUE)
z_lad    <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
# Height at which cumulative plant material (from ground) reaches fraction q.
rh_quantile <- function(d, z, q) {
  d[is.na(d)] <- 0
  tot <- sum(d); if (tot <= 0) return(NA_real_)
  cf  <- cumsum(d) / tot                 # fraction at/below top of each layer
  idx <- which(cf >= q)[1]
  if (is.na(idx)) return(max(z))
  if (idx == 1L) return(z[1] * q / cf[1])          # interp 0 -> z[1]
  z[idx - 1L] + (q - cf[idx - 1L]) / (cf[idx] - cf[idx - 1L]) * (z[idx] - z[idx - 1L])
}
canopy_ratio <- function(d) {
  rh25 <- rh_quantile(d, z_lad, 0.25)
  rh98 <- rh_quantile(d, z_lad, 0.98)
  if (is.na(rh25) || is.na(rh98) || rh98 <= 0) return(NA_real_)
  (rh98 - rh25) / rh98
}
df_hobo[, CR := sapply(seq_len(.N), function(i)
  canopy_ratio(as.numeric(.SD[i, ..lad_cols])))]

metrics <- df_hobo[, .(id_plot, LAI, VCI, fCover, Hmax, CR)]

# ---- Cluster assignment (P1..P4) -------------------------------------------
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[
  , .(x, y, Cluster = as.character(Cluster))]
hobo_xy <- sf::st_coordinates(hobo_pts)
df_hobo_clu <- data.table(id_plot = hobo_pts$id_plot,
                          x = hobo_xy[, "X"], y = hobo_xy[, "Y"])
df_hobo_clu[, Cluster := sapply(seq_len(.N), function(i) {
  d2 <- (df_forest$x - x[i])^2 + (df_forest$y - y[i])^2
  df_forest$Cluster[which.min(d2)]
})]
df_hobo_clu[, Cluster := relabel_cluster(Cluster)]

# ---- Observed HOBO log(slope) per plot --------------------------------------
hobo_raw <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a",
         as.Date(datetime) %in% CFG$date_seq,
         !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(id_plot, time) %>%
  summarise(t_hobo = mean(t_hobo, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()
hobo_obs <- merge(hobo_raw, era5_hourly, by = "time")[, {
  fit <- lm(t_hobo ~ Tair_era5)
  .(log_slope_obs = log(as.numeric(coef(fit)[2])))
}, by = id_plot]

# ---- Simulated REF (1111) log(slope) per plot -------------------------------
ref_dir <- {
  v10 <- here::here("out_files/musica_hobo_v10_fcovmean/1111")
  v9  <- here::here("out_files/musica_hobo_v9/1111")
  if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9
}
sim_rows <- list()
for (f in list.files(ref_dir, pattern = "\\.nc$", full.names = TRUE)) {
  id  <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
  res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) next
  sim_rows[[id]] <- data.table(id_plot = id, log_slope_sim = log(res$slope))
}
sim_ref <- rbindlist(sim_rows)

# ---- Merge ------------------------------------------------------------------
D <- Reduce(function(a, b) merge(a, b, by = "id_plot"),
            list(metrics, df_hobo_clu[, .(id_plot, Cluster)], hobo_obs, sim_ref))
cli_alert_success("Merged {nrow(D)} plots ; CR range [{round(min(D$CR,na.rm=T),2)}, {round(max(D$CR,na.rm=T),2)}]")

# ---- Correlations & added explanatory power ---------------------------------
cor_tab <- rbindlist(lapply(c("LAI", "VCI", "fCover", "Hmax", "CR"), function(v) {
  x <- D[[v]]; y <- D$log_slope_obs; ok <- is.finite(x) & is.finite(y)
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.table(metric = v, rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}))
cli_h2("Spearman |rho| with observed log(slope)")
print(cor_tab[order(-abs(rho))])

# Nested linear models: does CR add beyond LAI? beyond VCI? beyond LAI+VCI?
addpow <- function(base_vars) {
  d <- D[is.finite(CR) & is.finite(log_slope_obs)]
  f0 <- as.formula(paste("log_slope_obs ~", paste(base_vars, collapse = " + ")))
  f1 <- as.formula(paste("log_slope_obs ~", paste(c(base_vars, "CR"), collapse = " + ")))
  m0 <- lm(f0, d); m1 <- lm(f1, d)
  an <- anova(m0, m1)
  data.table(base = paste(base_vars, collapse = "+"),
             R2_base = summary(m0)$r.squared,
             R2_plusCR = summary(m1)$r.squared,
             dR2 = summary(m1)$r.squared - summary(m0)$r.squared,
             F = an$F[2], p_CR = an$`Pr(>F)`[2])
}
added <- rbindlist(lapply(list("LAI", "VCI", c("LAI","VCI")), addpow))
cli_h2("Added explanatory power of CR (nested-model F test)")
print(added)

fwrite(cor_tab, file.path(OUT, "tab_canopy_ratio_metrics.csv"))
fwrite(added,   file.path(OUT, "tab_canopy_ratio_added_power.csv"))
fwrite(D,       file.path(OUT, "tab_canopy_ratio_perplot.csv"))

# ---- Figure : CR vs buffering (obs + sim), coloured by cluster --------------
cr_rho <- cor_tab[metric == "CR"]
lng <- rbind(
  D[, .(id_plot, Cluster, CR, value = log_slope_obs, panel = "Observed (HOBO)")],
  D[, .(id_plot, Cluster, CR, value = log_slope_sim, panel = "Simulated (MuSICA REF)")])
lng[, panel := factor(panel, levels = c("Observed (HOBO)", "Simulated (MuSICA REF)"))]

ann <- lng[, {
  ok <- is.finite(CR) & is.finite(value)
  ct <- suppressWarnings(cor.test(CR[ok], value[ok], method = "spearman", exact = FALSE))
  .(lab = sprintf("rho==%+.2f~~(italic(p)==%.3f)~~n==%d",
                  unname(ct$estimate), ct$p.value, sum(ok)))
}, by = panel]

pal <- c(P1 = "#D7191C", P2 = "#FDAE61", P3 = "#74C476", P4 = "#1A9850")
p <- ggplot(lng, aes(CR, value)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#FFB400", fill = "#FFB400", alpha = 0.2, aes(group = 1)) +
  geom_point(aes(colour = Cluster, shape = Cluster), size = 3.2, alpha = 0.9) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), parse = TRUE,
            hjust = 1.05, vjust = 1.6, size = 5, inherit.aes = FALSE) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = pal) +
  labs(x = expression("Canopy Ratio  CR = (" * RH[98] - RH[25] * ") / " * RH[98]),
       y = expression(log(slope[micro/macro])),
       caption = "log(slope)<0 = buffering. High CR = understory-rich (RH25 low). Negative rho => understory-rich = stronger buffering.") +
  theme_bw(base_size = 18) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 11, colour = "grey40"),
        legend.position = "right")

ggsave(file.path(OUT, "fig_canopy_ratio_vs_buffering.png"),
       p, width = 13, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_canopy_ratio_vs_buffering.pdf"),
       p, width = 13, height = 6.2, device = cairo_pdf)
cli_alert_success("Saved fig_canopy_ratio_vs_buffering (png + pdf)")
