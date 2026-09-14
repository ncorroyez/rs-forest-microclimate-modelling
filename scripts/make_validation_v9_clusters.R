# ==============================================================================
# V9 HOBO validation figures with cluster context :
#   - Fig A : daily ΔTmax sim vs obs, all days, faceted Normal/Heatwave,
#             coloured/shape by cluster
#   - Fig B : daily ΔTmax sim vs obs, faceted by cluster (P1..P4)
#   - Fig C : per-plot mean ΔTmax sim vs obs (compact, for slide)
#
# Uses V9 REF (bit 1111) sims with the LEGACY MuSICA binary.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
  library(sf); library(terra)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/cluster_relabel.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cli_h1("V9 validation figures with cluster context")

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---------------------------------------------------------------------------
# 1. HOBO observations (daily Delta_obs)
# ---------------------------------------------------------------------------
df_hobo_obs <- as.data.table(read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq, df_macro,
                                              CFG$ids_to_remove))
cli_alert("HOBO obs daily rows : {nrow(df_hobo_obs)} | "
          ~ "{length(unique(df_hobo_obs$id_plot))} plots")

# ---------------------------------------------------------------------------
# 2. V9 sims daily Delta_sim from out_files/musica_hobo_v9/1111/
# ---------------------------------------------------------------------------
sim_dir <- here::here("out_files/musica_hobo_v9/1111")
nc_files <- list.files(sim_dir, pattern = "musica_out_HOBO_.+\\.nc$", full.names = TRUE)
cli_alert("V9 REF NCs : {length(nc_files)}")

sim_rows <- list()
for (f in nc_files) {
  id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
  res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                        z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) next
  sim_rows[[id]] <- as.data.table(res)[, .(id_plot = id, date, Tmax_macro,
                                            Tmax_sim = Tmax_micro,
                                            Delta_sim = Delta_Tmax)]
}
DT_sim <- rbindlist(sim_rows)
cli_alert("V9 sim daily rows : {nrow(DT_sim)} | "
          ~ "Delta_sim mean={round(mean(DT_sim$Delta_sim, na.rm=TRUE), 2)} °C")

# ---------------------------------------------------------------------------
# 3. Cluster assignment per HOBO plot via nearest forest cell
# ---------------------------------------------------------------------------
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
cli_alert("HOBO → cluster :")
print(table(df_hobo_clu$Cluster))

# ---------------------------------------------------------------------------
# 4. Merge sim+obs+cluster, build long daily table
# ---------------------------------------------------------------------------
DT <- merge(DT_sim[, .(id_plot, date, Delta_sim, Tmax_macro, Tmax_sim)],
            df_hobo_obs[, .(id_plot = id_plot, date, Delta_obs)],
            by = c("id_plot", "date"))
DT <- merge(DT, df_hobo_clu, by = "id_plot")
# Cluster from clhs_sample_floor05_v2.rds = RAW k-means code; relabel to P1..P4 by
# ascending LAI (raw 3->P1 sparsest .. 1->P4 densest). Direct paste0("P",.) was WRONG.
DT[, Cluster := relabel_cluster(Cluster)]
DT[, regime  := factor(ifelse(Tmax_macro >= 30, "Heatwave (Tmax ≥ 30 °C)", "Normal"),
                        levels = c("Normal", "Heatwave (Tmax ≥ 30 °C)"))]
cli_alert("Merged rows : {nrow(DT)}")
cli_alert("Per-cluster rows :"); print(DT[, .N, by = .(Cluster, regime)])

# ---------------------------------------------------------------------------
# Fig A : daily ΔTmax sim vs obs, faceted Normal/Heatwave, color by cluster
# ---------------------------------------------------------------------------
cli_h2("Fig A — daily by regime (Normal vs Heatwave)")
lims_A <- range(c(DT$Delta_obs, DT$Delta_sim), na.rm = TRUE)
pad_A  <- diff(lims_A) * 0.05; lims_A <- lims_A + c(-pad_A, pad_A)

stats_A <- DT[, {
  ok <- !is.na(Delta_obs) & !is.na(Delta_sim)
  x <- Delta_obs[ok]; y <- Delta_sim[ok]
  r <- suppressWarnings(cor(x, y))
  .(r = r, RMSE = sqrt(mean((y - x)^2)), MAE = mean(abs(y - x)),
    bias = mean(y - x), N = length(x))
}, by = regime]
print(stats_A)

ann_A <- stats_A[, .(regime,
                      lab = sprintf("italic(r)==%+.2f~~RMSE==%.2f~degree*C~~bias==%+.2f", r, RMSE, bias))]

p_A <- ggplot(DT, aes(x = Delta_obs, y = Delta_sim,
                       colour = Cluster, shape = Cluster)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.0, aes(group = 1), colour = "grey25") +
  geom_text(data = ann_A, aes(x = -Inf, y = Inf, label = lab),
            parse = TRUE, hjust = -0.05, vjust = 1.5, size = 4.5,
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ regime) +
  coord_fixed(xlim = lims_A, ylim = lims_A) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A")) +
  scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)"),
       colour = "Cluster", shape = "Cluster") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_HOBO_validation_v9_daily_byregime.png"),
        p_A, width = 12, height = 6, dpi = 300)
cli_alert_success("Saved fig_HOBO_validation_v9_daily_byregime.png")

# ---------------------------------------------------------------------------
# Fig B : daily ΔTmax sim vs obs, faceted by cluster
# ---------------------------------------------------------------------------
cli_h2("Fig B — daily by cluster")
stats_B <- DT[, {
  ok <- !is.na(Delta_obs) & !is.na(Delta_sim)
  x <- Delta_obs[ok]; y <- Delta_sim[ok]
  r <- suppressWarnings(cor(x, y))
  .(r = r, RMSE = sqrt(mean((y - x)^2)), bias = mean(y - x), N = length(x))
}, by = Cluster]
print(stats_B)
ann_B <- stats_B[, .(Cluster,
                      lab = sprintf("italic(r)==%+.2f~~RMSE==%.2f", r, RMSE))]

p_B <- ggplot(DT, aes(x = Delta_obs, y = Delta_sim, colour = Cluster, shape = Cluster)) +
  geom_point(alpha = 0.5, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, colour = "grey25") +
  geom_text(data = ann_B, aes(x = -Inf, y = Inf, label = lab),
            parse = TRUE, hjust = -0.05, vjust = 1.5, size = 4.2,
            inherit.aes = FALSE, colour = "grey15") +
  facet_wrap(~ Cluster, nrow = 1) +
  coord_fixed(xlim = lims_A, ylim = lims_A) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A")) +
  scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)")) +
  guides(colour = "none", shape = "none") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_HOBO_validation_v9_daily_bycluster.png"),
        p_B, width = 14, height = 4, dpi = 300)
cli_alert_success("Saved fig_HOBO_validation_v9_daily_bycluster.png")

# ---------------------------------------------------------------------------
# Fig C : per-plot mean ΔTmax sim vs obs (slide-ready)
# ---------------------------------------------------------------------------
cli_h2("Fig C — per-plot mean (slide-ready)")
DT_mean <- DT[, .(Delta_sim_mean = mean(Delta_sim, na.rm = TRUE),
                  Delta_obs_mean = mean(Delta_obs, na.rm = TRUE)),
              by = .(id_plot, Cluster)]
fit_C <- lm(Delta_sim_mean ~ Delta_obs_mean, data = DT_mean)
stats_C <- data.table(
  r    = cor(DT_mean$Delta_obs_mean, DT_mean$Delta_sim_mean),
  RMSE = sqrt(mean((DT_mean$Delta_sim_mean - DT_mean$Delta_obs_mean)^2)),
  MAE  = mean(abs(DT_mean$Delta_sim_mean - DT_mean$Delta_obs_mean)),
  bias = mean(DT_mean$Delta_sim_mean - DT_mean$Delta_obs_mean))
cli_alert("Per-plot fit (REF, daily mean) :"); print(stats_C)
lims_C <- range(c(DT_mean$Delta_obs_mean, DT_mean$Delta_sim_mean), na.rm = TRUE)
pad_C  <- diff(lims_C) * 0.08; lims_C <- lims_C + c(-pad_C, pad_C)

p_C <- ggplot(DT_mean, aes(x = Delta_obs_mean, y = Delta_sim_mean,
                            colour = Cluster, shape = Cluster)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
  geom_point(size = 4, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, aes(group = 1),
              colour = "#FFB400", fill = "#FFB400", alpha = 0.22) +
  annotate("text", x = -Inf, y = Inf,
            label = sprintf("italic(r) == %+.2f", stats_C$r),
            parse = TRUE, hjust = -0.1, vjust = 1.4, size = 6, colour = "grey15") +
  annotate("text", x = -Inf, y = Inf,
            label = sprintf("RMSE == %.2f * ' °C'", stats_C$RMSE),
            parse = TRUE, hjust = -0.1, vjust = 2.8, size = 6, colour = "grey15") +
  annotate("text", x = -Inf, y = Inf,
            label = sprintf("bias == %+.2f * ' °C'", stats_C$bias),
            parse = TRUE, hjust = -0.1, vjust = 4.2, size = 6, colour = "grey15") +
  coord_fixed(xlim = lims_C, ylim = lims_C) +
  scale_colour_manual(values = c(P1 = "#1B9E77", P2 = "#D95F02",
                                  P3 = "#7570B3", P4 = "#E7298A")) +
  scale_shape_manual(values = c(P1 = 16, P2 = 17, P3 = 15, P4 = 18)) +
  labs(x = expression("Field " ~ Delta * T[max] ~ "(°C)"),
       y = expression("Simulated " ~ Delta * T[max] ~ "(°C)"),
       colour = "Cluster", shape = "Cluster") +
  theme_bw(base_size = 15) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_HOBO_validation_v9_perplot_mean.png"),
        p_C, width = 8, height = 7, dpi = 300)
cli_alert_success("Saved fig_HOBO_validation_v9_perplot_mean.png")

# Save merged data for downstream use
saveRDS(list(DT = DT, DT_mean = DT_mean,
              stats_regime = stats_A, stats_cluster = stats_B,
              stats_perplot = stats_C, df_hobo_clu = df_hobo_clu),
         here::here("outputs/v9_validation_clusters.rds"))
cli_alert_success("Saved outputs/v9_validation_clusters.rds")
