# ==============================================================================
# Per-plot LOO Δv boxplots (53 HOBO plots) — companion to the LOO heatmap.
# Δv_plot = log(slope_REF) − log(slope_LOO-v)  for v in {LAI,Hmax,fCover,LAD}.
#   Δv < 0  => removing v reduces buffering  => v CONTRIBUTES buffering.
# Shows the distribution (median + spread) the heatmap's single archetype value hides.
#
# Outputs:
#   fig_loo_boxplot_pooled.{png,pdf}        (1 box / trait, 53 plots)
#   fig_loo_boxplot_by_cluster.{png,pdf}    (faceted P1→P4)
#   tab_loo_deltav_perplot.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(tidyverse); library(lubridate); library(musica.tools)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

REF_BIT  <- "1111"
LOO_BITS <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

# ---- per-plot log(slope) for one coalition ---------------------------------
slope_dir <- function(bit) {
  v10 <- here::here("out_files/musica_hobo_v10_fcovmean", bit)
  v9  <- here::here("out_files/musica_hobo_v9", bit)
  if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9
}
collect_logslope <- function(bit) {
  rows <- list()
  for (f in list.files(slope_dir(bit), pattern = "\\.nc$", full.names = TRUE)) {
    id  <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                              z_target = CFG$tair_target_height),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    rows[[id]] <- data.table(id_plot = id, log_slope = log(res$slope))
  }
  rbindlist(rows)
}

cli_h1("Collecting per-plot log(slope) for REF + 4 LOO coalitions")
S <- list(REF = collect_logslope(REF_BIT))
for (v in names(LOO_BITS)) S[[v]] <- collect_logslope(LOO_BITS[v])
cli_alert("REF plots: {nrow(S$REF)}")

# ---- Δv per plot ------------------------------------------------------------
dv <- rbindlist(lapply(names(LOO_BITS), function(v) {
  m <- merge(S$REF[, .(id_plot, ref = log_slope)],
             S[[v]][, .(id_plot, loo = log_slope)], by = "id_plot")
  m[, .(id_plot, trait = v, delta = ref - loo)]
}))

# ---- cluster labels ---------------------------------------------------------
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[
  , .(x, y, Cluster = as.character(Cluster))]
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
hxy <- sf::st_coordinates(hobo_pts)
clu <- data.table(id_plot = hobo_pts$id_plot, x = hxy[, "X"], y = hxy[, "Y"])
clu[, Cluster := sapply(seq_len(.N), function(i) {
  d2 <- (df_forest$x - x[i])^2 + (df_forest$y - y[i])^2
  df_forest$Cluster[which.min(d2)] })]
clu[, Cluster := relabel_cluster(Cluster)]
dv <- merge(dv, clu[, .(id_plot, Cluster)], by = "id_plot")

# trait order: strongest buffering contributor (most negative median) first
ord <- dv[, .(m = median(delta, na.rm = TRUE)), by = trait][order(m), trait]
trait_lab <- c(LAI = "LAI", Hmax = "italic(H)[max]", fCover = "fCover", LAD = "LAD")
dv[, trait := factor(trait, levels = ord)]
dv[, trait_lab := factor(trait_lab[as.character(trait)], levels = trait_lab[ord])]

fwrite(dv, file.path(OUT, "tab_loo_deltav_perplot.csv"))
cli_h2("Pooled median Δv per trait"); print(dv[, .(median = round(median(delta),3),
  IQR = round(IQR(delta),3), n = .N), by = trait][order(median)])

pal <- c(P1 = "#D7191C", P2 = "#FDAE61", P3 = "#74C476", P4 = "#1A9850")
trait_fill <- c(LAI = "#4575B4", Hmax = "#91BFDB", fCover = "#FC8D59", LAD = "#D73027")

base_theme <- theme_bw(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", size = 20),
        plot.caption = element_text(size = 12, colour = "grey40"),
        legend.position = "none")

# ---- (A) pooled -------------------------------------------------------------
pA <- ggplot(dv, aes(trait_lab, delta, fill = as.character(trait))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.45, colour = "grey20") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_fill_manual(values = trait_fill) +
  labs(x = NULL, y = expression(Delta[v]^"LOO" ~ "(log " * slope * ")"),
       title = "Per-plot LOO attribution (53 HOBO plots)",
       caption = expression(Delta[v] < 0 ~ "= trait contributes buffering")) +
  base_theme
ggsave(file.path(OUT, "fig_loo_boxplot_pooled.png"), pA, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_loo_boxplot_pooled.pdf"), pA, width = 10, height = 7, device = cairo_pdf)
cli_alert_success("Saved fig_loo_boxplot_pooled")

# ---- (B) faceted by cluster -------------------------------------------------
pB <- ggplot(dv, aes(trait_lab, delta, fill = as.character(trait))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.5, colour = "grey20") +
  facet_wrap(~ Cluster, nrow = 1) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_fill_manual(values = trait_fill) +
  labs(x = NULL, y = expression(Delta[v]^"LOO" ~ "(log " * slope * ")"),
       caption = expression(Delta[v] < 0 ~ "= trait contributes buffering. P1 sparse → P4 dense.")) +
  base_theme +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 15, angle = 30, hjust = 1))
ggsave(file.path(OUT, "fig_loo_boxplot_by_cluster.png"), pB, width = 16, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_loo_boxplot_by_cluster.pdf"), pB, width = 16, height = 7, device = cairo_pdf)
cli_alert_success("Saved fig_loo_boxplot_by_cluster")
