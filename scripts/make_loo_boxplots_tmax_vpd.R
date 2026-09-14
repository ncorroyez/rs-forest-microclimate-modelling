# ==============================================================================
# Per-plot LOO Δv boxplots on ΔTmax and ΔVPDmax at 1 m (nair==1), for the whole
# period and for the 10% hottest days. Companion / alternative to the heatmap.
#
#   Δv = metric_micro(REF) − metric_micro(LOO-v)   for v in {LAI,Hmax,fCover,LAD}
#   metric = daily max Tair (°C) or daily max VPD (kPa) at 1 m, then averaged
#   over the period. Macro term cancels in the LOO difference, so Δv on the
#   micro maximum == Δv on the offset ΔTmax / ΔVPDmax.
#   Δv < 0  => removing v raises the max  => v cools / dries-down => buffers.
#
# Hottest days = top 10% of macro daily Tmax (ERA5).
# Outputs: fig_loo_boxplot_tmax_vpd.{png,pdf} + tab_loo_deltav_tmax_vpd.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(tidyverse); library(lubridate); library(musica.tools)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT     <- here::here("outputs/figs_MEB2026_final")
P_HPA   <- 1013
REF_BIT  <- "1111"
LOO_BITS <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))

# ---- hottest days (top 10% of macro daily Tmax) -----------------------------
df_macro <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq))
thr_hot  <- quantile(df_macro$Tmax_macro, 0.90, na.rm = TRUE)
hot_days <- df_macro[Tmax_macro >= thr_hot, date]
cli_alert("Hottest-day threshold (macro Tmax) = {round(thr_hot,1)}°C ; {length(hot_days)}/{nrow(df_macro)} days")

# ---- per-plot daily Tmax & VPDmax at 1 m (nair==1) --------------------------
daily_1m <- function(nc_path) {
  nc <- try(nc_open(nc_path), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc))
  tu  <- ncatt_get(nc, "time", "units")$value
  t0  <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec<- floor_date(t0 + dhours(ncvar_get(nc, "time")), "hour")
  Tk  <- ncvar_get(nc, "Tair_z")[1, ]                 # nair==1 (lowest ~1 m)
  wmr <- ncvar_get(nc, "wair_z")[1, ]
  Tc  <- Tk - 273.15
  e   <- (wmr / (1 + wmr)) * P_HPA
  vpd <- pmax(esat_hpa(Tc) - e, 0) / 10               # kPa
  data.table(date = as.Date(tvec), Tc = Tc, vpd = vpd)[
    date %in% CFG$date_seq, .(Tmax = max(Tc, na.rm = TRUE),
                              VPDmax = max(vpd, na.rm = TRUE)), by = date]
}
collect_metric <- function(bit) {
  dir <- { v10 <- here::here("out_files/musica_hobo_v10_fcovmean", bit)
           v9  <- here::here("out_files/musica_hobo_v9", bit)
           if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9 }
  rows <- list()
  for (f in list.files(dir, pattern = "\\.nc$", full.names = TRUE)) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    d  <- daily_1m(f); if (is.null(d) || nrow(d) == 0) next
    rows[[id]] <- data.table(
      id_plot = id,
      Tmax_all = mean(d$Tmax, na.rm = TRUE),
      Tmax_hot = mean(d[date %in% hot_days, Tmax], na.rm = TRUE),
      VPD_all  = mean(d$VPDmax, na.rm = TRUE),
      VPD_hot  = mean(d[date %in% hot_days, VPDmax], na.rm = TRUE))
  }
  rbindlist(rows)
}

cli_h1("Collecting per-plot Tmax/VPDmax @1m for REF + 4 LOO coalitions")
S <- list(REF = collect_metric(REF_BIT))
for (v in names(LOO_BITS)) S[[v]] <- collect_metric(LOO_BITS[v])
cli_alert("REF plots: {nrow(S$REF)}")

# ---- Δv per plot (REF − LOO) for each metric × period -----------------------
mcols <- c("Tmax_all","Tmax_hot","VPD_all","VPD_hot")
dv <- rbindlist(lapply(names(LOO_BITS), function(v) {
  m <- merge(S$REF, S[[v]], by = "id_plot", suffixes = c("_ref","_loo"))
  out <- m[, .(id_plot, trait = v)]
  for (mc in mcols) out[[mc]] <- m[[paste0(mc,"_ref")]] - m[[paste0(mc,"_loo")]]
  out
}))
long <- melt(dv, id.vars = c("id_plot","trait"), measure.vars = mcols,
             variable.name = "key", value.name = "delta")
long[, metric := ifelse(grepl("^Tmax", key), "Delta*T[max]~(degree*C)",
                                              "Delta*VPD[max]~(kPa)")]
long[, period := ifelse(grepl("hot$", key), "10% hottest days", "All period")]
long[, period := factor(period, levels = c("All period","10% hottest days"))]

# cluster labels (for table)
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[, .(x,y,Cluster=as.character(Cluster))]
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
hxy <- sf::st_coordinates(hobo_pts)
clu <- data.table(id_plot = hobo_pts$id_plot, x = hxy[,"X"], y = hxy[,"Y"])
clu[, Cluster := relabel_cluster(sapply(seq_len(.N), function(i)
  df_forest$Cluster[which.min((df_forest$x-x[i])^2 + (df_forest$y-y[i])^2)]))]
long <- merge(long, clu[, .(id_plot, Cluster)], by = "id_plot")
fwrite(long, file.path(OUT, "tab_loo_deltav_tmax_vpd.csv"))

cli_h2("Median Δv per trait × metric × period")
print(long[, .(median = round(median(delta, na.rm=TRUE),3)), by=.(metric,period,trait)][order(metric,period,median)])

# ---- ordering & labels ------------------------------------------------------
library(patchwork)
trait_lab  <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
ord <- c("LAI","fCover","LAD","Hmax")
long[, trait_lab := factor(trait_lab[trait], levels = trait_lab[ord])]
trait_fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")

# robust y-limits per metric (zoom on the IQR story; sparse-plot amplification
# outliers extend far above and are clipped — see by-cluster view for those).
ylim_of <- function(d) {
  q <- quantile(d, c(0.02, 0.98), na.rm = TRUE)
  pad <- diff(q) * 0.12; c(q[1] - pad, q[2] + pad)
}
n_clip <- function(d, yl) sum(d < yl[1] | d > yl[2], na.rm = TRUE)

panel <- function(metric_key, ylab_expr) {
  d  <- long[metric == metric_key]
  yl <- ylim_of(d$delta)
  nc <- n_clip(d$delta, yl)
  ggplot(d, aes(trait_lab, delta, fill = trait)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.0, alpha = 0.35, colour = "grey20") +
    facet_wrap(~ period) +
    coord_cartesian(ylim = yl) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_fill_manual(values = trait_fill) +
    labs(x = NULL, y = ylab_expr,
         subtitle = if (nc > 0) sprintf("(%d sparse-plot amplification points clipped above)", nc) else NULL) +
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 10, colour = "grey45"),
          legend.position = "none")
}

pT <- panel("Delta*T[max]~(degree*C)", expression(Delta[v]^"LOO" ~ "on " * Delta*T[max] ~ "(" * degree * "C)"))
pV <- panel("Delta*VPD[max]~(kPa)",    expression(Delta[v]^"LOO" ~ "on " * Delta*VPD[max] ~ "(kPa)"))

p <- (pT / pV) +
  plot_annotation(
    caption = "1 m above ground (nair==1). Δv<0 = trait lowers the daily max (cools / dries-down) => buffers. fCover dominates and strengthens on hot days.",
    theme = theme(plot.caption = element_text(size = 11, colour = "grey40")))

ggsave(file.path(OUT, "fig_loo_boxplot_tmax_vpd.png"), p, width = 13, height = 10, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_loo_boxplot_tmax_vpd.pdf"), p, width = 13, height = 10, device = cairo_pdf)
cli_alert_success("Saved fig_loo_boxplot_tmax_vpd (png + pdf)")
