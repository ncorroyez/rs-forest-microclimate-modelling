# ==============================================================================
# Sampling-height comparison: LOO Δv on ΔTmax & ΔVPDmax at nair==1 (legacy node,
# 0.0285*Hmax => variable 0.1-1.1 m) VS interpolated fixed 1 m (HOBO height).
# Tests whether the LAI sparse-amplification / dense-buffering sign-flip is a
# structural effect or partly a sampling-height artifact.
#
#   Δv = metric_micro(REF) − metric_micro(LOO-v)   v in {LAI,Hmax,fCover,LAD}
#   Δv < 0  => trait lowers the daily max => buffers.
# Outputs: fig_loo_height_compare_tmax.{png,pdf}, fig_loo_height_compare_vpd.{png,pdf}
#          tab_loo_deltav_height_compare.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(tidyverse); library(lubridate)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT      <- here::here("outputs/figs_MEB2026_final")
P_HPA    <- 1013
Z_FIX    <- 1.0
REF_BIT  <- "1111"
LOO_BITS <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110")
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))

df_macro <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq))
thr_hot  <- quantile(df_macro$Tmax_macro, 0.90, na.rm = TRUE)
hot_days <- df_macro[Tmax_macro >= thr_hot, date]

# ---- hourly Tair & wair at both heights, then daily max ---------------------
daily_both <- function(nc_path) {
  nc <- try(nc_open(nc_path), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc))
  tu  <- ncatt_get(nc, "time", "units")$value
  t0  <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec<- floor_date(t0 + dhours(ncvar_get(nc, "time")), "hour")
  Tk  <- ncvar_get(nc, "Tair_z")          # nair x time
  wmr <- ncvar_get(nc, "wair_z")
  rh  <- ncvar_get(nc, "relative_height")
  vh  <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl  <- rh * vh                          # absolute layer heights (m)

  vpd_of <- function(Tc, w) pmax(esat_hpa(Tc) - (w/(1+w))*P_HPA, 0) / 10  # kPa

  # nair==1
  Tc_n1 <- Tk[1, ] - 273.15;  vpd_n1 <- vpd_of(Tc_n1, wmr[1, ])
  # interp to fixed 1 m
  if (Z_FIX <= zl[1]) { ilo <- 1L; ihi <- 1L; w <- 0 }
  else if (Z_FIX >= zl[length(zl)]) { ilo <- length(zl); ihi <- ilo; w <- 0 }
  else { ilo <- max(which(zl <= Z_FIX)); ihi <- ilo + 1L
         w <- (Z_FIX - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  Tc_z1 <- ((1-w)*Tk[ilo, ] + w*Tk[ihi, ]) - 273.15
  w_z1  <- (1-w)*wmr[ilo, ] + w*wmr[ihi, ]
  vpd_z1<- vpd_of(Tc_z1, w_z1)

  dt <- data.table(date = as.Date(tvec),
                   Tn1 = Tc_n1, Vn1 = vpd_n1, Tz1 = Tc_z1, Vz1 = vpd_z1)[
                   date %in% CFG$date_seq]
  dt[, .(Tmax_n1 = max(Tn1), VPD_n1 = max(Vn1),
         Tmax_z1 = max(Tz1), VPD_z1 = max(Vz1)), by = date]
}

collect <- function(bit) {
  dir <- { v10 <- here::here("out_files/musica_hobo_v10_fcovmean", bit)
           v9  <- here::here("out_files/musica_hobo_v9", bit)
           if (dir.exists(v10) && length(list.files(v10, "\\.nc$")) >= 53) v10 else v9 }
  rows <- list()
  for (f in list.files(dir, pattern = "\\.nc$", full.names = TRUE)) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    d  <- daily_both(f); if (is.null(d) || nrow(d) == 0) next
    agg <- function(col, days = NULL) mean(if (is.null(days)) d[[col]] else d[date %in% days][[col]], na.rm = TRUE)
    rows[[id]] <- data.table(
      id_plot = id,
      Tmax_n1_all = agg("Tmax_n1"), Tmax_n1_hot = agg("Tmax_n1", hot_days),
      Tmax_z1_all = agg("Tmax_z1"), Tmax_z1_hot = agg("Tmax_z1", hot_days),
      VPD_n1_all  = agg("VPD_n1"),  VPD_n1_hot  = agg("VPD_n1", hot_days),
      VPD_z1_all  = agg("VPD_z1"),  VPD_z1_hot  = agg("VPD_z1", hot_days))
  }
  rbindlist(rows)
}

cli_h1("Collecting both heights for REF + 4 LOO")
S <- list(REF = collect(REF_BIT)); for (v in names(LOO_BITS)) S[[v]] <- collect(LOO_BITS[v])
cli_alert("REF plots: {nrow(S$REF)}")

mcols <- setdiff(names(S$REF), "id_plot")
dv <- rbindlist(lapply(names(LOO_BITS), function(v) {
  m <- merge(S$REF, S[[v]], by = "id_plot", suffixes = c("_ref","_loo"))
  out <- m[, .(id_plot, trait = v)]
  for (mc in mcols) out[[mc]] <- m[[paste0(mc,"_ref")]] - m[[paste0(mc,"_loo")]]
  out
}))
long <- melt(dv, id.vars = c("id_plot","trait"), variable.name = "key", value.name = "delta")
long[, metric := fifelse(grepl("^Tmax", key), "Delta*T[max]~(degree*C)", "Delta*VPD[max]~(kPa)")]
long[, height := fifelse(grepl("_n1_", key), "nair==1 (0.1-1.1 m)", "1 m (interp., HOBO)")]
long[, period := fifelse(grepl("_hot$", key), "10% hottest", "All period")]
long[, height := factor(height, levels = c("nair==1 (0.1-1.1 m)", "1 m (interp., HOBO)"))]

# cluster
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[, .(x,y,Cluster=as.character(Cluster))]
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
hxy <- sf::st_coordinates(hobo_pts)
clu <- data.table(id_plot = hobo_pts$id_plot, x = hxy[,"X"], y = hxy[,"Y"])
clu[, Cluster := relabel_cluster(sapply(seq_len(.N), function(i)
  df_forest$Cluster[which.min((df_forest$x-x[i])^2 + (df_forest$y-y[i])^2)]))]
long <- merge(long, clu[, .(id_plot, Cluster)], by = "id_plot")
long[, Cluster := factor(Cluster, levels = paste0("P",1:4))]
fwrite(long, file.path(OUT, "tab_loo_deltav_height_compare.csv"))

# ---- comparison figures (all-period; height × cluster) ----------------------
trait_lab <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
ord <- c("LAI","fCover","LAD","Hmax")
long[, trait_lab := factor(trait_lab[trait], levels = trait_lab[ord])]
trait_fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")

make_cmp <- function(metric_key, ylab_expr, fname) {
  d <- long[metric == metric_key & period == "All period"]
  q <- quantile(d$delta, c(0.04, 0.96), na.rm = TRUE); pad <- diff(q)*0.10
  yl <- c(q[1]-pad, q[2]+pad); nc <- sum(d$delta < yl[1] | d$delta > yl[2], na.rm = TRUE)
  p <- ggplot(d, aes(trait_lab, delta, fill = trait)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.4, colour = "grey20") +
    facet_grid(height ~ Cluster) +
    coord_cartesian(ylim = yl) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_fill_manual(values = trait_fill) +
    labs(x = NULL, y = ylab_expr, title = "All period — sampling-height comparison",
         caption = sprintf("Top: legacy nair==1 (height scales with Hmax). Bottom: fixed 1 m (HOBO). Δv<0 = buffers. P1 sparse→P4 dense.%s",
                           if (nc>0) sprintf(" %d points clipped.", nc) else "")) +
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold", size = 15),
          axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
          plot.caption = element_text(size = 10.5, colour = "grey40"),
          legend.position = "none")
  ggsave(file.path(OUT, paste0(fname, ".png")), p, width = 16, height = 8.5, dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(fname, ".pdf")), p, width = 16, height = 8.5, device = cairo_pdf)
  cli_alert_success("Saved {fname}")
}

cli_h1("Height-comparison figures")
make_cmp("Delta*T[max]~(degree*C)",
         expression(Delta[v]^"LOO" ~ "on " * Delta*T[max] ~ "(" * degree * "C)"),
         "fig_loo_height_compare_tmax")
make_cmp("Delta*VPD[max]~(kPa)",
         expression(Delta[v]^"LOO" ~ "on " * Delta*VPD[max] ~ "(kPa)"),
         "fig_loo_height_compare_vpd")

cli_h2("LAI median Δv: does the sign-flip survive the height change? (all period)")
print(long[trait == "LAI" & period == "All period",
           .(median = round(median(delta, na.rm=TRUE),3)), by=.(metric, height, Cluster)][order(metric, Cluster, height)])
