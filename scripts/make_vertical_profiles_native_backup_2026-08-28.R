# ==============================================================================
# Vertical profiles of air temperature, wind, RH and VPD inside/above the canopy,
# NATIVE20 simulations, per cluster (P1->P4), Real LAD vs Uniform LAD.
#
# Reconstructed from the existing forward-coalition sims (NO new MuSICA run):
#   Real LAD    = coalition 1111 (all four traits real)
#   Uniform LAD = coalition 1110 (LAD bit off, LAI/Hmax/fCover real)
# Bit order over traits = [LAI, Hmax, fCover, LAD]; bit=1 -> trait REAL.
# The 1111-vs-1110 contrast isolates the pure LAD vertical-profile-shape effect
# at fixed LAI/Hmax/fCover, per pixel (mirrors legacy r_r_r_r vs r_r_r_u).
#
# VPD/RH are DERIVED from the profile variables Tair_z (K), wair_z (mol/mol
# mixing ratio), wind_z (m/s), exactly as scripts/make_vertical_profiles_wind_rh_vpd.R.
#
# Summer JJAS (months 6-9), daytime window (hours 10-16). NB: time is parsed as
# UTC and the 10-16 h window is applied on UTC (labelled "local") to stay
# byte-for-byte comparable with the legacy reference figure — France summer is
# UTC+2, so this is a deliberate simplification carried over from the reference.
#
# relative_height = z / h_canopy (in-canopy levels are fixed fractions and align
# across pixels; above-canopy levels vary with each pixel's canopy height, so we
# aggregate by LEVEL INDEX and average relative_height across pixels per level).
#
# Cluster mean = unweighted mean over the cluster's pixels (equal weight per
# pixel, regardless of each pixel's valid-timestep count).
#
# Output : outputs/figures_chap1/fig_vertical_profiles_native_norm.{png,pdf}
#          outputs/figures_chap1/fig_vertical_profiles_native_metres.{png,pdf}
# Reads :
#         out_files/musica_native20_forward/{1111,1110}/musica_out_HOBO_<id>.nc
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv      (stage A2)
# Writes: outputs/figures_chap1/fig_vertical_profiles_native_{norm,metres}.{png,pdf}
#         outputs/figures_chap1/tab_subcanopy_native.csv
#         out_files/Chapter1/tables/tab_diurnal_cycle_1m.csv  (the 10-16 h window justification)
#   Rscript scripts/make_vertical_profiles_native.R
# ==============================================================================

suppressMessages({
  library(ncdf4); library(here); library(cli)
  library(tidyverse); library(lubridate)
})
source(here::here("scripts/_article_style.R"))   # PAL_CLUSTER, theme_article (style partagé, anti-dérive)

OUT <- here::here("outputs/figures_chap1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- config ------------------------------------------------------------------
P_HPA    <- 1013            # surface pressure approximation (hPa)
MONTHS   <- 6:9             # June - September (JJAS)
HOURS    <- 10:16           # daytime midday window (UTC, labelled local; see header)
CLUSTERS <- c("P1", "P2", "P3", "P4")
SHAPES   <- c("Real" = "1111", "Uniform LAD" = "1110")   # LAD real vs uniform
FWD_DIR  <- here::here("out_files/musica_native20_forward")
CLU_CSV  <- here::here("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")

# --- helpers -----------------------------------------------------------------
#' Tetens saturation vapour pressure (hPa). Local copy: pipeline/00_config.R defines the
#' same name, and this script does not source it.
#' @param Tc air temperature, degrees C
#' @return saturation vapour pressure, hPa
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))   # Tetens, hPa

# Per-pixel, per-level (15) JJAS-daytime mean profile. Returns one row per level.
#' JJAS midday mean vertical profile of T, wind, RH and VPD for one run
#' @param nc_path path to a MuSICA output NetCDF
#' @return tibble, one row per air level: level, rel_height, vhtop, Temp, Wind, RH, VPD
read_profiles <- function(nc_path) {
  nc  <- nc_open(nc_path)
  on.exit(nc_close(nc))
  rh   <- ncvar_get(nc, "relative_height")              # nair (z/h_canopy)
  vh   <- median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)  # canopy top (m)
  tu   <- ncatt_get(nc, "time", "units")$value
  t0   <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec <- t0 + lubridate::dhours(ncvar_get(nc, "time"))
  # ncdf4 returns dims fastest-first -> Tair_z(time,nair) becomes [nair, time]
  Tk   <- ncvar_get(nc, "Tair_z")                       # nair x time (K)
  wmr  <- ncvar_get(nc, "wair_z")                       # nair x time (mol/mol)
  wnd  <- ncvar_get(nc, "wind_z")                       # nair x time (m/s)

  keep <- month(tvec) %in% MONTHS & hour(tvec) %in% HOURS
  Tc   <- Tk[, keep] - 273.15                           # nair x ntime
  e    <- (wmr[, keep] / (1 + wmr[, keep])) * P_HPA     # actual vapour press, hPa
  es   <- esat_hpa(Tc)
  vpd  <- es - e;  vpd[vpd < 0]   <- 0                  # hPa (keep matrix dims)
  relh <- 100 * e / es;  relh[relh > 100] <- 100        # RH % (keep matrix dims)
  wd   <- wnd[, keep]

  tibble(
    level      = seq_along(rh),
    rel_height = as.numeric(rh),
    vhtop      = vh,
    Temp       = as.numeric(rowMeans(Tc,   na.rm = TRUE)),  # air temperature (°C), per level
    Wind       = as.numeric(rowMeans(wd,   na.rm = TRUE)),  # per level
    RH         = as.numeric(rowMeans(relh, na.rm = TRUE)),  # %
    VPD        = as.numeric(rowMeans(vpd,  na.rm = TRUE)) / 10  # hPa -> kPa
  )
}

# --- cluster assignment ------------------------------------------------------
cli_h1("Vertical profiles : Temp / Wind / RH / VPD (NATIVE20, per cluster)")
clu <- readr::read_csv(CLU_CSV, show_col_types = FALSE) %>%
  transmute(id = id_plot, P = factor(P, levels = CLUSTERS))
n_by_cluster <- clu %>% count(P, name = "n_pixels")
cli_h2("Per-cluster pixel counts (from cluster table)")
print(as.data.frame(n_by_cluster))

# --- gather : per pixel, then average over the cluster's pixels ---------------
pixel_rows <- list()
missing    <- character(0)
for (sh in names(SHAPES)) for (cl in CLUSTERS) {
  ids <- clu$id[clu$P == cl]
  for (id in ids) {
    f <- file.path(FWD_DIR, SHAPES[[sh]], sprintf("musica_out_HOBO_%s.nc", id))
    if (!file.exists(f)) { missing <- c(missing, f); next }
    d <- read_profiles(f)
    d$Profile <- cl; d$Shape <- sh; d$id <- id
    pixel_rows[[length(pixel_rows) + 1L]] <- d
  }
}
if (length(missing)) {
  cli_alert_warning("Missing {length(missing)} file(s):")
  for (m in missing) cli_alert_warning("  {m}")
}
pix <- bind_rows(pixel_rows)

# actual pixel counts contributing per cluster/shape (after file.exists)
used_counts <- pix %>% distinct(Shape, Profile, id) %>% count(Shape, Profile, name = "n_used")
cli_h2("Pixels used per cluster x shape (after file.exists)")
print(as.data.frame(used_counts))

# cluster mean = unweighted mean over pixels, per level index
df <- pix %>%
  group_by(Profile, Shape, level) %>%
  summarise(
    rel_height = mean(rel_height, na.rm = TRUE),
    vhtop      = median(vhtop, na.rm = TRUE),   # per-cluster canopy top (m)
    Temp       = mean(Temp, na.rm = TRUE),
    Wind       = mean(Wind, na.rm = TRUE),
    RH         = mean(RH,   na.rm = TRUE),
    VPD        = mean(VPD,  na.rm = TRUE),
    .groups    = "drop"
  )

# per-cluster canopy top (median veget_height_top). rel_height=1 lands exactly
# here under the z = rel_height * vhtop transform; equals the cluster's Hmax.
canopy_top <- df %>% distinct(Profile, vhtop) %>%
  mutate(Profile = factor(Profile, levels = CLUSTERS)) %>% arrange(Profile)
cli_h2("Per-cluster canopy top (median veget_height_top, m)")
print(as.data.frame(canopy_top))

# --- reshape -----------------------------------------------------------------
#' Pivot the four profile variables long and set the factor levels for faceting
#' @param d wide profile table
#' @param yvar name of the column to use as the vertical axis
#' @return long tibble with var, Profile, Shape and y
mk_long <- function(d, yvar) {
  d %>%
    pivot_longer(c(Temp, Wind, RH, VPD), names_to = "var", values_to = "value") %>%
    mutate(
      var     = factor(var, levels = c("Temp", "Wind", "RH", "VPD"),
                       labels = c("Air~temperature~(degree*C)",
                                  "Wind~(m~s^{-1})",
                                  "RH~('%')",
                                  "VPD~(kPa)")),
      Profile = factor(Profile, levels = CLUSTERS),
      Shape   = factor(Shape, levels = c("Real", "Uniform LAD")),
      y       = .data[[yvar]]
    )
}

pal <- PAL_CLUSTER

#' Shared geometry, scales and theme for both profile panels
#' @param p a ggplot with the aesthetics already mapped
#' @return the same ggplot with the common layers added
base_layers <- function(p) {
  p +
    geom_path(linewidth = 1.1) +
    geom_point(aes(shape = Shape), size = 1.8) +
    facet_wrap(~ var, nrow = 1, scales = "free_x",
               labeller = labeller(var = label_parsed)) +
    scale_colour_manual(values = pal, name = "Cluster") +
    scale_linetype_manual(values = c("Real" = "solid", "Uniform LAD" = "22"),
                          name = "LAD profile") +
    scale_shape_manual(values = c("Real" = 16, "Uniform LAD" = 1),
                       name = "LAD profile") +
    theme_article(16) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
}

# --- (a) normalized y-axis : relative height z / h_canopy --------------------
long_norm <- mk_long(df %>% mutate(y_norm = rel_height), "y_norm")
p_norm <- base_layers(
  ggplot(long_norm, aes(x = value, y = y,
                        colour = Profile, linetype = Shape,
                        group = interaction(Profile, Shape)))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = Inf, y = 1.02, label = "canopy top",
           hjust = 1.05, vjust = 0, size = 4, colour = "grey40") +
  scale_y_continuous(breaks = seq(0, 2, 0.25)) +
  labs(x = NULL, y = "Relative height, z / canopy top")

# --- (b) absolute metres : z = rel_height * per-cluster median veget_height_top
df_m <- df %>% mutate(z_m = rel_height * vhtop)
long_m <- mk_long(df_m, "z_m")
p_metres <- base_layers(
  ggplot(long_m, aes(x = value, y = y,
                     colour = Profile, linetype = Shape,
                     group = interaction(Profile, Shape)))) +
  # each cluster's own canopy top (median veget_height_top = z at rel_height 1)
  geom_hline(data = canopy_top,
             aes(yintercept = vhtop, colour = Profile),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  labs(x = NULL, y = "Height above ground  z (m)")

# --- save --------------------------------------------------------------------
ggsave_article(file.path(OUT, "fig_vertical_profiles_native_norm"),
               p_norm, width = 16, height = 5.6)
ggsave_article(file.path(OUT, "fig_vertical_profiles_native_metres"),
               p_metres, width = 16, height = 5.6)
cli_alert_success("Saved fig_vertical_profiles_native_norm + _metres (png + pdf) to {OUT}")

# --- sub-canopy (level index 1, rel_height ~0.0285) sanity table --------------
sub <- df %>% filter(level == 1L) %>%
  arrange(Profile, Shape) %>%
  transmute(Profile, Shape,
            rel_height = round(rel_height, 4),
            z_m = round(rel_height * vhtop, 2),
            Temp = round(Temp, 2), Wind = round(Wind, 2),
            RH = round(RH, 2), VPD = round(VPD, 2))
cli_h2("Sub-canopy level (lowest, ~0.0285 z/h) JJAS midday means")
print(as.data.frame(sub))
write_csv(sub, file.path(OUT, "tab_subcanopy_native.csv"))

# ---- justification of the 10:00-16:00 window --------------------------------
# Added 2026-07-31. Fig. B2's caption says this window "brackets the 15:00 peak", a
# justification that no script computed. Derive it here from the full-model runs, so the
# window is defended by a number the reader can regenerate rather than by an assertion.
suppressPackageStartupMessages({library(lubridate); library(data.table)})
.vf <- list.files("out_files/musica_hobo_native20/1111", pattern = "\\.nc$", full.names = TRUE)
if (length(.vf)) {
  .d <- rbindlist(lapply(head(.vf, 20), function(f) {
    nc <- nc_open(f); tu <- ncatt_get(nc, "time", "units")$value
    tv <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC") + ncvar_get(nc, "time") * 3600
    rh <- ncvar_get(nc, "relative_height"); ht <- ncvar_get(nc, "veget_height_top")
    Ta <- ncvar_get(nc, "Tair_z"); nc_close(nc)
    k <- month(tv) %in% MONTHS                       # same summer as the profiles above
    i <- which.min(abs(rh * median(ht, na.rm = TRUE) - 1))   # air layer nearest 1 m
    data.table(h = hour(tv[k]), T = Ta[i, k] - 273.15)
  }))
  .c <- .d[, .(T = mean(T, na.rm = TRUE)), by = h][order(h)]
  cat(sprintf("\n1 m sub-canopy diurnal peak: %02d:00 on the forcing clock (%.2f degC); window used %02d-%02d\n",
              .c[which.max(T), h], .c[, max(T)], min(HOURS), max(HOURS)))
  fwrite(.c, "out_files/Chapter1/tables/tab_diurnal_cycle_1m.csv")
}
