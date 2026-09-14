# ==============================================================================
# Vertical profiles of air temperature, wind, RH and VPD inside/above the canopy,
# REF coalition, per cluster (P1->P4), LEGACY binary (article reference).
#
# Reads existing archetype REF outputs (no re-run). VPD/RH are DERIVED from the
# profile variables Tair_z (K), wair_z (mol/mol mixing ratio), wind_z (m/s),
# which MuSICA does NOT output directly.
#
# Summer 2021, daytime window (10-16 h local). relative_height = z / h_canopy,
# so 1.0 = canopy top; nair==1 (~1 m) = sub-canopy.
#
# Output : outputs/figs_MEB2026_final/fig_vertical_profiles_legacy.{png,pdf}
# ==============================================================================

suppressMessages({
  library(ncdf4); library(here); library(cli)
  library(tidyverse); library(lubridate)
})
source(here::here("scripts/_article_style.R"))   # PAL_CLUSTER, theme_article (style partagé, anti-dérive)

OUT <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- config ------------------------------------------------------------------
P_HPA       <- 1013            # surface pressure approximation (hPa)
MONTHS      <- 6:9             # June - September
HOURS       <- 10:16          # daytime midday window (local)
CLUSTERS    <- 1:4            # Arch_C1..4  ->  P1..P4
VERSIONS <- list(
  V9  = "out_files/H1_archetypes_v9"        # legacy binary only (article reference)
)

# --- helpers -----------------------------------------------------------------
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))   # Tetens, hPa

read_profiles <- function(nc_path) {
  nc  <- nc_open(nc_path)
  on.exit(nc_close(nc))
  rh   <- ncvar_get(nc, "relative_height")              # nair
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
    rel_height = as.numeric(rh),
    Temp       = as.numeric(rowMeans(Tc,   na.rm = TRUE)),  # air temperature (°C), per level
    Wind       = as.numeric(rowMeans(wd,   na.rm = TRUE)),  # mean over time, per level
    RH         = as.numeric(rowMeans(relh, na.rm = TRUE)),
    VPD        = as.numeric(rowMeans(vpd,  na.rm = TRUE)) / 10  # hPa -> kPa
  )
}

# --- gather ------------------------------------------------------------------
cli_h1("Vertical profiles : Temp / Wind / RH / VPD (legacy)")
# Arch_C%d dirs use the RAW k-means codes; map to display labels via cluster_relabel
# (raw 3→P1 sparsest, 2→P2, 4→P3, 1→P4 densest). Direct C%d→P%d was WRONG (scrambled).
RELAB <- c("1" = "P4", "2" = "P2", "3" = "P1", "4" = "P3")
# LAD-profile shape contrast at fixed LAI / Hmax / fCover: real (r_r_r_r) vs vertically
# uniform (r_r_r_u). The solid-minus-dashed gap is the pure vertical-profile (H1) effect,
# isolated from the density difference between types.
SHAPES <- c("Real" = "r_r_r_r", "Uniform LAD" = "r_r_r_u")
rows <- list()
for (vn in names(VERSIONS)) for (cl in CLUSTERS) for (sh in names(SHAPES)) {
  f <- here::here(VERSIONS[[vn]],
                  sprintf("Arch_C%d", cl),
                  sprintf("musica_out_Arch_C%d_ARCH_%s.nc", cl, SHAPES[[sh]]))
  if (!file.exists(f)) { cli_alert_warning("missing {f}"); next }
  d <- read_profiles(f)
  d$Version <- vn
  d$Profile <- unname(RELAB[as.character(cl)])
  d$Shape   <- sh
  rows[[length(rows) + 1L]] <- d
  cli_alert_success("{vn} Arch_C{cl}/{sh} -> {RELAB[as.character(cl)]} : {nrow(d)} levels")
}
df <- bind_rows(rows)

long <- df %>%
  pivot_longer(c(Temp, Wind, RH, VPD), names_to = "var", values_to = "value") %>%
  mutate(
    var     = factor(var, levels = c("Temp", "Wind", "RH", "VPD"),
                     labels = c("Air~temperature~(degree*C)",
                                "Wind~(m~s^{-1})",
                                "RH~('%')",
                                "VPD~(kPa)")),
    Version = factor(Version, levels = "V9",
                     labels = "V9 (legacy)"),
    Profile = factor(Profile, levels = c("P1", "P2", "P3", "P4")),
    Shape   = factor(Shape, levels = c("Real", "Uniform LAD"))
  )

# --- colours : P1 sparse (warm) -> P4 dense (cool/green) — palette partagée ----
pal <- PAL_CLUSTER

p <- ggplot(long, aes(x = value, y = rel_height,
                      colour = Profile, linetype = Shape,
                      group = interaction(Profile, Shape))) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  annotate("text", x = Inf, y = 1.02, label = "canopy top",
           hjust = 1.05, vjust = 0, size = 4, colour = "grey40") +
  geom_path(linewidth = 1.1) +
  geom_point(aes(shape = Shape), size = 1.8) +
  facet_wrap(~ var, nrow = 1, scales = "free_x",
             labeller = labeller(var = label_parsed)) +
  scale_colour_manual(values = pal, name = "Cluster") +
  scale_linetype_manual(values = c("Real" = "solid", "Uniform LAD" = "22"),
                        name = "LAD profile") +
  scale_shape_manual(values = c("Real" = 16, "Uniform LAD" = 1),
                     name = "LAD profile") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.25)) +
  labs(x = NULL, y = "Relative height  z / h[canopy]") +
  theme_article(16) +                                     # style maison partagé
  theme(legend.position = "bottom", legend.direction = "horizontal")

# save to the legacy figs dir AND the pipeline annex
ANNEX <- here::here("outputs/figures_pipeline/annex")
dir.create(ANNEX, recursive = TRUE, showWarnings = FALSE)
for (d in c(OUT, ANNEX)) {
  ggsave(file.path(d, "fig_vertical_profiles_legacy.png"),
         p, width = 16, height = 5.6, dpi = 300, bg = "white")
  ggsave(file.path(d, "fig_vertical_profiles_legacy.pdf"),
         p, width = 16, height = 5.6, device = cairo_pdf)
}
cli_alert_success("Saved fig_vertical_profiles_legacy (png + pdf) to figs + pipeline annex")

# --- sub-canopy (nair==1) summary table --------------------------------------
sub <- df %>% group_by(Version, Profile, Shape) %>% slice_min(rel_height, n = 1) %>%
  ungroup() %>% arrange(Version, Profile, Shape)
cli_h2("Sub-canopy level (nair==1, ~1 m) midday means")
print(as.data.frame(sub %>% mutate(across(c(Temp, Wind, RH, VPD), ~round(.x, 2)))))
write_csv(sub, file.path(OUT, "tab_subcanopy_legacy.csv"))
