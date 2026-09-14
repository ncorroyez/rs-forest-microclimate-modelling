# ==============================================================================
# Vertical profiles of air temperature, wind, RH and VPD inside/above the canopy,
# REF coalition, per archetype (P1->P4), Real vs Uniform LAD. JJAS daytime.
#
# Produces TWO y-axis variants:
#   * NORMALIZED : z / h_canopy  (1.0 = canopy top)                -> _norm
#   * DENORMALIZED (metres) : z = relative_height * h_canopy       -> _metres
#     where h_canopy is the per-archetype median veget_height_top
#     (== archetype Hmax: C1=34, C2=15, C3=21, C4=29 m; identical for real/uniform).
#
# ------------------------------------------------------------------------------
# IMPORTANT PROVENANCE NOTE (read before quoting this figure):
#   The rest of Chapter 1 runs on the NATIVE20 pipeline (traits computed natively
#   at 20 m squares). NATIVE per-archetype REAL-vs-UNIFORM profile sims DO NOT
#   EXIST yet — the native20 rebuild (2026-07-17) produced only per-HOBO-pixel
#   forward-inclusion coalitions (out_files/musica_native20_forward/) and tiles.
#   So this figure still reads the LEGACY archetype sims (H1_archetypes_v9,
#   legacy binary, May 2026). Re-pointing to native would require a native
#   archetype ARCH_ re-run (see report). The relabel below (relabel_cluster())
#   is numerically IDENTICAL to the old hand map c("1"=P4,"2"=P2,"3"=P1,"4"=P3):
#   the labels were NOT scrambled; only the sim SOURCE is legacy.
#
# VPD/RH are DERIVED from Tair_z (K), wair_z (mol/mol), wind_z (m/s), which
# MuSICA does NOT output directly.
#
# Output : outputs/figures_chap1/fig_vertical_profiles_{norm,metres}.{png,pdf}
# ==============================================================================

suppressMessages({
  library(ncdf4); library(here); library(cli)
  library(tidyverse); library(lubridate)
})
source(here::here("R/cluster_relabel.R"))       # relabel_cluster() (single source of truth for labels)
source(here::here("scripts/_article_style.R"))  # PAL_CLUSTER, theme_article — palette authority (sourced last)

OUT <- here::here("outputs/figures_chap1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- config ------------------------------------------------------------------
P_HPA    <- 1013            # surface pressure approximation (hPa)
MONTHS   <- 6:9             # June - September
HOURS    <- 10:16           # daytime midday window (local)
CLUSTERS <- 1:4             # Arch_C1..4 (raw k-means IDs) -> relabel_cluster -> P1..P4
SIM_DIR  <- "out_files/H1_archetypes_v9"   # LEGACY (native archetype sims absent, see note)

# LAD-profile shape contrast at fixed LAI / Hmax / fCover: real (r_r_r_r) vs
# vertically uniform (r_r_r_u). Solid-minus-dashed gap = pure vertical-profile (H1)
# effect, isolated from the density difference between archetypes.
SHAPES <- c("Real" = "r_r_r_r", "Uniform LAD" = "r_r_r_u")

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
  hcan <- median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)  # canopy top (m)

  keep <- month(tvec) %in% MONTHS & hour(tvec) %in% HOURS
  Tc   <- Tk[, keep] - 273.15                           # nair x ntime
  e    <- (wmr[, keep] / (1 + wmr[, keep])) * P_HPA     # actual vapour press, hPa
  es   <- esat_hpa(Tc)
  vpd  <- es - e;  vpd[vpd < 0]   <- 0                  # hPa (keep matrix dims)
  relh <- 100 * e / es;  relh[relh > 100] <- 100        # RH % (keep matrix dims)
  wd   <- wnd[, keep]

  tibble(
    rel_height = as.numeric(rh),
    h_canopy   = hcan,                                     # per-archetype scalar (recycled)
    Temp       = as.numeric(rowMeans(Tc,   na.rm = TRUE)),
    Wind       = as.numeric(rowMeans(wd,   na.rm = TRUE)),
    RH         = as.numeric(rowMeans(relh, na.rm = TRUE)),
    VPD        = as.numeric(rowMeans(vpd,  na.rm = TRUE)) / 10  # hPa -> kPa
  )
}

# --- gather ------------------------------------------------------------------
cli_h1("Vertical profiles : Temp / Wind / RH / VPD (legacy archetype sims)")
cli_alert_info("SIM source: {SIM_DIR}  (native archetype real-vs-uniform sims absent)")
rows <- list()
for (cl in CLUSTERS) for (sh in names(SHAPES)) {
  f <- here::here(SIM_DIR, sprintf("Arch_C%d", cl),
                  sprintf("musica_out_Arch_C%d_ARCH_%s.nc", cl, SHAPES[[sh]]))
  if (!file.exists(f)) { cli_alert_warning("missing {f}"); next }
  d <- read_profiles(f)
  d$Profile <- as.character(relabel_cluster(cl))
  d$Shape   <- sh
  rows[[length(rows) + 1L]] <- d
  cli_alert_success("Arch_C{cl}/{sh} -> {d$Profile[1]} : {nrow(d)} levels, h_canopy={round(d$h_canopy[1],1)} m")
}
df <- bind_rows(rows) %>%
  mutate(z_m = rel_height * h_canopy)          # denormalized height (metres)

# per-archetype canopy top (m) for the metres reference lines (real == uniform)
cap <- df %>% distinct(Profile, h_canopy) %>%
  mutate(Profile = factor(Profile, levels = c("P1","P2","P3","P4")))

long <- df %>%
  pivot_longer(c(Temp, Wind, RH, VPD), names_to = "var", values_to = "value") %>%
  mutate(
    var     = factor(var, levels = c("Temp", "Wind", "RH", "VPD"),
                     labels = c("Air~temperature~(degree*C)",
                                "Wind~(m~s^{-1})",
                                "RH~('%')",
                                "VPD~(kPa)")),
    Profile = factor(Profile, levels = c("P1", "P2", "P3", "P4")),
    Shape   = factor(Shape, levels = c("Real", "Uniform LAD"))
  )

pal <- PAL_CLUSTER   # P1 sparse (warm) -> P4 dense (cool/green), palette partagée

# ---- builder shared by both variants ----------------------------------------
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

# ---- NORMALIZED variant (z / h_canopy) --------------------------------------
p_norm <- ggplot(long, aes(x = value, y = rel_height,
                           colour = Profile, linetype = Shape,
                           group = interaction(Profile, Shape))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  annotate("text", x = Inf, y = 1.02, label = "canopy top",
           hjust = 1.05, vjust = 0, size = 4, colour = "grey40")
p_norm <- base_layers(p_norm) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.25)) +
  labs(x = NULL, y = "Relative height  z / h[canopy]")

# ---- DENORMALIZED variant (metres) ------------------------------------------
# canopy top is per-archetype (its own Hmax), so use coloured per-archetype
# segments — NOT a single line at a fixed height (that would be wrong here).
p_metres <- ggplot(long, aes(x = value, y = z_m,
                             colour = Profile, linetype = Shape,
                             group = interaction(Profile, Shape))) +
  geom_hline(data = cap, aes(yintercept = h_canopy, colour = Profile),
             linetype = "dotted", linewidth = 0.5, alpha = 0.6,
             inherit.aes = FALSE)
p_metres <- base_layers(p_metres) +
  labs(x = NULL, y = "Height above ground  z (m)")

# --- save --------------------------------------------------------------------
save_both <- function(p, stem) {
  ggsave(file.path(OUT, paste0(stem, ".png")), p,
         width = 16, height = 5.6, dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p,
         width = 16, height = 5.6, device = cairo_pdf)
  cli_alert_success("Saved {stem} (png + pdf)")
}
save_both(p_norm,   "fig_vertical_profiles_norm")
save_both(p_metres, "fig_vertical_profiles_metres")

# --- denormalization sanity table --------------------------------------------
cli_h2("Denormalization landmarks (per archetype)")
sanity <- df %>% group_by(Profile) %>%
  summarise(h_canopy_m = first(h_canopy),
            z_min_m    = min(z_m),
            z_top_m    = 1 * first(h_canopy),
            z_max_m    = max(z_m), .groups = "drop") %>%
  arrange(Profile)
print(as.data.frame(sanity %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
