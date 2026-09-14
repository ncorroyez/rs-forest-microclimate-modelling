# ==============================================================================
# SUPERSEDED 2026-08-02 — no longer the source of Figure B3.
# Replaced by scripts/make_trait_vertical_gradient_perplot.R, which builds the same
# four panels WITHOUT a common baseline: every one of the 400 cLHS plots is perturbed
# around its own measured canopy by the fixed native steps of Section 2.6. The baseline
# construction below made each increment's size depend on the distance from the design
# mean, and in P1 the step labelled dLAI actually REMOVED two thirds of the canopy, so
# it cooled while raising sub-canopy wind. On the per-plot design the same step adds
# foliage and the two signs are coherent (+0.036 degC, -0.059 m/s). Kept for provenance.
# ==============================================================================
# Marginal effect of each LiDAR trait on the WITHIN-CANOPY VERTICAL GRADIENT of
# four microclimate variables (Tair, Wind, RH, VPD), per archetype P1-P4.
#
# Committee request: "increase the focus on how the traits drive the vertical
# gradient". SHIPPED as Figure B3 of the long article (C3 of the short one); runner
# stage B10, synced by run_chapter1.R. No longer exploratory.
#
# HOW IT MUST BE READ, and why the caption says so (settled 2026-08-01). These are
# increments from ONE COMMON MEAN CANOPY, not per-plot sensitivities, so the size of a
# trait's increment depends on how far that archetype's measured value sits from the
# baseline. The profile appears to lead in P2 and P3 largely because their measured leaf
# areas (3.73 and 3.67) sit almost on the baseline (3.42), leaving the dLAI step near
# zero there, whereas the same step strips two thirds of the canopy in P1. This figure is
# therefore NOT a sensitivity ranking: Figure 4, built from the per-plot native-step
# perturbations, is the reference for that.
#
# DATA  : forward-coalition native sims
#   out_files/musica_native20_forward/<coalition>/musica_out_HOBO_<id>_*.nc
#   <coalition> = 4-bit string over [LAI, Hmax, fCover, LAD], 1=real, 0=baseline.
#   Bit convention VERIFIED against R/h1_shapley_archetypes.R (.ARCH_COAL_MAP):
#     0000 = ARCH_m_m_a_u (all baseline)   1000 = ARCH_r_m_a_u (LAI real) ...
#     LAI|Hmax|fCover|LAD, left bit = LAI.  1=real, 0=baseline (mean/uniform).
#
# INCREMENT DEFINITION (clean forward path, quantity-first order):
#   dLAI    = profile(1000) - profile(0000)
#   dHmax   = profile(1100) - profile(1000)
#   dfCover = profile(1110) - profile(1100)
#   dLAD    = profile(1111) - profile(1110)
#   Differencing is done PER PIXEL and PER VERTICAL LEVEL (nair index), then
#   averaged over the cluster's pixels. RH/VPD are DERIVED per coalition BEFORE
#   differencing (they are non-linear in Tair/wair).
#
# VERTICAL AXIS:
#   relative_height (z/h_canopy, nair index) is IDENTICAL across coalitions per
#   pixel (verified max spread 0.016, only above canopy) -> index-differencing
#   is exact.  Metres z = rel_height * veget_height_top.  NOTE: the Hmax bit
#   changes veget_height_top (baseline mean ~24 m -> real Hmax), so for the
#   dHmax increment ONLY the two coalitions live on different metre grids; we
#   plot dHmax against the per-pixel MEAN canopy height of the pair (documented
#   caveat). The normalized axis is the canonical one.
#
# nair==1 (rel_height ~ 0.0285) = sub-canopy near-ground level.
#
# Averaging window : JJAS daytime (months 6-9, hours 10-16, tz UTC as reference).
#
# OUTPUT :
#   outputs/figures_chap1/fig_trait_vertical_gradient_effects_norm.{png,pdf}
#   outputs/figures_chap1/fig_trait_vertical_gradient_effects_metres.{png,pdf}
#   + compact table (cluster x variable -> largest sub-canopy trait increment).
# ==============================================================================

suppressMessages({
  library(ncdf4); library(here); library(cli)
  library(tidyverse); library(lubridate); library(patchwork)
})
source(here::here("scripts/_article_style.R"))   # theme_article, ggsave_article

OUT <- here::here("outputs/figures_chap1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- config ------------------------------------------------------------------
P_HPA    <- 1013            # surface pressure approximation (hPa)
MONTHS   <- 6:9            # June - September
HOURS    <- 10:16          # daytime midday window (local, tz UTC as reference)
FWD_DIR  <- "out_files/musica_native20_forward"
CLUST_CSV<- "out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv"

# forward path pairs: name -> c(on-coalition, off-coalition)
INCREMENTS <- list(
  "dLAI"    = c(on = "1000", off = "0000"),
  "dHmax"   = c(on = "1100", off = "1000"),
  "dfCover" = c(on = "1110", off = "1100"),
  "dLAD"    = c(on = "1111", off = "1110")
)
COALITIONS <- unique(unlist(INCREMENTS))   # the 5 coalitions we must read

# increment palette (Okabe-Ito subset; colour = trait increment, NOT cluster)
PAL_TRAIT <- c(dLAI    = "#0072B2",   # blue
               dHmax   = "#E69F00",   # orange
               dfCover = "#009E73",   # green
               dLAD    = "#CC79A7")   # pink

# --- helpers -----------------------------------------------------------------
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))   # Tetens, hPa

# Per-level (nair index) JJAS-daytime mean profile for one nc file.
read_profiles <- function(nc_path) {
  nc  <- nc_open(nc_path); on.exit(nc_close(nc))
  rh   <- ncvar_get(nc, "relative_height")              # nair
  tu   <- ncatt_get(nc, "time", "units")$value
  t0   <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec <- t0 + lubridate::dhours(ncvar_get(nc, "time"))
  Tk   <- ncvar_get(nc, "Tair_z")                       # nair x time (K)
  wmr  <- ncvar_get(nc, "wair_z")                       # nair x time (mol/mol)
  wnd  <- ncvar_get(nc, "wind_z")                       # nair x time (m/s)
  hcan <- median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)  # canopy top (m)

  keep <- month(tvec) %in% MONTHS & hour(tvec) %in% HOURS
  Tc   <- Tk[, keep] - 273.15
  e    <- (wmr[, keep] / (1 + wmr[, keep])) * P_HPA     # actual vapour press, hPa
  es   <- esat_hpa(Tc)
  vpd  <- es - e;  vpd[vpd < 0]     <- 0                # hPa
  relh <- 100 * e / es;  relh[relh > 100] <- 100        # RH %

  tibble(
    nair   = seq_along(rh),
    rh     = as.numeric(rh),
    h_can  = hcan,
    Temp   = as.numeric(rowMeans(Tc,   na.rm = TRUE)),  # air temperature (degC)
    Wind   = as.numeric(rowMeans(wnd[, keep], na.rm = TRUE)),
    RH     = as.numeric(rowMeans(relh, na.rm = TRUE)),
    VPD    = as.numeric(rowMeans(vpd,  na.rm = TRUE)) / 10   # hPa -> kPa
  )
}

# --- cluster map -------------------------------------------------------------
cli_h1("Trait -> vertical-gradient marginal effects (forward coalitions)")
cm <- read_csv(here::here(CLUST_CSV), show_col_types = FALSE) %>%
  transmute(id = as.character(id_plot), Profile = P)
cli_alert_info("cluster map: {nrow(cm)} pixels; per-cluster: " ,
               paste(names(table(cm$Profile)), table(cm$Profile),
                     sep = "=", collapse = "  "))

# --- read every needed coalition x clustered pixel ---------------------------
rows <- list()
for (co in COALITIONS) {
  for (id in cm$id) {
    f <- here::here(FWD_DIR, co, sprintf("musica_out_HOBO_%s.nc", id))
    if (!file.exists(f)) { cli_alert_warning("missing {co}/{id}"); next }
    d <- read_profiles(f); d$coalition <- co; d$id <- id
    rows[[length(rows) + 1L]] <- d
  }
}
prof <- bind_rows(rows) %>% left_join(cm, by = "id")
# guard: equal pixel count per coalition within each cluster
chk <- prof %>% filter(nair == 1) %>% count(Profile, coalition) %>%
  pivot_wider(names_from = coalition, values_from = n)
print(chk)

# guard: index-differencing is sound even when the Hmax bit moves canopy height.
# 41_39 flips 24 m (baseline) -> ~3 m (real) between 1000 and 1100. VERIFIED:
# WITHIN-CANOPY levels (rh <= 1, nair 1..10) keep a byte-identical rh grid across
# coalitions (sub-canopy nair==1 is always 0.0285); only the ABOVE-CANOPY levels
# (rh > 1) diverge, because they sit at fixed metre offsets above a canopy top
# that the Hmax bit moves. Index-differencing is therefore EXACT within canopy
# (the task's focus); above-canopy dHmax mixes slightly different rh -> read with
# care. Assert the exact within-canopy alignment on this Hmax-flipping pixel.
.a <- prof %>% filter(id == "41_39", coalition == "1000") %>% arrange(nair)
.b <- prof %>% filter(id == "41_39", coalition == "1100") %>% arrange(nair)
if (nrow(.a) && nrow(.b)) {
  wc <- .a$rh <= 1                                   # within-canopy levels
  stopifnot(max(abs(.a$rh[wc] - .b$rh[wc])) < 1e-3,  # rh grid identical in canopy
            max(abs(.a$h_can - .b$h_can)) > 5)       # Hmax genuinely flipped
}

# --- build per-pixel, per-level increments (difference on/off) ---------------
inc_rows <- list()
for (nm in names(INCREMENTS)) {
  on_co  <- INCREMENTS[[nm]]["on"]; off_co <- INCREMENTS[[nm]]["off"]
  don <- prof %>% filter(coalition == on_co)
  dof <- prof %>% filter(coalition == off_co)
  m <- inner_join(
    don %>% select(id, Profile, nair, rh_on = rh, h_on = h_can,
                   Temp, Wind, RH, VPD),
    dof %>% select(id, nair, rh_off = rh, h_off = h_can,
                   Temp0 = Temp, Wind0 = Wind, RH0 = RH, VPD0 = VPD),
    by = c("id", "nair"))
  m <- m %>% mutate(
    increment = nm,
    dTemp = Temp - Temp0, dWind = Wind - Wind0,
    dRH   = RH   - RH0,   dVPD  = VPD  - VPD0,
    rh    = (rh_on + rh_off) / 2,          # near-identical; mean of pair
    z_m   = rh * (h_on + h_off) / 2        # metres: mean canopy height of pair
  )
  inc_rows[[length(inc_rows) + 1L]] <- m %>%
    select(id, Profile, increment, nair, rh, z_m, dTemp, dWind, dRH, dVPD)
}
inc <- bind_rows(inc_rows)

# --- average over cluster pixels at each level -------------------------------
agg <- inc %>%
  group_by(Profile, increment, nair) %>%
  summarise(rh = mean(rh), z_m = mean(z_m),
            dTemp = mean(dTemp), dWind = mean(dWind),
            dRH = mean(dRH), dVPD = mean(dVPD), n = n(),
            .groups = "drop")

long <- agg %>%
  pivot_longer(c(dTemp, dWind, dRH, dVPD),
               names_to = "var", values_to = "dval") %>%
  mutate(
    var       = factor(var, levels = c("dTemp", "dWind", "dRH", "dVPD")),
    Profile   = factor(Profile, levels = c("P1", "P2", "P3", "P4")),
    increment = factor(increment, levels = names(INCREMENTS))
  )

# x-axis titles per variable (row labels)
VAR_LAB <- c(dTemp = "Delta*Air~temperature~(degree*C)",
             dWind = "Delta*Wind~(m~s^{-1})",
             dRH   = "Delta*RH~('%')",
             dVPD  = "Delta*VPD~(kPa)")
INC_LABS <- c(dLAI = expression(Delta*LAI),  dHmax   = expression(Delta*H[max]),
              dfCover = expression(Delta*fCover), dLAD = expression(Delta*LAD))

# --- one row per microclimate variable; facet_wrap frees x PER PANEL ---------
# (facet_grid free_x frees x by column only, forcing all 4 variables of a column
#  onto one shared x-scale and crushing VPD/Tair -> use per-variable facet_wrap.)
make_row <- function(v, yvar, ylab, show_x_strip) {
  d <- long %>% filter(var == v)
  ggplot(d, aes(x = dval, y = .data[[yvar]],
                colour = increment, group = increment)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey55", linewidth = 0.5) +
    { if (yvar == "rh")
        geom_hline(yintercept = 1, linetype = "dotted",
                   colour = "grey40", linewidth = 0.5) } +
    geom_path(linewidth = 0.9) +
    geom_point(size = 1.1) +
    facet_wrap(~ Profile, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = PAL_TRAIT, name = "Trait increment",
                        labels = INC_LABS) +
    labs(x = parse(text = VAR_LAB[[v]])[[1]], y = ylab) +
    theme_article(12) +
    theme(strip.text = if (show_x_strip) element_text(face = "bold")
                       else element_blank())
}

build_fig <- function(yvar, ylab) {
  vars <- c("dTemp", "dWind", "dRH", "dVPD")
  rows <- lapply(seq_along(vars), function(i)
    make_row(vars[i], yvar, if (i == 2) ylab else "", show_x_strip = (i == 1)))
  wrap_plots(rows, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
}

p_norm   <- build_fig("rh",  "Height  z / h_canopy  (1 = canopy top)")
p_metres <- build_fig("z_m", "Height above ground (m)")

ggsave_article(file.path(OUT, "fig_trait_vertical_gradient_effects_norm"),
               p_norm, width = 11, height = 11)
ggsave_article(file.path(OUT, "fig_trait_vertical_gradient_effects_metres"),
               p_metres, width = 11, height = 11)
# also satisfy the literally-requested base path (== normalized variant)
ggsave_article(file.path(OUT, "fig_trait_vertical_gradient_effects"),
               p_norm, width = 11, height = 11)

# --- compact table: largest sub-canopy (nair==1) increment per cluster x var -
sub <- agg %>% filter(nair == 1) %>%
  pivot_longer(c(dTemp, dWind, dRH, dVPD),
               names_to = "var", values_to = "dval") %>%
  group_by(Profile, var) %>%
  slice_max(order_by = abs(dval), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Profile, var = recode(var, dTemp = "Tair(degC)", dWind = "Wind(m/s)",
                                  dRH = "RH(%)", dVPD = "VPD(kPa)"),
            largest_trait = increment, signed_effect = round(dval, 3), n_pix = n) %>%
  arrange(Profile, var)

cli_h2("Largest sub-canopy (nair==1) trait increment, per cluster x variable")
print(as.data.frame(sub), row.names = FALSE)
write_csv(sub, file.path(OUT, "tab_trait_vertical_gradient_subcanopy.csv"))

cli_alert_success("Figures + table written to {OUT}")
