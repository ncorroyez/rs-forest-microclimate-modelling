# ==============================================================================
# MEB 2026 v8 — overnight batch
#
#   1. Re-aggregate 16-coalition cache for HOBO (add missing 0000/1000/1100/1110)
#   2. Map of Blois forest with 400 cLHS + 53 HOBO locations
#   3. Forward-selection HOBO figure (LOO ranking : fCover→LAI→Hmax→LAD)
#   4. Gradient heatmap LOO × 4 archetypes
#   5. Linear regression overlay added to HOBO Spearman figure
#
# No file deletion ; all outputs land in outputs/figs_MEB2026_final/ as v8 names.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot); library(ncdf4); library(sf)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
OUT  <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

VARS         <- c("LAI","Hmax","fCover","LAD")
.CLUSTER_SHORT <- c("1"="C1","2"="C2","3"="C3","4"="C4")
.PAL_CLUSTER   <- c("C1"="#E69F00","C2"="#0072B2","C3"="#009E73","C4"="#CC79A7")

# ============================================================================
# (1) Re-aggregate missing 16-coalition rows for HOBO
# ============================================================================
reaggregate_hobo <- function() {
  cli_h1("Re-aggregate HOBO 16-coalition cache (proper Delta_Tmax via df_macro)")
  suppressMessages({
    library(musica.tools); library(dplyr); library(lubridate); library(tidyr)
  })
  source(here::here("R/config.R"), local = FALSE)
  source(here::here("R/io.R"),     local = FALSE)
  source(here::here("R/musica.R"), local = FALSE)
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

  problem_ids <- c("41_04","41_09","41_10","41_11","41_17","41_26","41_36",
                    "41_39","41_40","41_42","41_54","41_59")
  h1f_map <- list(
    "H1f_0_Null_baseline"   = "0000",
    "H1f_1_LAI_only"        = "1000",
    "H1f_2_LAI_Hmax"        = "1100",
    "H1f_3_LAI_Hmax_fCover" = "1110",
    "H1f_4_Full_real"       = "1111")

  new_rows <- list()
  for (sid in problem_ids) {
    for (h1f_dir in names(h1f_map)) {
      bc <- h1f_map[[h1f_dir]]
      f <- file.path("out_files/musica_hobo_validation", h1f_dir,
                      paste0("musica_out_HOBO_", sid, ".nc"))
      if (!file.exists(f)) next
      daily <- tryCatch(
        extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                z_target = CFG$tair_target_height),
        error = function(e) NULL)
      if (is.null(daily) || nrow(daily) == 0) next
      new_rows[[length(new_rows) + 1L]] <- data.table(
        id_plot = sid, bit_code = bc,
        Tmax_mean = mean(daily$Delta_Tmax, na.rm = TRUE),
        Tmax_P90  = quantile(daily$Delta_Tmax, 0.90, na.rm = TRUE, names = FALSE),
        n_days = nrow(daily))
    }
  }
  new_dt <- rbindlist(new_rows)
  cli_alert("Extracted {nrow(new_dt)} new (sensor × coalition) rows from H1f NCs")
  print(new_dt[id_plot == "41_04"])

  agg <- readRDS(file.path(DATA, "DT_agg_HOBO_floor05_16coalitions.rds"))
  agg[, key := paste(id_plot, bit_code, sep = "_")]
  new_dt[, key := paste(id_plot, bit_code, sep = "_")]
  # Replace any spurious K-valued rows then add new
  agg <- agg[!key %in% new_dt$key]
  agg_full <- rbind(agg, new_dt)
  agg_full[, key := NULL]
  setkey(agg_full, id_plot, bit_code)
  saveRDS(agg_full, file.path(DATA, "DT_agg_HOBO_floor05_16coalitions_v2.rds"))
  cnt <- agg_full[, .N, by = id_plot]
  cli_alert_success("Full coverage : {sum(cnt$N == 16)} / 53 sensors with 16 coalitions")
  cli_alert("NULL Tmax_mean range now : [{round(min(agg_full[bit_code=='0000']$Tmax_mean),2)}, {round(max(agg_full[bit_code=='0000']$Tmax_mean),2)}]")
  invisible(agg_full)
}

# ============================================================================
# (2) Map of cLHS + HOBO locations
# ============================================================================
build_map <- function() {
  cli_h1("Map of Blois forest plots")
  clhs <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  hobo_meta <- sf::st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
  hobo_meta <- sf::st_drop_geometry(hobo_meta)
  hobo_cache <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  hobo_sub <- hobo_meta[hobo_meta$id_plot %in% hobo_cache$id_plot, ]
  cli_alert("cLHS plots : {nrow(clhs)} ; HOBO sensors : {nrow(hobo_sub)}")

  pal_pts <- c("cLHS plots (n = 400)" = "#56B4E9",
                 "HOBO sensors (n = 53)" = "#D55E00")
  p <- ggplot() +
    geom_point(data = clhs,
                 aes(x = x, y = y, colour = "cLHS plots (n = 400)"),
                 size = 1.5, alpha = 0.65, shape = 16) +
    geom_point(data = hobo_sub,
                 aes(x = coord_x_utm31n, y = coord_y_utm31n,
                       colour = "HOBO sensors (n = 53)"),
                 size = 3.2, shape = 17) +
    scale_colour_manual(values = pal_pts, name = NULL) +
    coord_equal() +
    labs(x = "Easting (UTM 31N, m)", y = "Northing (UTM 31N, m)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
           legend.text = element_text(size = 13),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_map_blois_locations.png")
  ggsave(fig_path, p, width = 8, height = 8, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (3) Forward-selection HOBO figure
# ============================================================================
# LOO ranking on cLHS (averaged MAE) : fCover=1, LAI=2, Hmax=3, LAD=4
# Forward path :
#   NULL (0000) -> +fCover (0010) -> +fCover+LAI (1010)
#                -> +fCover+LAI+Hmax (1110) -> REF (1111)
# At each step, compare simulated Tmax_mean (per sensor) to observed Tmacro - dTmax
# We display sim_Tmax vs observed_Tmax (= Tmacro - HOBO_dTmax).
# Reference Tmacro comes from the ERA5 forcing or a stored mean — we use
# Tmax_NULL + dTmax_NULL as proxy.  Simpler : use observed buffering ΔT (= dTmax)
# vs simulated buffering ΔT_sim = T_macro - Tmax_sim.  Since T_macro is constant
# across coalitions for the same plot, the comparison is equivalent to comparing
# Tmax_sim to (Tmax_macro - dTmax_obs).  We need T_macro per sensor.
# ----------------------------------------------------------------------------
build_forward_selection <- function() {
  cli_h1("Forward-selection HOBO figure")
  agg <- readRDS(file.path(DATA, "DT_agg_HOBO_floor05_16coalitions_v2.rds"))
  cg <- readRDS("outputs/lovb/data/DT_cross_gam_targets.rds")
  # Convention : ΔTmax = Tmicro - Tmacro  (negative = canopy buffers)
  hobo_obs <- as.data.table(cg$HOBO)[, .(id_plot = id, dT_obs = -dTmax_mean)]

  # Forward path : signed Avg Δ_v on the 4 profiles, from most anti-buffer
  # (positive Avg) to strongest buffer (negative Avg).
  # Avg Δ_v : fCover = +0.07 > LAD = +0.03 > Hmax = -0.03 > LAI = -0.12
  # Bit positions are (LAI, Hmax, fCover, LAD) = (1,2,3,4)
  # NULL=0000, +fCover=0010, +fCover+LAD=0011, +fCover+LAD+Hmax=0111, REF=1111
  steps <- data.table(
    step    = 1:5,
    bit     = c("0000","0010","0011","0111","1111"),
    label   = c("Baseline",
                  "+ fCover",
                  "+ fCover + LAD",
                  "+ fCover + LAD + Hmax",
                  "+ fCover + LAD + Hmax + LAI"),
    short   = c("NULL","+fCover","+LAD","+H[max]","+LAI"))
  steps[, label := factor(label, levels = label)]

  # Cache convention : Tmax_mean = mean(Tmicro_sim - Tmacro)
  # We use the same Tmicro - Tmacro convention (negative = buffer) on both axes.
  rows <- list()
  for (i in seq_len(nrow(steps))) {
    bc <- steps$bit[i]
    sub <- agg[bit_code == bc, .(id_plot, dT_sim = Tmax_mean)]
    sub <- merge(sub, hobo_obs, by = "id_plot")
    sub[, step := i]
    sub[, label := steps$label[i]]
    rows[[i]] <- sub
  }
  long <- rbindlist(rows)
  long[, label := factor(label, levels = levels(steps$label))]

  # Per-step fit stats : Pearson r, R² (from linear reg), RMSE and MAE.
  # New convention : x = observed, y = simulated
  fits <- long[, {
    ok <- !is.na(dT_sim) & !is.na(dT_obs)
    x <- dT_obs[ok]; y <- dT_sim[ok]
    r   <- suppressWarnings(cor(x, y, method = "pearson"))
    fit <- lm(y ~ x)
    co <- coef(fit); R2 <- summary(fit)$r.squared
    rmse <- sqrt(mean((y - x)^2))
    mae  <- mean(abs(y - x))
    .(n = length(x), r = r,
        slope = co[2], intercept = co[1],
        R2 = R2, RMSE = rmse, MAE = mae)
  }, by = .(step, label)]
  cli_alert("Forward-selection fit per step :"); print(fits)
  fwrite(fits, file.path(OUT, "tab_forward_selection_fits.csv"))

  # Common axis range
  xy_lim <- range(c(long$dT_sim, long$dT_obs), na.rm = TRUE)
  xy_pad <- diff(xy_lim) * 0.06
  xy_lim <- c(xy_lim[1] - xy_pad, xy_lim[2] + xy_pad)

  ann <- fits[, .(label,
                    lab_r     = sprintf("italic(r) == %+.2f", r),
                    lab_rmse  = sprintf("RMSE == %.2f * ' °C'", RMSE),
                    lab_mae   = sprintf("MAE == %.2f * ' °C'", MAE))]

  p <- ggplot(long, aes(x = dT_obs, y = dT_sim)) +
    geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
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
    labs(x = bquote("Field " ~ Delta * T[max] ~ "(°C)"),
          y = bquote("Simulated " ~ Delta * T[max] ~ "(°C)")) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 13),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_forward_selection_HOBO.png")
  ggsave(fig_path, p, width = 18, height = 5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")

  # Also save a compact "R² evolution" line plot
  # R² and RMSE evolution
  fits_evo <- melt(fits[, .(step, label, R2, RMSE)], id.vars = c("step","label"),
                      variable.name = "metric", value.name = "value")
  fits_evo[, metric_lab := ifelse(metric == "R2", "R²", "RMSE (°C)")]
  fits_evo[, metric_lab := factor(metric_lab, levels = c("R²","RMSE (°C)"))]
  p_evo <- ggplot(fits_evo, aes(x = step, y = value)) +
    geom_line(linewidth = 1.1, colour = "#2C5F2D") +
    geom_point(size = 4, colour = "#2C5F2D") +
    facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
    scale_x_continuous(breaks = 1:5,
                        labels = c("baseline","+fCover","+LAI","+Hmax","REF")) +
    labs(x = NULL, y = NULL,
          title = "Forward-selection — predictive skill on HOBO buffering") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 13),
           axis.text.x = element_text(face = "bold", size = 11),
           strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 12),
           panel.grid.minor = element_blank())
  fig_path2 <- file.path(OUT, "fig_forward_selection_evolution.png")
  ggsave(fig_path2, p_evo, width = 11, height = 5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path2}}")
}

# ============================================================================
# (4) Heatmap LOO × 4 archetypes (continuous gradient)
# ============================================================================
build_heatmap_archetypes <- function() {
  cli_h1("Heatmap LOO × 4 archetypes (gradient)")
  arch <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
  cli_alert("Archetype cache cols : {paste(names(arch), collapse=', ')}")
  cli_alert("Archetype cache rows : {nrow(arch)}")
  # Expect cols Tmax_mean_LOVB_<var> etc per archetype
  print(arch)

  # Build long table : archetype × var × |Δ_v_LOVB|
  long_rows <- list()
  for (v in VARS) {
    col_lovb <- paste0("Tmax_mean_LOVB_", v)
    if (!col_lovb %in% names(arch)) {
      # try alt name : Delta_<v>_mean
      delta_col <- paste0("Delta_", v, "_mean")
      if (delta_col %in% names(arch)) {
        long_rows[[length(long_rows) + 1L]] <- data.table(
          Cluster = arch$Cluster, variable = v,
          delta = arch[[delta_col]])
        next
      }
    }
    if (col_lovb %in% names(arch) && "Tmax_mean_REF" %in% names(arch)) {
      long_rows[[length(long_rows) + 1L]] <- data.table(
        Cluster = arch$Cluster, variable = v,
        delta = arch[[col_lovb]] - arch$Tmax_mean_REF)
    }
  }
  long <- rbindlist(long_rows)
  long[, Cluster_lab := paste0("C", as.character(Cluster))]
  long[, Cluster_lab := factor(Cluster_lab, levels = paste0("C", 1:4))]
  long[, variable := factor(variable, levels = VARS)]
  long[, var_lab := ifelse(variable == "Hmax", "H[max]", as.character(variable))]
  long[, var_lab := factor(var_lab, levels = c("LAI","H[max]","fCover","LAD"))]

  cli_alert("Δ_v range : [{round(min(long$delta),2)}, {round(max(long$delta),2)}]")

  # Map sign to convention C : positive = v contributes to buffering
  # Δ_v = Tmax_LOVB - Tmax_REF should be > 0 when removing v warms canopy (= v buffers)
  # We use signed value for the gradient.

  pal_gradient <- c("#2166AC","#67A9CF","#D1E5F0","#FDDBC7","#EF8A62","#B2182B")
  max_abs <- max(abs(long$delta), na.rm = TRUE) * 1.05

  p <- ggplot(long, aes(x = Cluster_lab, y = var_lab, fill = delta)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f", delta)),
               fontface = "bold", size = 5,
               colour = ifelse(abs(long$delta) > max_abs * 0.55,
                                 "white", "grey15")) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                            midpoint = 0, limits = c(-max_abs, max_abs),
                            name = expression(Delta[v]^"LOO" ~ "(°C)")) +
    scale_y_discrete(limits = rev,
                       labels = function(x) parse(text = as.character(x))) +
    labs(x = "Archetype", y = NULL) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(face = "bold", size = 14),
           axis.text.y = element_text(face = "bold", size = 14),
           panel.grid = element_blank(),
           legend.position = "right")
  fig_path <- file.path(OUT, "fig_heatmap_LOO_archetypes_gradient.png")
  ggsave(fig_path, p, width = 8, height = 6, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (5) HOBO figure with both Spearman AND linear regression annotations
# ============================================================================
build_hobo_with_linreg <- function() {
  cli_h1("HOBO figure with Spearman + linear regression overlay")
  DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  cg <- readRDS("outputs/lovb/data/DT_cross_gam_targets.rds")
  hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT_h <- merge(DT_h, hobo_t, by = "id_plot")
  for (v in VARS) {
    DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - Tmax_mean_NULL]
    DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
  }
  DT_shap <- readRDS(file.path(DATA, "DT_shapley_per_HOBO_floor05.rds"))
  shap_w <- dcast(DT_shap, id_plot ~ variable, value.var = "phi")
  setnames(shap_w, c("LAI","Hmax","fCover","LAD"),
            paste0("Shapley_", c("LAI","Hmax","fCover","LAD")))
  DT_h <- merge(DT_h, shap_w, by = "id_plot", all.x = TRUE)

  long <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
    rbindlist(lapply(VARS, function(v) data.table(
      method = m, Variable = v, dT = DT_h$dT, y = DT_h[[paste0(m, "_", v)]])))))
  long[, method_disp := ifelse(method == "LOVB", "LOO", method)]
  long[, method_disp := factor(method_disp, levels = c("LOO","LVA","Shapley"))]
  long[, Variable_lab := ifelse(Variable == "Hmax", "H[max]", Variable)]
  long[, Variable_lab := factor(Variable_lab,
                                    levels = c("LAI","H[max]","fCover","LAD"))]
  long <- long[!is.na(y) & !is.na(dT)]

  spear_boot <- function(x, y, B = 1000L, seed = 42L) {
    n <- length(x); if (n < 3) return(list(rho=NA,p=NA,ci_lo=NA,ci_hi=NA))
    rho_hat <- suppressWarnings(cor(x, y, method = "spearman"))
    p_val <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
    set.seed(seed); rb <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      rb[b] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
    }
    ci <- quantile(rb, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    list(rho = rho_hat, p = p_val, ci_lo = ci[1], ci_hi = ci[2])
  }
  lin_fit <- function(x, y) {
    fit <- lm(y ~ x); co <- coef(fit); s <- summary(fit)
    list(slope = co[2], intercept = co[1],
          R2 = s$r.squared,
          p_slope = s$coefficients[2, 4])
  }
  ann <- long[, {sp <- spear_boot(dT, y); lr <- lin_fit(dT, y)
                  .(rho = sp$rho, p = sp$p, ci_lo = sp$ci_lo, ci_hi = sp$ci_hi,
                     slope = lr$slope, R2 = lr$R2, p_slope = lr$p_slope)},
              by = .(method_disp, Variable_lab)]
  ann[, rho_lab   := sprintf("rho == %+.2f", rho)]
  ann[, R2_lab    := sprintf("R^2 == %.2f", R2)]
  ann[, slope_lab := sprintf("slope == %+.2f", slope)]

  p <- ggplot(long, aes(x = dT, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.8, alpha = 0.20) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = rho_lab),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4.6,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = R2_lab),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 4.2,
               inherit.aes = FALSE, colour = "grey25") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = slope_lab),
               parse = TRUE, hjust = -0.08, vjust = 4.6, size = 4.0,
               inherit.aes = FALSE, colour = "grey35") +
    facet_grid(method_disp ~ Variable_lab, scales = "free", switch = "y",
                 labeller = labeller(Variable_lab = label_parsed)) +
    labs(x = bquote(Delta * T[max] ~ "(°C)"),
          y = bquote("Simulated  " * Delta[v] ~ "(°C)")) +
    theme_bw(base_size = 16) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 15),
           strip.placement = "outside",
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_HOBO_final_linreg.png")
  ggsave(fig_path, p, width = 16, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  fwrite(ann, file.path(OUT, "tab_HOBO_spearman_linreg.csv"))
}

# ============================================================================
if (sys.nframe() == 0) {
  reaggregate_hobo()
  build_map()
  build_forward_selection()
  build_heatmap_archetypes()
  build_hobo_with_linreg()
  cli_alert_success("v8 overnight batch complete.")
}
