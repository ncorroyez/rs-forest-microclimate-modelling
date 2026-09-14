# ==============================================================================
# Regenerate HOBO Spearman figures with new sign convention :
#   x = Tmicro_max - Tmacro_max  (negative = buffered)
#   y = Tmax_LOVB/LVA/Shapley - Tmax_REF  (positive = v contributes,
#       same as existing Delta_v/LVA_v/phi_v columns)
# Units : "°C" replaces "K" (cosmetic).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

# ---- Spearman + bootstrap helper ---------------------------------------------
spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rho_hat <- suppressWarnings(cor(x, y, method = "spearman"))
  p_val   <- suppressWarnings(cor.test(x, y, method = "spearman",
                                          exact = FALSE)$p.value)
  set.seed(seed)
  rb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    rb[b] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }
  ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
  list(rho = rho_hat, p = p_val, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

# ---- Generic figure builder ---------------------------------------------------
build_hobo_fig <- function(DT_long, vars, ann, x_label, y_label, title,
                              out_path, point_colour = "#2C5F2D",
                              loess_colour = "#FFB400", caption_extra = "") {
  p <- ggplot(DT_long, aes(x = dT_signed, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = point_colour, alpha = 0.7, size = 1.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = loess_colour, fill = loess_colour,
                  linewidth = 0.8, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.2, size = 3.1, fontface = "bold",
               inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = title,
      subtitle = "n = 53 sensors  |  Spearman rho + bootstrap CI95 (1000 reps)  |  LOESS span = 0.9",
      x = x_label, y = y_label,
      caption = paste0("Non-parametric rank correlation. LOESS shown for visual reference only.",
                         if (nchar(caption_extra) > 0L) paste0("\n", caption_extra) else "")
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 11),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
           panel.grid.minor = element_blank())
  ggsave(out_path, p, width = 11.4, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
}

# ---- Load HOBO sensor-level data + flip x sign --------------------------------
load_hobo_dt <- function() {
  DT_h <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  cg   <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
  DT <- merge(DT_h, hobo, by = "id_plot")
  # FLIP : dT_signed = Tmicro - Tmacro = -(macro - obs)
  DT[, dT_signed := -dTmax_obs]
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  # LVA per sensor (raw temp diff Tmax_LVA - Tmax_NULL)
  for (v in vars) {
    DT[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) -
                                  get("Tmax_mean_NULL")]
  }
  # Aliases :  LOVB_v = Delta_v_mean (= Tmax_LOVB - Tmax_REF in raw temps)
  for (v in vars) {
    DT[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
    DT[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) +
                                       get(paste0("LOVB_", v))) / 2]
  }
  DT[]
}

# ---- Builders for each method -------------------------------------------------
make_long <- function(DT, prefix) {
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  rbindlist(lapply(vars, function(v)
    data.table(Variable = v,
                dT_signed = DT$dT_signed,
                y         = DT[[paste0(prefix, "_", v)]])))[
                  , Variable := factor(Variable, levels = vars)][]
}

make_ann <- function(long, vars) {
  rows <- lapply(vars, function(v) {
    sub <- long[Variable == v]
    r <- spearman_boot_ci(sub$dT_signed, sub$y, B = 1000L)
    data.table(Variable = v, rho = r$rho, p = r$p,
                ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  })
  ann <- rbindlist(rows)
  ann[, label := sprintf("rho = %+.2f  p = %s\n[%+.2f ; %+.2f]",
                           rho, ifelse(p < 1e-3, "<0.001",
                                          sprintf("%.3f", p)),
                           ci_lo, ci_hi)]
  ann[, Variable := factor(Variable, levels = vars)]
  ann
}

run_all <- function() {
  cli_h1("HOBO sensor-level figures — new sign convention (x = Tmicro - Tmacro)")
  DT <- load_hobo_dt()
  vars <- c("LAI", "Hmax", "fCover", "LAD")

  x_label_signed <- "Observed  T_micro,max − T_macro,max  (°C)"

  # --- LOVB ---
  long_lovb <- make_long(DT, "LOVB")
  ann_lovb  <- make_ann(long_lovb, vars)
  cli_alert("LOVB Spearman :"); print(ann_lovb[, .(Variable, rho, p, ci_lo, ci_hi)])
  build_hobo_fig(long_lovb, vars, ann_lovb,
                  x_label = x_label_signed,
                  y_label = "Simulated  T_max,LOVB_v − T_max,REF  (°C)",
                  title   = "HOBO validation — LOVB",
                  out_path = here::here("outputs/lovb/figures/hobo/fig_HOBO_LOVB_spearman_dTmax.png"))
  fwrite(ann_lovb[, .(Variable, rho_S = rho, p_value = p,
                        ci_lo, ci_hi, x_metric = "Tmicro-Tmacro")],
          here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))

  # --- LVA --- (NB : new sign convention applied to LVA too)
  long_lva <- make_long(DT, "LVA")
  ann_lva  <- make_ann(long_lva, vars)
  cli_alert("LVA Spearman :"); print(ann_lva[, .(Variable, rho, p, ci_lo, ci_hi)])
  build_hobo_fig(long_lva, vars, ann_lva,
                  x_label = x_label_signed,
                  y_label = "Simulated  T_max,NULL − T_max,LVA_v  (°C)",
                  title   = "HOBO validation — LVA",
                  out_path = here::here("outputs/lovb/figures/hobo/fig_HOBO_LVA_spearman_dTmax.png"))

  # --- Shapley (Owen approx) ---
  long_sh <- make_long(DT, "Shapley")
  ann_sh  <- make_ann(long_sh, vars)
  cli_alert("Shapley Spearman :"); print(ann_sh[, .(Variable, rho, p, ci_lo, ci_hi)])
  build_hobo_fig(long_sh, vars, ann_sh,
                  x_label = x_label_signed,
                  y_label = "Simulated  phi_v Shapley  (°C)",
                  title   = "HOBO validation — Shapley (Owen approx)",
                  out_path = here::here("outputs/lovb/figures/hobo/fig_HOBO_Shapley_spearman_dTmax.png"),
                  caption_extra = paste0(
                    "phi_v approximated as (LVA_v + LOVB_v) / 2 per sensor : ",
                    "averages position-1 (LVA singleton) and position-4 (LOVB) marginal contributions ; ",
                    "positions 2-3 (size-2 coalitions) not available at HOBO sensors (10/16 coalitions cached)."
                  ))

  cli_alert_success("All sensor-level figures regenerated with new sign convention.")
}

if (sys.nframe() == 0) run_all()
