# ==============================================================================
# Silvilaser v3bis — 4 additional figures completing RESULT 2 (archetypes) and
# RESULT 3 (HOBO) for LVA and Shapley.
#
# Sign conventions (verified in R/lovb_03_contrib.R, unchanged) :
#   LOVB    Delta_v = Tmax_REF - Tmax_LOVB_v  (positive = removing v reduces buffering)
#   LVA     Delta_v = Tmax_LVA_v - Tmax_NULL  (positive = adding v alone increases buffering)
#   Shapley phi_v   = exact average over 16 coalitions (positive = v contributes)
#
# Fig 4 caveat : per-sensor exact Shapley needs all 16 coalitions ; HOBO has
# only 10 cached. We use the Owen-style approximation phi_v ≈ (LVA + LOVB)/2
# which corresponds to averaging position-1 and position-4 marginals from
# the 24 permutations (positions 2-3 missing). Documented in figure caption.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

OUT_DIR <- here::here("outputs/figs_silvilaser_v3bis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Color palette matching existing LOVB archetype figure
.VAR_ORDER    <- c("LAI", "fCover", "Hmax", "LAD")
.PAL_VAR      <- c(LAI = "#440154", fCover = "#31688e",
                    Hmax = "#35b779", LAD = "#fde725")
.CLUSTER_LABS <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Bas couvert", "Arch_C4" = "C4 Inter")

# ---- Figure 1 : LVA contribution by archetype --------------------------------
build_fig_arch_LVA <- function() {
  cli_h1("Figure 1 — LVA contribution by archetype (MEAN)")
  DT_a <- readRDS(here::here("outputs/lovb/data/DT_contrib_archetypes.rds"))

  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lva_col <- paste0("LVA_", v, "_mean")
    data.table(archetype = DT_a$archetype,
                variable  = v,
                LVA_v     = DT_a[[lva_col]])
  }))
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, archetype_lab := .CLUSTER_LABS[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_LABS)]

  cli_alert("LVA by archetype :")
  print(dcast(long, archetype_lab ~ variable, value.var = "LVA_v"))

  ymax <- max(abs(long$LVA_v), na.rm = TRUE) * 1.2

  p <- ggplot(long, aes(x = variable, y = LVA_v, fill = variable)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_col(colour = "grey20", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", LVA_v),
                    vjust = ifelse(LVA_v >= 0, -0.4, 1.4)),
               fontface = "bold", size = 3.6) +
    facet_wrap(~ archetype_lab, nrow = 1L) +
    scale_fill_manual(values = .PAL_VAR, guide = "none") +
    scale_y_continuous(limits = c(-ymax, ymax),
                        expand = expansion(mult = c(0.06, 0.06))) +
    labs(
      title    = "LVA contribution by archetype  (MEAN)",
      subtitle = "Delta_v = Tmax_LVA_v_mean - Tmax_NULL_mean  |  positive = activating v alone increases buffering vs null canopy",
      x = NULL,
      y = "Delta_v  Tmax (summer mean)  [degree C]",
      caption = "Bit code NULL=0000  |  LVA_LAI=1000  LVA_Hmax=0100  LVA_fCover=0010  LVA_LAD=0001"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 12),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 1),
           panel.grid.minor = element_blank())

  fig_path <- file.path(OUT_DIR, "fig_archetypes_LVA_mean.png")
  ggsave(fig_path, p, width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Figure 2 : Shapley contribution by archetype ----------------------------
build_fig_arch_Shapley <- function() {
  cli_h1("Figure 2 — Shapley contribution by archetype (MEAN)")
  sh <- fread(here::here("outputs/h1/shapley_archetypes_values.csv"))
  sh <- sh[archetype %in% c("Arch_C1", "Arch_C2", "Arch_C3", "Arch_C4")]
  sh[, archetype_lab := .CLUSTER_LABS[archetype]]
  sh[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_LABS)]
  sh[, variable := factor(variable, levels = .VAR_ORDER)]

  cli_alert("Shapley phi by archetype :")
  print(dcast(sh, archetype_lab ~ variable, value.var = "phi"))

  ymax <- max(abs(c(sh$phi, sh$ci95_lo, sh$ci95_hi)), na.rm = TRUE) * 1.15

  p <- ggplot(sh, aes(x = variable, y = phi, fill = variable)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_col(colour = "grey20", linewidth = 0.3, width = 0.7) +
    geom_errorbar(aes(ymin = ci95_lo, ymax = ci95_hi),
                    width = 0.18, colour = "grey25", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.2f", phi),
                    vjust = ifelse(phi >= 0, -0.5, 1.6)),
               fontface = "bold", size = 3.6) +
    facet_wrap(~ archetype_lab, nrow = 1L) +
    scale_fill_manual(values = .PAL_VAR, guide = "none") +
    scale_y_continuous(limits = c(-ymax, ymax),
                        expand = expansion(mult = c(0.06, 0.06))) +
    labs(
      title    = "Shapley contribution by archetype  (MEAN)",
      subtitle = "phi_v = average marginal contribution across 16 coalitions  |  positive = v contributes to buffering",
      x = NULL,
      y = "phi_v   Tmax (summer mean)  [degree C]",
      caption = "Error bars : bootstrap CI95 (n=200 reps). Exact 2^4 enumeration ; 24 permutations per archetype."
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 12),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 1),
           panel.grid.minor = element_blank())

  fig_path <- file.path(OUT_DIR, "fig_archetypes_Shapley_mean.png")
  ggsave(fig_path, p, width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Figures 3 & 4 : HOBO validation LVA + Shapley ---------------------------
# Helper : Spearman + bootstrap CI95 (mirror previous logic)
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

# Generic HOBO figure builder
build_hobo_fig <- function(DT_pairs, y_var, y_label, title_tag, fig_name,
                              vars = c("LAI", "Hmax", "fCover", "LAD"),
                              point_colour = "#2C5F2D",
                              loess_colour = "#FFB400",
                              caption_extra = "") {
  long <- rbindlist(lapply(seq_along(vars), function(i)
    data.table(Variable = vars[i],
                dTmax_obs = DT_pairs$dTmax_obs,
                y         = DT_pairs[[paste0(y_var, "_", vars[i])]])))
  long[, Variable := factor(Variable, levels = vars)]

  ann_rows <- vector("list", length(vars))
  for (i in seq_along(vars)) {
    sub <- long[Variable == vars[i]]
    r <- spearman_boot_ci(sub$dTmax_obs, sub$y, B = 1000L)
    ann_rows[[i]] <- data.table(Variable = vars[i], rho = r$rho, p = r$p,
                                  ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  }
  ann <- rbindlist(ann_rows)
  ann[, label := sprintf("rho = %+.2f  p = %s\n[%+.2f ; %+.2f]",
                           rho,
                           ifelse(p < 1e-3, "<0.001", sprintf("%.3f", p)),
                           ci_lo, ci_hi)]
  ann[, Variable := factor(Variable, levels = vars)]
  cli_alert("Spearman per variable :")
  print(ann)

  p <- ggplot(long, aes(x = dTmax_obs, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = point_colour, alpha = 0.7, size = 1.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = loess_colour, fill = loess_colour,
                  linewidth = 0.8, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.2, size = 3.1, fontface = "bold",
               inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = paste0("HOBO validation (dTmax_obs) — ", title_tag),
      subtitle = "n = 53 sensors  |  Spearman rho + bootstrap CI95 (1000 reps)  |  LOESS span = 0.9",
      x = "Observed buffering mean Tmax_macro - Tmax_obs (K)",
      y = y_label,
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
  # HOBO figs : route to outputs/lovb/figures/hobo/ ; arch figs stay in OUT_DIR
  fig_path <- if (grepl("^hobo/", fig_name))
                  here::here("outputs/lovb/figures", fig_name)
                else file.path(OUT_DIR, fig_name)
  dir.create(dirname(fig_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(fig_path, p, width = 11.4, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Build the HOBO dataset (sensor-level pairs) -----------------------------
build_hobo_pairs <- function() {
  DT_h    <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  cg      <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo_t  <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
  DT <- merge(DT_h, hobo_t, by = "id_plot")

  # LVA per sensor and variable : Tmax_LVA_v - Tmax_NULL
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  for (v in vars) {
    DT[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) -
                                  get("Tmax_mean_NULL")]
  }
  # LOVB per sensor : Delta_v already available as Delta_LAI_mean etc.
  # Owen-style Shapley approximation per sensor : (LVA + LOVB) / 2
  for (v in vars) {
    DT[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) +
                                     get(paste0("Delta_", v, "_mean"))) / 2]
    # Also alias LOVB consistently
    DT[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
  }
  DT[]
}

build_fig_HOBO_LVA <- function() {
  cli_h1("Figure 3 — HOBO validation LVA (dTmax_obs)")
  DT <- build_hobo_pairs()
  build_hobo_fig(DT, y_var = "LVA",
                  y_label = "Simulated LVA delta (K, mean over 122 d)",
                  title_tag = "LVA",
                  fig_name = "hobo/fig_HOBO_LVA_spearman_dTmax.png")
}

build_fig_HOBO_Shapley <- function() {
  cli_h1("Figure 4 — HOBO validation Shapley (dTmax_obs)")
  DT <- build_hobo_pairs()
  build_hobo_fig(DT, y_var = "Shapley",
                  y_label = "Simulated phi_v Shapley (K, mean over 122 d)",
                  title_tag = "Shapley (Owen approx)",
                  fig_name = "hobo/fig_HOBO_Shapley_spearman_dTmax.png",
                  caption_extra = paste0(
                    "phi_v approximated as (LVA_v + LOVB_v) / 2 per sensor : ",
                    "averages position-1 (LVA singleton) and position-4 (LOVB) marginal contributions ; ",
                    "positions 2-3 (size-2 coalitions) not available at HOBO sensors (10/16 coalitions cached)."
                  ))
}

if (sys.nframe() == 0) {
  build_fig_arch_LVA()
  build_fig_arch_Shapley()
  build_fig_HOBO_LVA()
  build_fig_HOBO_Shapley()
  cli_alert_success("All 4 figures saved in {.path {OUT_DIR}}")
}
