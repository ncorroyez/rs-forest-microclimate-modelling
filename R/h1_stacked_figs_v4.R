# ==============================================================================
# v3quater : facet-based stacked figures with proper symbols and bigger fonts.
#  - ΔTmax symbol on HOBO x-axis (instead of T_micro - T_macro raw)
#  - ρ rendered as Greek symbol in annotations (via parse=TRUE)
#  - φ Shapley renamed Δ_v^Shapley for cohérence visuelle
#  - facet_grid(method ~ variable) for HOBO and Type A : single col title
#  - Type B keeps patchwork (LAD violin needs separate geom), but with col
#    titles only at top row
#  - base_size = 14 throughout
#  - Type A converted from bars to SCATTER (Tmax_REF vs Tmax_method_v)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

OUT_DIR <- here::here("outputs/figs_silvilaser_v3quater")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.VAR_ORDER    <- c("LAI", "fCover", "Hmax", "LAD")
.PAL_VAR      <- c(LAI = "#440154", fCover = "#31688e",
                    Hmax = "#35b779", LAD = "#fde725")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")

# Symbol expressions
method_yexpr <- function(method) switch(method,
  "LOVB"    = bquote(Delta[v]^"LOVB" ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA"  ~ "(°C)"),
  "Shapley" = bquote(Delta[v]^"Shapley" ~ "(°C)")
)
method_full_label <- function(method) switch(method,
  "LOVB"    = "LOVB :  Δ_v^LOVB = T_max,LOVB_v − T_max,REF",
  "LVA"     = "LVA :  Δ_v^LVA = T_max,NULL − T_max,LVA_v",
  "Shapley" = "Shapley :  Δ_v^Shapley = exact attribution from 16 coalitions"
)

x_expr_dT <- bquote(Delta * T[max]^"obs" ~ "(°C)")
x_expr_TmaxREF <- bquote(T[max*","*REF] ~ "(°C)")

# ==============================================================================
# Type B  (per-plot Δ_v vs trait + LAD violin)
# ==============================================================================
build_long_typeB <- function(method) {
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")
  if (method == "LOVB") {
    rows <- lapply(.VAR_ORDER, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("Delta_", v, "_mean")]]))
  } else if (method == "LVA") {
    rows <- lapply(.VAR_ORDER, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("LVA_", v, "_mean")]]))
  } else if (method == "Shapley") {
    DT_phi <- readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
    DT_phi[, trait := NA_real_]
    for (v in .VAR_ORDER) DT_phi[variable == v, trait := get(trait_map[v])]
    rows <- list(DT_phi[, .(Cluster, variable, trait, delta = phi)])
  }
  long <- rbindlist(rows)
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  long[]
}

make_typeB_row <- function(method, show_x = TRUE, show_col_titles = FALSE) {
  long <- build_long_typeB(method)
  yr <- range(long$delta, na.rm = TRUE)
  y_lim <- c(min(yr[1], 0) * 1.05, max(yr[2], 0) * 1.05)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  scatter_panel <- function(v) {
    sub <- long[variable == v]
    ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.8) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      scale_y_continuous(limits = y_lim) +
      labs(x = if (show_x) paste(v) else NULL,
            y = if (v == "LAI") method_yexpr(method) else NULL,
            title = if (show_col_titles) v else NULL) +
      theme_bw(base_size = 14) +
      theme(plot.title       = element_text(face = "bold", hjust = 0.5),
             legend.position  = "none",
             panel.grid.minor = element_blank(),
             axis.text.y      = if (v == "LAI") element_text() else element_blank(),
             axis.ticks.y     = if (v == "LAI") element_line() else element_blank(),
             axis.title.y     = element_text(size = 12))
  }
  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(alpha = 0.45, colour = "grey30", scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.5,
                  size = 1.0) +
    geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                  fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = y_lim) +
    labs(x = if (show_x) "Cluster" else NULL,
          y = NULL,
          title = if (show_col_titles) "LAD (distrib.)" else NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5),
           legend.position  = "none",
           axis.text.x      = if (show_x)
                                  element_text(angle = 25, hjust = 1, size = 10)
                                else element_blank(),
           axis.ticks.x     = if (show_x) element_line() else element_blank(),
           axis.text.y      = element_blank(),
           axis.ticks.y     = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())

  (scatter_panel("LAI") + scatter_panel("fCover") +
   scatter_panel("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

build_typeB_facet <- function() {
  cli_h1("Type B facet (3 methods × 4 vars, col title only at top)")
  rows <- list(
    LOVB    = make_typeB_row("LOVB",    show_x = FALSE, show_col_titles = TRUE),
    LVA     = make_typeB_row("LVA",     show_x = FALSE, show_col_titles = FALSE),
    Shapley = make_typeB_row("Shapley", show_x = TRUE,  show_col_titles = FALSE)
  )
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  leg_plot <- ggplot(data.table(x=1:4, y=1:4, c=factor(.CLUSTER_LABS,
                                                         levels=.CLUSTER_LABS)),
                       aes(x, y, colour = c)) +
    geom_point(size = 3) + scale_colour_manual(values = pal, name = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
           legend.text = element_text(size = 13))
  legend <- cowplot::get_legend(leg_plot)
  stacked <- (rows$LOVB / rows$LVA / rows$Shapley) +
    plot_annotation(
      title    = bquote("Type B — Per-plot attribution"
                          ~ Delta[v] ~ "vs LiDAR traits (n = 400 cLHS)"),
      subtitle = "Rows : LOVB / LVA / Shapley     |     Cols 1-3 : scatter + LOESS span 0.75 with 95%CI     |     Col 4 : LAD distribution by archetype",
      caption  = "Sign convention : Δ_v positive  ⇒  v contributes to buffering (reduces sub-canopy warming).",
      theme = theme(plot.title    = element_text(face = "bold", size = 16),
                     plot.subtitle = element_text(colour = "grey25", size = 12),
                     plot.caption  = element_text(size = 10, colour = "grey35",
                                                    hjust = 0))
    )
  final <- cowplot::plot_grid(stacked, legend, ncol = 1,
                                rel_heights = c(1, 0.045))
  fig_path <- file.path(OUT_DIR, "fig_typeB_facet_v3quater.png")
  ggsave(fig_path, final, width = 16, height = 14, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ==============================================================================
# HOBO  (facet_grid method ~ variable)
# ==============================================================================
build_hobo_pairs <- function() {
  DT_h <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  cg   <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT <- merge(DT_h, hobo, by = "id_plot")
  for (v in .VAR_ORDER) {
    DT[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) -
                                  get("Tmax_mean_NULL")]
    DT[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
    DT[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) +
                                       get(paste0("LOVB_", v))) / 2]
  }
  long <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
    rbindlist(lapply(.VAR_ORDER, function(v) data.table(
      method = m, Variable = v,
      dT = DT$dT,
      y  = DT[[paste0(m, "_", v)]])))))
  long[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
  long[, Variable := factor(Variable, levels = .VAR_ORDER)]
  long[]
}

spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rho_hat <- suppressWarnings(cor(x, y, method = "spearman"))
  p_val   <- suppressWarnings(cor.test(x, y, method = "spearman",
                                          exact = FALSE)$p.value)
  set.seed(seed); rb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    rb[b] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }
  ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
  list(rho = rho_hat, p = p_val, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

build_hobo_facet <- function() {
  cli_h1("HOBO validation facet — 3 methods x 4 vars")
  long <- build_hobo_pairs()
  # Compute Spearman per (method, variable)
  ann <- long[, {
    r <- spearman_boot_ci(dT, y)
    .(rho = r$rho, p = r$p, ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  }, by = .(method, Variable)]
  # parseable label : Greek rho via expression
  ann[, rho_lab := sprintf("rho == %+.2f", rho)]
  ann[, p_lab   := ifelse(p < 1e-3, "italic(p) < 0.001",
                            sprintf("italic(p) == %.3f", p))]
  ann[, ci_lab  := sprintf("CI[95] *' '* '[' * %+.2f * '; ' * %+.2f * ']'",
                             ci_lo, ci_hi)]

  p <- ggplot(long, aes(x = dT, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.7, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = rho_lab),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = p_lab),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.8,
               inherit.aes = FALSE, colour = "grey25") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = ci_lab),
               parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.3,
               inherit.aes = FALSE, colour = "grey35") +
    facet_grid(method ~ Variable, scales = "free_y", switch = "y") +
    labs(
      title    = bquote("HOBO sensor-level validation (n = 53)"
                          ~ " — " ~ rho[Spearman] ~ "between observed"
                          ~ Delta * T[max] ~ "and simulated"
                          ~ Delta[v] ),
      subtitle = "Rows : LOVB / LVA / Shapley     |     LOESS span = 0.9     |     Spearman ρ + bootstrap CI95 (1000 reps)",
      x = x_expr_dT,
      y = bquote("Simulated  " * Delta[v] ~ "(°C)"),
      caption = "Sign convention : Δ_v positive  ⇒  v contributes to buffering. Positive ρ expected (well-buffered plots = negative ΔT_obs)."
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 13),
           strip.placement  = "outside",
           plot.title       = element_text(face = "bold", size = 16),
           plot.subtitle    = element_text(colour = "grey25", size = 12),
           plot.caption     = element_text(size = 10, colour = "grey35", hjust = 0),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT_DIR, "fig_HOBO_facet_v3quater.png")
  ggsave(fig_path, p, width = 16, height = 13, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ==============================================================================
# Type A — scatter (Tmax_REF vs Tmax_method_v) facet_grid
# ==============================================================================
build_typeA_facet <- function() {
  cli_h1("Type A scatter facet — Tmax_REF vs Tmax_method_v")
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  DT_phi <- readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))

  # phi long per plot, wide it to match DT_c
  phi_w <- dcast(DT_phi, x + y ~ variable, value.var = "phi",
                   fun.aggregate = mean)
  setnames(phi_w, c("LAI","Hmax","fCover","LAD"),
            c("phi_LAI","phi_Hmax","phi_fCover","phi_LAD"))
  DT_c <- merge(DT_c, phi_w, by = c("x","y"), all.x = TRUE)
  # Shapley imputed sim : Tmax_REF - phi_v  (positive phi = v contributes  ⇒
  # without v the buffering would be lower  ⇒  T_max_Shapley_v ≈ Tmax_REF - phi_v)
  for (v in .VAR_ORDER)
    DT_c[, paste0("Tmax_Shapley_", v) := get("Tmax_mean_REF") -
                                              get(paste0("phi_", v))]

  long <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m) {
    rbindlist(lapply(.VAR_ORDER, function(v) {
      ycol <- if (m == "LOVB")    paste0("Tmax_mean_LOVB_", v)
               else if (m == "LVA") paste0("Tmax_mean_LVA_", v)
               else                  paste0("Tmax_Shapley_", v)
      data.table(method = m, Variable = v,
                  Cluster = DT_c$Cluster,
                  Tmax_REF = DT_c$Tmax_mean_REF,
                  y        = DT_c[[ycol]])
    }))
  }))
  long[, method      := factor(method, levels = c("LOVB","LVA","Shapley"))]
  long[, Variable    := factor(Variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  # RMSE annotation per facet
  ann <- long[, .(rmse = round(sqrt(mean((y - Tmax_REF)^2, na.rm = TRUE)), 2)),
                by = .(method, Variable)]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = y, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 3.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_grid(method ~ Variable, scales = "free", switch = "y") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type A — Tmax simulated (method) vs Tmax simulated (REF)"
                          ~ ", n = 400 cLHS"),
      subtitle = "Rows : LOVB / LVA / Shapley (imputed = T_max,REF − Δ_v^Shapley)     |     Dashed line = 1:1 identity     |     Deviation = v's contribution to buffering",
      x = bquote(T[max] ~ "buffering REF (mean over 122 d)" ~ "(°C)"),
      y = bquote(T[max] ~ "buffering method v" ~ "(°C)"),
      caption = "RMSE per facet = magnitude of v's perturbation to the REF simulation."
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 13),
           strip.placement  = "outside",
           plot.title       = element_text(face = "bold", size = 16),
           plot.subtitle    = element_text(colour = "grey25", size = 12),
           plot.caption     = element_text(size = 10, colour = "grey35", hjust = 0),
           legend.position  = "bottom",
           legend.text      = element_text(size = 12),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT_DIR, "fig_typeA_facet_v3quater.png")
  ggsave(fig_path, p, width = 16, height = 13, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_typeB_facet()
  build_hobo_facet()
  build_typeA_facet()
  cli_alert_success("All v3quater figures in {.path {OUT_DIR}}")
}
