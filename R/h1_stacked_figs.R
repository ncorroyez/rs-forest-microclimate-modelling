# ==============================================================================
# Stacked figures + symbol harmonization for Silvilaser deck v3ter.
# Produces :
#   - 3 individual Type B figures (LOVB, LVA, Shapley) with proper Δ, ρ, φ
#   - 1 stacked Type B (3 methods x 4 cols, 12 panels total)
#   - 1 stacked HOBO validation (3 methods x 4 cols)
#   - 1 stacked Type A archetype bars (3 methods x 4 archetypes, 12 panels)
# All outputs in outputs/figs_silvilaser_v3ter/.
# Symbols : Δ, ρ, φ via plotmath expressions ; subscripts properly rendered.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

OUT_DIR <- here::here("outputs/figs_silvilaser_v3ter")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Palettes ----------------------------------------------------------------
.VAR_ORDER    <- c("LAI", "fCover", "Hmax", "LAD")
.PAL_VAR      <- c(LAI = "#440154", fCover = "#31688e",
                    Hmax = "#35b779", LAD = "#fde725")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")
.CLUSTER_ARCH <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Low cover", "Arch_C4" = "C4 Inter")

# Method symbol per row (plotmath expressions)
method_symbol <- function(method) switch(method,
  "LOVB"    = bquote(Delta[v]^LOVB),
  "LVA"     = bquote(Delta[v]^LVA),
  "Shapley" = bquote(varphi[v])
)
method_y_axis <- function(method) switch(method,
  "LOVB"    = bquote(Delta[v]^LOVB == T[max*","*LOVB[v]] - T[max*","*REF] ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^LVA == T[max*","*NULL] - T[max*","*LVA[v]] ~ "(°C)"),
  "Shapley" = bquote(varphi[v] ~ "Shapley (°C)")
)

# ---- Build long format data for a given method (Type B) ----------------------
build_long_typeB <- function(method) {
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  vars <- .VAR_ORDER
  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")
  if (method == "LOVB") {
    rows <- lapply(vars, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("Delta_", v, "_mean")]]
    ))
  } else if (method == "LVA") {
    rows <- lapply(vars, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("LVA_", v, "_mean")]]
    ))
  } else if (method == "Shapley") {
    DT_phi <- readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
    DT_phi[, trait := NA_real_]
    for (v in vars) DT_phi[variable == v, trait := get(trait_map[v])]
    rows <- list(DT_phi[, .(Cluster, variable, trait, delta = phi)])
  }
  long <- rbindlist(rows)
  long[, variable    := factor(variable, levels = vars)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  long[]
}

# ---- One row of 4 panels (Type B) --------------------------------------------
make_typeB_row <- function(method, y_label_expr,
                              y_lim = NULL, show_x_axis = TRUE,
                              show_legend = FALSE, row_title = NULL) {
  long <- build_long_typeB(method)
  if (is.null(y_lim)) {
    yr <- range(long$delta, na.rm = TRUE)
    y_lim <- c(min(yr[1], 0) * 1.05, max(yr[2], 0) * 1.05)
  }
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  make_scatter <- function(v, show_y_axis = TRUE) {
    sub <- long[variable == v]
    p <- ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                  linewidth = 0.5) +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.5) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      scale_y_continuous(limits = y_lim) +
      labs(x = if (show_x_axis) paste(v, "value") else NULL,
            y = if (show_y_axis) y_label_expr else NULL,
            title = if (!is.null(row_title) && v == "LAI") row_title else NULL) +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 11),
             legend.position  = "none",
             panel.grid.minor = element_blank(),
             axis.title.x     = element_text(size = 9),
             axis.title.y     = element_text(size = 9))
    if (!show_y_axis) p <- p + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank())
    p
  }

  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                linewidth = 0.5) +
    geom_violin(alpha = 0.45, colour = "grey30", linewidth = 0.4,
                  scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.5,
                  size = 0.7) +
    geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                  fill = "white", outlier.shape = NA, linewidth = 0.35) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = y_lim) +
    labs(x = if (show_x_axis) "Cluster" else NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "none",
           axis.text.x      = if (show_x_axis)
                                  element_text(angle = 25, hjust = 1, size = 7.5)
                                else element_blank(),
           axis.ticks.x     = if (show_x_axis) element_line() else element_blank(),
           axis.text.y      = element_blank(),
           axis.ticks.y     = element_blank(),
           axis.title.x     = element_text(size = 9),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())

  p1 <- make_scatter("LAI",    show_y_axis = TRUE)
  p2 <- make_scatter("fCover", show_y_axis = FALSE)
  p3 <- make_scatter("Hmax",   show_y_axis = FALSE)
  (p1 + p2 + p3 + p_lad) + plot_layout(nrow = 1, widths = c(1,1,1,1))
}

# ---- Type B : single-method figures (regenerate v3ter) ----------------------
build_typeB_single <- function(method) {
  cli_h1(sprintf("Type B single — %s", method))
  row_plot <- make_typeB_row(method,
                                y_label_expr = method_y_axis(method),
                                show_x_axis = TRUE)
  long <- build_long_typeB(method)
  lad_med <- median(long[variable == "LAD", delta], na.rm = TRUE)
  lad_rng <- range(long[variable == "LAD", delta], na.rm = TRUE)

  # legend
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  leg_plot <- ggplot(long[variable == "LAI"],
                       aes(x = trait, y = delta, colour = Cluster_lab)) +
    geom_point(size = 2) + scale_colour_manual(values = pal, name = NULL) +
    theme_bw() + theme(legend.position = "bottom",
                          legend.text = element_text(size = 10))
  legend <- cowplot::get_legend(leg_plot)

  title_sym <- bquote("Contribution" ~ .(method_symbol(method))
                        ~ "vs LiDAR trait (MEAN)" ~ "—" ~ .(method))
  sub_text <- "LAI / fCover / Hmax = scatter + LOESS span 0.75 with 95%CI    |    LAD = distribution by archetype (functional object, no scalar projection)"
  caption_text <- sprintf("LAD distribution : range [%.2f, %.2f] °C, median %+.2f °C. 4 panels share Y-axis for direct visual comparison.",
                            lad_rng[1], lad_rng[2], lad_med)

  combined <- row_plot +
    plot_annotation(
      title    = title_sym,
      subtitle = sub_text,
      caption  = caption_text,
      theme = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey25", size = 9),
                     plot.caption  = element_text(size = 8, colour = "grey35",
                                                    hjust = 0))
    )
  final <- cowplot::plot_grid(combined, legend, ncol = 1,
                                rel_heights = c(1, 0.06))
  fig_path <- file.path(OUT_DIR,
                          sprintf("fig_typeB_%s_v3ter.png", method))
  ggsave(fig_path, final, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ---- STACKED Type B (3 methods x 4 cols, 12 panels) -------------------------
build_typeB_stack <- function() {
  cli_h1("Stacked Type B — 3 methods")
  # Compute global Y limits per method (rows independent)
  rows_p <- list()
  for (i in seq_along(c("LOVB", "LVA", "Shapley"))) {
    method <- c("LOVB", "LVA", "Shapley")[i]
    is_last <- (method == "Shapley")
    rows_p[[method]] <- make_typeB_row(method,
                                            y_label_expr = method_y_axis(method),
                                            show_x_axis = is_last)
  }
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  leg_plot <- ggplot(data.table(x=1:4, y=1:4,
                                   c=factor(.CLUSTER_LABS, levels=.CLUSTER_LABS)),
                       aes(x, y, colour = c)) +
    geom_point(size = 3) + scale_colour_manual(values = pal, name = NULL) +
    theme_bw() + theme(legend.position = "bottom",
                          legend.text = element_text(size = 11))
  legend <- cowplot::get_legend(leg_plot)

  stacked <- (rows_p[["LOVB"]] / rows_p[["LVA"]] / rows_p[["Shapley"]]) +
    plot_annotation(
      title    = bquote("Type B" ~ "—" ~ "Per-plot attribution"
                          ~ Delta[v] * "/" * varphi[v]
                          ~ "vs LiDAR traits (n = 400 cLHS)"),
      subtitle = "Rows : LOVB / LVA / Shapley     |     Cols 1-3 : scatter Δ_v vs trait + LOESS     |     Col 4 : LAD distribution by archetype",
      caption  = "Sign convention : positive Δ / φ = v contributes to buffering (reduces sub-canopy warming). LAD column is a violin (distribution of Δ_LAD over 400 plots, grouped by archetype).",
      theme = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey25", size = 9),
                     plot.caption  = element_text(size = 8, colour = "grey35",
                                                    hjust = 0))
    )
  final <- cowplot::plot_grid(stacked, legend, ncol = 1,
                                rel_heights = c(1, 0.04))
  fig_path <- file.path(OUT_DIR, "fig_typeB_STACKED_v3ter.png")
  ggsave(fig_path, final, width = 14, height = 13, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ---- HOBO validation per method ----------------------------------------------
build_hobo_long <- function() {
  DT_h <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  cg   <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
  DT <- merge(DT_h, hobo, by = "id_plot")
  DT[, dT_signed := -dTmax_obs]
  vars <- .VAR_ORDER
  for (v in vars) {
    DT[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) -
                                  get("Tmax_mean_NULL")]
    DT[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
    DT[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) +
                                       get(paste0("LOVB_", v))) / 2]
  }
  DT[]
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

make_hobo_row <- function(method, y_label_expr, DT, show_x_axis = TRUE) {
  vars <- .VAR_ORDER
  long <- rbindlist(lapply(vars, function(v)
    data.table(Variable = v, dT = DT$dT_signed,
                y = DT[[paste0(method, "_", v)]])))
  long[, Variable := factor(Variable, levels = vars)]

  ann_rows <- vector("list", length(vars))
  for (i in seq_along(vars)) {
    sub <- long[Variable == vars[i]]
    r <- spearman_boot_ci(sub$dT, sub$y, B = 1000L)
    ann_rows[[i]] <- data.table(Variable = vars[i],
                                  rho = r$rho, p = r$p,
                                  ci_lo = r$ci_lo, ci_hi = r$ci_hi)
  }
  ann <- rbindlist(ann_rows)
  ann[, label := sprintf("rho = %+.2f  p = %s\n[%+.2f ; %+.2f]",
                           rho, ifelse(p < 1e-3, "<0.001",
                                          sprintf("%.3f", p)),
                           ci_lo, ci_hi)]
  ann[, Variable := factor(Variable, levels = vars)]

  make_panel <- function(v, show_y_axis = TRUE) {
    sub <- long[Variable == v]
    asub <- ann[Variable == v]
    p <- ggplot(sub, aes(x = dT, y = y)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
      geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                    colour = "#FFB400", fill = "#FFB400",
                    linewidth = 0.7, alpha = 0.18, span = 0.9) +
      geom_text(data = asub, aes(x = -Inf, y = Inf, label = label),
                 hjust = -0.05, vjust = 1.15, size = 2.8, fontface = "bold",
                 inherit.aes = FALSE, colour = "grey15") +
      labs(x = if (show_x_axis)
                  bquote(T[micro*","*max] - T[macro*","*max] ~ "(°C)") else NULL,
            y = if (show_y_axis) y_label_expr else NULL,
            title = v) +
      theme_bw(base_size = 10) +
      theme(plot.title       = element_text(face = "bold", size = 11),
             panel.grid.minor = element_blank(),
             axis.title.x     = element_text(size = 9),
             axis.title.y     = element_text(size = 9))
    if (!show_y_axis) p <- p + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank())
    p
  }
  ps <- lapply(seq_along(vars), function(i) make_panel(vars[i], i == 1L))
  (ps[[1]] + ps[[2]] + ps[[3]] + ps[[4]]) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

# ---- STACKED HOBO validation (3 methods x 4 cols) ----------------------------
build_HOBO_stack <- function() {
  cli_h1("Stacked HOBO validation — 3 methods")
  DT <- build_hobo_long()
  methods <- c("LOVB", "LVA", "Shapley")
  rows_p <- list()
  for (i in seq_along(methods)) {
    method <- methods[i]
    is_last <- (method == "Shapley")
    rows_p[[method]] <- make_hobo_row(method,
                                          y_label_expr = method_y_axis(method),
                                          DT = DT,
                                          show_x_axis = is_last)
  }
  stacked <- (rows_p[["LOVB"]] / rows_p[["LVA"]] / rows_p[["Shapley"]]) +
    plot_annotation(
      title    = bquote("HOBO sensor-level validation (n = 53)" ~ "—"
                          ~ rho[Spearman] ~ "vs observed"
                          ~ (T[micro*","*max] - T[macro*","*max])),
      subtitle = "Rows : LOVB / LVA / Shapley     |     LOESS span = 0.9     |     Spearman + bootstrap CI95 (1000 reps)",
      caption  = "Sign convention : positive Δ / φ = v contributes to buffering. Positive rho expected as well-buffered plots (negative x) correspond to small Δ_v while exposed plots show large Δ_v.",
      theme = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey25", size = 9),
                     plot.caption  = element_text(size = 8, colour = "grey35",
                                                    hjust = 0))
    )
  fig_path <- file.path(OUT_DIR, "fig_HOBO_STACKED_v3ter.png")
  ggsave(fig_path, stacked, width = 14, height = 12, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ---- STACKED Type A (3 methods x 4 archetypes x 4 vars) ----------------------
make_typeA_row <- function(method, y_label_expr, show_x_axis = TRUE) {
  DT_a <- readRDS(here::here("outputs/lovb/data/DT_contrib_archetypes.rds"))
  if (method == "LOVB") {
    long <- rbindlist(lapply(.VAR_ORDER, function(v) data.table(
      archetype = DT_a$archetype, variable = v,
      val = DT_a[[paste0("Delta_", v, "_mean")]])))
  } else if (method == "LVA") {
    long <- rbindlist(lapply(.VAR_ORDER, function(v) data.table(
      archetype = DT_a$archetype, variable = v,
      val = DT_a[[paste0("LVA_", v, "_mean")]])))
  } else if (method == "Shapley") {
    sh <- fread(here::here("outputs/h1/shapley_archetypes_values.csv"))
    sh <- sh[archetype %in% paste0("Arch_C", 1:4)]
    long <- data.table(archetype = sh$archetype, variable = sh$variable,
                         val = sh$phi)
  }
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  ymax <- max(abs(long$val), na.rm = TRUE) * 1.2

  p <- ggplot(long, aes(x = variable, y = val, fill = variable)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(colour = "grey20", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", val),
                    vjust = ifelse(val >= 0, -0.3, 1.3)),
               fontface = "bold", size = 2.7) +
    facet_wrap(~ archetype_lab, nrow = 1L) +
    scale_fill_manual(values = .PAL_VAR, guide = "none") +
    scale_y_continuous(limits = c(-ymax, ymax),
                        expand = expansion(mult = c(0.06, 0.06))) +
    labs(x = if (show_x_axis) NULL else NULL,
          y = y_label_expr) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 10),
           plot.title       = element_text(face = "bold"),
           axis.title.y     = element_text(size = 9),
           axis.text.x      = if (show_x_axis)
                                  element_text(size = 8.5, face = "bold")
                                else element_blank(),
           axis.ticks.x     = if (show_x_axis) element_line() else element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())
  p
}

build_TypeA_stack <- function() {
  cli_h1("Stacked Type A archetype bars — 3 methods")
  methods <- c("LOVB", "LVA", "Shapley")
  rows_p <- list()
  for (i in seq_along(methods)) {
    method <- methods[i]
    is_last <- (method == "Shapley")
    rows_p[[method]] <- make_typeA_row(method,
                                            y_label_expr = method_y_axis(method),
                                            show_x_axis = is_last)
  }
  stacked <- (rows_p[["LOVB"]] / rows_p[["LVA"]] / rows_p[["Shapley"]]) +
    plot_annotation(
      title    = bquote("Type A archetype attribution" ~ "—"
                          ~ "exact" ~ Delta[v] * "/" * varphi[v]
                          ~ "at 4 K-means centroids"),
      subtitle = "Rows : LOVB / LVA / Shapley     |     4 facets per row = 4 archetype centroids     |     4 bars per facet = LAI / fCover / Hmax / LAD",
      caption  = "Sign convention : positive Δ / φ = v contributes to buffering. Exact attribution from 16 coalitions per archetype (Shapley uses CI95 not shown here).",
      theme = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "grey25", size = 9),
                     plot.caption  = element_text(size = 8, colour = "grey35",
                                                    hjust = 0))
    )
  fig_path <- file.path(OUT_DIR, "fig_typeA_STACKED_v3ter.png")
  ggsave(fig_path, stacked, width = 14, height = 12, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  for (m in c("LOVB", "LVA", "Shapley")) build_typeB_single(m)
  build_typeB_stack()
  build_HOBO_stack()
  build_TypeA_stack()
  cli_alert_success("All figures saved in {.path {OUT_DIR}}")
}
