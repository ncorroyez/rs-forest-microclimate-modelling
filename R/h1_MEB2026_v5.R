# ==============================================================================
# MEB 2026 v5 : enforce consistent axis bounds.
#   Type B : per-variable x-bounds (same across 3 methods), per-method y-bounds
#            (same across 4 variables).
#   HOBO   : facet_grid scales = "free" → x free per col, y free per row (auto)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
OUT  <- here::here("outputs/figs_MEB2026_final")

VARS         <- c("LAI","fCover","Hmax","LAD")
.CLUSTER_SHORT <- c("1"="C1","2"="C2","3"="C3","4"="C4")
.PAL_CLUSTER   <- c("C1"="#E69F00","C2"="#0072B2","C3"="#009E73","C4"="#CC79A7")

method_yexpr <- function(m) switch(m,
  "LOVB"    = bquote(Delta[v]^"LOVB"    ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA"     ~ "(°C)"),
  "Shapley" = bquote(varphi[v]          ~ "(°C)"))

# ============================================================================
# Type B with explicit per-variable x_lim AND per-method y_lim
# ============================================================================
build_long_typeB <- function(method) {
  DT <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
  tm <- c(LAI="LAI", fCover="fCover", Hmax="Hmax", LAD="FPC1")
  if (method == "LOVB") {
    rows <- lapply(VARS, function(v) data.table(
      Cluster=DT$Cluster, variable=v,
      trait=DT[[tm[v]]], delta=DT[[paste0("Delta_",v,"_mean")]]))
  } else if (method == "LVA") {
    rows <- lapply(VARS, function(v) data.table(
      Cluster=DT$Cluster, variable=v,
      trait=DT[[tm[v]]], delta=DT[[paste0("LVA_",v,"_mean")]]))
  } else {
    DT_phi <- readRDS(file.path(DATA, "DT_shapley_per_plot_floor05.rds"))
    DT_phi[, trait := NA_real_]
    for (v in VARS) DT_phi[variable == v, trait := get(tm[v])]
    rows <- list(DT_phi[, .(Cluster, variable, trait, delta = phi)])
  }
  long <- rbindlist(rows)
  long[, variable := factor(variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_SHORT[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_SHORT)]
  long
}

# Compute per-variable x ranges and per-method y ranges (across all methods/vars)
compute_bounds <- function() {
  longs <- list()
  for (m in c("LOVB","LVA","Shapley")) {
    longs[[m]] <- build_long_typeB(m)
    longs[[m]][, method := m]
  }
  all_long <- rbindlist(longs, fill = TRUE)
  # Per-variable x bounds
  x_bounds <- all_long[!is.na(trait), .(xmin = min(trait), xmax = max(trait)),
                          by = variable]
  x_bounds[, xrange := xmax - xmin]
  x_bounds[, xmin_pad := xmin - 0.03 * xrange]
  x_bounds[, xmax_pad := xmax + 0.03 * xrange]
  # Per-method y bounds (across all 4 vars)
  y_bounds <- all_long[!is.na(delta), .(ymin = min(delta), ymax = max(delta)),
                          by = method]
  y_bounds[, ymin_pad := pmin(ymin, 0) * 1.05]
  y_bounds[, ymax_pad := pmax(ymax, 0) * 1.05]
  list(x = x_bounds, y = y_bounds, all = all_long)
}

make_typeB_row <- function(method, x_bounds, y_lim,
                              show_x_axis = TRUE, show_col_titles = FALSE) {
  long <- build_long_typeB(method)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)

  scatter <- function(v) {
    sub <- long[variable == v]
    xb  <- x_bounds[variable == v]
    ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.8) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      coord_cartesian(xlim = c(xb$xmin_pad, xb$xmax_pad), ylim = y_lim) +
      labs(x = if (show_x_axis) v else NULL,
            y = if (v == "LAI") method_yexpr(method) else NULL,
            title = if (show_col_titles) v else NULL) +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
             legend.position = "none",
             axis.title = element_text(size = 14),
             axis.text  = element_text(size = 12),
             axis.text.y = if (v == "LAI") element_text() else element_blank(),
             axis.ticks.y = if (v == "LAI") element_line() else element_blank(),
             panel.grid.minor = element_blank())
  }
  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(alpha = 0.45, colour = "grey30", scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.5, size = 1.0) +
    geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                  fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    coord_cartesian(ylim = y_lim) +
    labs(x = if (show_x_axis) "Cluster" else NULL, y = NULL,
          title = if (show_col_titles) "LAD" else NULL) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
           legend.position = "none",
           axis.text.x = if (show_x_axis) element_text(size = 12, face = "bold")
                          else element_blank(),
           axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
           axis.text.y = element_blank(), axis.ticks.y = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())
  (scatter("LAI") + scatter("fCover") + scatter("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

build_typeB <- function() {
  cli_h1("Type B with consistent bounds")
  bds <- compute_bounds()
  cli_alert("Per-variable x bounds :"); print(bds$x[, .(variable, xmin = round(xmin,2), xmax = round(xmax,2))])
  cli_alert("Per-method y bounds :"); print(bds$y[, .(method, ymin = round(ymin,2), ymax = round(ymax,2))])

  rows <- list()
  for (m in c("LOVB","LVA","Shapley")) {
    yb <- bds$y[method == m]
    y_lim <- c(yb$ymin_pad, yb$ymax_pad)
    is_last <- (m == "Shapley")
    is_first <- (m == "LOVB")
    rows[[m]] <- make_typeB_row(m, bds$x, y_lim,
                                    show_x_axis = is_last, show_col_titles = is_first)
  }
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_SHORT, levels=.CLUSTER_SHORT)),
            aes(x, x, colour = c)) + geom_point(size = 3.5) +
      scale_colour_manual(values = pal, name = NULL) +
      theme_bw(base_size = 15) + theme(legend.position = "bottom",
                                           legend.text = element_text(size = 13)))
  stacked <- (rows$LOVB / rows$LVA / rows$Shapley)
  final <- cowplot::plot_grid(stacked, leg, ncol = 1, rel_heights = c(1, 0.04))
  fig_path <- file.path(OUT, "fig_typeB_final.png")
  ggsave(fig_path, final, width = 14, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# HOBO with facet_grid scales = "free" (per-col x, per-row y, automatic)
# ============================================================================
spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 3) return(list(rho = NA, p = NA, ci_lo = NA, ci_hi = NA))
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

build_hobo <- function() {
  cli_h1("HOBO with scales = 'free' (per-col x, per-row y)")
  DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT_h <- merge(DT_h, hobo_t, by = "id_plot")
  VARS_HOBO <- c("LAI","Hmax","fCover","LAD")
  for (v in VARS_HOBO) {
    DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
    DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
  }
  DT_shap <- readRDS(file.path(DATA, "DT_shapley_per_HOBO_floor05.rds"))
  shap_w <- dcast(DT_shap, id_plot ~ variable, value.var = "phi")
  setnames(shap_w, c("LAI","Hmax","fCover","LAD"),
            paste0("Shapley_", c("LAI","Hmax","fCover","LAD")))
  DT_h <- merge(DT_h, shap_w, by = "id_plot", all.x = TRUE)

  long <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
    rbindlist(lapply(VARS_HOBO, function(v) data.table(
      method = m, Variable = v,
      dT = DT_h$dT, y = DT_h[[paste0(m, "_", v)]])))))
  long[, method := factor(method, levels = c("LOVB","LVA","Shapley"))]
  long[, Variable := factor(Variable, levels = VARS_HOBO)]
  long <- long[!is.na(y) & !is.na(dT)]

  ann <- long[, {r <- spearman_boot_ci(dT, y);
                  .(rho = r$rho, p = r$p, ci_lo = r$ci_lo, ci_hi = r$ci_hi)},
              by = .(method, Variable)]
  ann[, rho_lab := sprintf("rho == %+.2f", rho)]
  ann[, p_lab   := ifelse(p < 1e-3, "italic(p) < 0.001",
                            sprintf("italic(p) == %.3f", p))]
  ann[, ci_lab  := sprintf("CI[95] *' '* '[' * %+.2f * '; ' * %+.2f * ']'",
                              ci_lo, ci_hi)]

  p <- ggplot(long, aes(x = dT, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.7, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = rho_lab),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = p_lab),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 4.2,
               inherit.aes = FALSE, colour = "grey25") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = ci_lab),
               parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.7,
               inherit.aes = FALSE, colour = "grey35") +
    # IMPORTANT : scales = "free" gives x per col (= same Variable → same x bounds),
    # y per row (= same method → same y bounds)
    facet_grid(method ~ Variable, scales = "free", switch = "y") +
    labs(x = bquote(Delta * T[max] ~ "(°C)"),
          y = bquote("Simulated  " * Delta[v] ~ "(°C)")) +
    theme_bw(base_size = 16) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 15),
           strip.placement = "outside",
           axis.title = element_text(size = 15),
           axis.text = element_text(size = 12),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_HOBO_final.png")
  ggsave(fig_path, p, width = 16, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_typeB()
  build_hobo()
  cli_alert_success("MEB v5 updates done")
}
