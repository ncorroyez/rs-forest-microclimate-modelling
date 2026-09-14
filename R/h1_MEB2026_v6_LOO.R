# ==============================================================================
# MEB 2026 v6 : rename LOVB -> LOO in 3 figures (display only, data unchanged).
# Regenerate fig_typeB_final.png, fig_HOBO_final.png, fig_heatmap_MAE_final.png
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

# Display label : LOO instead of LOVB
method_yexpr <- function(m) switch(m,
  "LOVB"    = bquote(Delta[v]^"LOO" ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA" ~ "(°C)"),
  "Shapley" = bquote(varphi[v]       ~ "(°C)"))

# ============================================================================
# Type B (LOO label on row 1)
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

compute_bounds <- function() {
  longs <- list()
  for (m in c("LOVB","LVA","Shapley")) {
    longs[[m]] <- build_long_typeB(m)
    longs[[m]][, method := m]
  }
  all_long <- rbindlist(longs, fill = TRUE)
  x_bounds <- all_long[!is.na(trait), .(xmin = min(trait), xmax = max(trait)),
                          by = variable]
  x_bounds[, xrange := xmax - xmin]
  x_bounds[, xmin_pad := xmin - 0.03 * xrange]
  x_bounds[, xmax_pad := xmax + 0.03 * xrange]
  y_bounds <- all_long[!is.na(delta), .(ymin = min(delta), ymax = max(delta)),
                          by = method]
  y_bounds[, ymin_pad := pmin(ymin, 0) * 1.05]
  y_bounds[, ymax_pad := pmax(ymax, 0) * 1.05]
  list(x = x_bounds, y = y_bounds)
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
  cli_h1("Type B v6 (LOO label)")
  bds <- compute_bounds()
  rows <- list()
  for (m in c("LOVB","LVA","Shapley")) {
    yb <- bds$y[method == m]
    y_lim <- c(yb$ymin_pad, yb$ymax_pad)
    rows[[m]] <- make_typeB_row(m, bds$x, y_lim,
                                    show_x_axis = (m == "Shapley"),
                                    show_col_titles = (m == "LOVB"))
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
# HOBO with LOO strip
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
  cli_h1("HOBO v6 (LOO strip)")
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
  # Use LOO label as factor level for display
  long[, method_disp := ifelse(method == "LOVB", "LOO", method)]
  long[, method_disp := factor(method_disp, levels = c("LOO","LVA","Shapley"))]
  long[, Variable := factor(Variable, levels = VARS_HOBO)]
  long <- long[!is.na(y) & !is.na(dT)]

  ann <- long[, {r <- spearman_boot_ci(dT, y);
                  .(rho = r$rho, p = r$p, ci_lo = r$ci_lo, ci_hi = r$ci_hi)},
              by = .(method_disp, Variable)]
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
    facet_grid(method_disp ~ Variable, scales = "free", switch = "y") +
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

# ============================================================================
# Heatmap with LOO column labels
# ============================================================================
build_heatmap <- function() {
  cli_h1("Heatmap v6 (LOO labels)")
  arch <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
  clhs <- fread(file.path(TAB, "tab_attribution_cLHS_mean_floor05.csv"))
  arch[, variable := factor(variable, levels = VARS)]
  clhs[, variable := factor(variable, levels = VARS)]
  methods <- c("LOVB_arch","LVA_arch","Shapley_arch",
                "LOVB_cLHS","LVA_cLHS","Shapley_cLHS")
  labels  <- c("LOO","LVA","Shapley","LOO","LVA","Shapley")    # <-- LOO instead of LOVB
  long_list <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    m <- methods[i]
    src <- if (grepl("_arch$", m)) arch else clhs
    long_list[[i]] <- src[, .(variable, method = labels[i], col_idx = i,
                                Rank = get(paste0(m, "_rank")))]
  }
  long <- rbindlist(long_list)
  cons <- long[, .(Rank = mean(Rank), col_idx = 7L), by = variable]
  cons[, method := "Averaged Score"]
  long_full <- rbind(long[, .(variable, method, col_idx, Rank)],
                      cons[, .(variable, method, col_idx, Rank)])
  long_full[, fill_rank := as.character(round(Rank))]
  pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
  x_labels <- c("LOO","LVA","Shapley","LOO","LVA","Shapley","Averaged\nScore")    # <-- LOO
  top_labels <- data.table(x = c(2, 5, 7),
                              label = c("n = 4 archetypes", "n = 400 cLHS", ""))
  sep_xs <- c(3.5, 6.5)

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.5) +
    geom_text(aes(label = ifelse(method == "Averaged Score",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                                 "white", "grey15"),
               fontface = "bold", size = 7) +
    geom_vline(xintercept = sep_xs, colour = "white", linewidth = 5) +
    geom_vline(xintercept = sep_xs, colour = "grey30", linewidth = 1.2) +
    geom_text(data = top_labels[label != ""],
               aes(x = x, y = 4.8, label = label),
               inherit.aes = FALSE, fontface = "italic", size = 4,
               colour = "grey30") +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    coord_cartesian(ylim = c(0.5, 5.05), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(face = "bold", size = 13),
           axis.text.y = element_text(face = "bold", size = 15),
           panel.grid = element_blank(),
           plot.margin = margin(25, 10, 10, 10),
           legend.position = "right",
           legend.title = element_text(face = "bold"))
  fig_path <- file.path(OUT, "fig_heatmap_MAE_final.png")
  ggsave(fig_path, p, width = 8, height = 8, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_typeB()
  build_hobo()
  build_heatmap()
  cli_alert_success("v6 done : LOVB -> LOO display in 3 figures")
}
