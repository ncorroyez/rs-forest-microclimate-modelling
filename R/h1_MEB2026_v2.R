# ==============================================================================
# MEB 2026 final v2 : minimal clean figures.
# Outputs (overwrite previous):
#   fig_typeB_final.png      : 3 rows (LOVB/LVA/Shapley) × 4 cols (LAI/fCover/Hmax/LAD violin)
#   fig_HOBO_final.png       : 3 rows × 4 cols (LAI/Hmax/fCover/LAD scatter)
#   fig_heatmap_MAE_final.png: 7 cols heatmap (no Level annotations, "Averaged Score" replaces "Consensus")
#   tab_MAE_values.csv       : kept
# Conventions: no title/subtitle, clean axes only. ΔTmax symbol for HOBO x.
# Cluster labels : C1, C2, C3, C4 (no descriptors).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
OUT  <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

VARS         <- c("LAI","fCover","Hmax","LAD")
VARS_HOBO    <- c("LAI","Hmax","fCover","LAD")
.CLUSTER_SHORT <- c("1"="C1","2"="C2","3"="C3","4"="C4")
.PAL_CLUSTER   <- c("C1"="#E69F00","C2"="#0072B2","C3"="#009E73","C4"="#CC79A7")

method_yexpr <- function(m) switch(m,
  "LOVB"    = bquote(Delta[v]^"LOVB"    ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA"     ~ "(°C)"),
  "Shapley" = bquote(varphi[v]          ~ "(°C)"))

# ============================================================================
# Type B merged figure : 3 rows × 4 cols (LAD = violin)
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

make_typeB_row <- function(method, show_x_axis = TRUE, show_col_titles = FALSE) {
  long <- build_long_typeB(method)
  yr <- range(long$delta, na.rm = TRUE)
  y_lim <- c(min(yr[1], 0)*1.05, max(yr[2], 0)*1.05)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)

  scatter <- function(v) {
    sub <- long[variable == v]
    ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.6) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      scale_y_continuous(limits = y_lim) +
      labs(x = if (show_x_axis) v else NULL,
            y = if (v == "LAI") method_yexpr(method) else NULL,
            title = if (show_col_titles) v else NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
             legend.position = "none",
             axis.text.y = if (v == "LAI") element_text() else element_blank(),
             axis.ticks.y = if (v == "LAI") element_line() else element_blank(),
             panel.grid.minor = element_blank())
  }
  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_violin(alpha = 0.45, colour = "grey30", scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.5, size = 0.9) +
    geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                  fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = y_lim) +
    labs(x = if (show_x_axis) "Cluster" else NULL, y = NULL,
          title = if (show_col_titles) "LAD" else NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
           legend.position = "none",
           axis.text.x = if (show_x_axis)
                            element_text(angle = 0, hjust = 0.5, size = 9, face = "bold")
                          else element_blank(),
           axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
           axis.text.y = element_blank(), axis.ticks.y = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())
  (scatter("LAI") + scatter("fCover") + scatter("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

build_typeB_final <- function() {
  cli_h1("Type B final (merged 3 rows × 4 cols)")
  rows <- list(
    LOVB    = make_typeB_row("LOVB",    show_x_axis = FALSE, show_col_titles = TRUE),
    LVA     = make_typeB_row("LVA",     show_x_axis = FALSE, show_col_titles = FALSE),
    Shapley = make_typeB_row("Shapley", show_x_axis = TRUE,  show_col_titles = FALSE)
  )
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_SHORT, levels=.CLUSTER_SHORT)),
            aes(x, x, colour = c)) + geom_point(size = 3) +
      scale_colour_manual(values = pal, name = NULL) +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom"))

  stacked <- (rows$LOVB / rows$LVA / rows$Shapley)
  final <- cowplot::plot_grid(stacked, leg, ncol = 1, rel_heights = c(1, 0.04))
  fig_path <- file.path(OUT, "fig_typeB_final.png")
  ggsave(fig_path, final, width = 14, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# HOBO merged figure : facet_grid(method ~ variable)
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

build_hobo_final <- function() {
  cli_h1("HOBO final (merged 3 rows × 4 cols)")
  DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT_h <- merge(DT_h, hobo_t, by = "id_plot")
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
  long[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
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
  cli_alert("Spearman (n=53) :"); print(ann[, .(method, Variable, rho = round(rho,3),
                                                    p = signif(p,3), ci_lo = round(ci_lo,3),
                                                    ci_hi = round(ci_hi,3))])

  p <- ggplot(long, aes(x = dT, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.6) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.7, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = rho_lab),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = p_lab),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.5,
               inherit.aes = FALSE, colour = "grey25") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = ci_lab),
               parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.1,
               inherit.aes = FALSE, colour = "grey35") +
    facet_grid(method ~ Variable, scales = "free_y", switch = "y") +
    labs(x = bquote(Delta * T[max] ~ "(°C)"),
          y = bquote("Simulated  " * Delta[v] ~ "(°C)")) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 13),
           strip.placement = "outside",
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_HOBO_final.png")
  ggsave(fig_path, p, width = 16, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# Heatmap MAE — minimal, no Level annotations, "Averaged Score" instead of Consensus
# ============================================================================
build_heatmap_final <- function() {
  cli_h1("Heatmap MAE final (minimal)")
  arch <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
  clhs <- fread(file.path(TAB, "tab_attribution_cLHS_mean_floor05.csv"))
  arch[, variable := factor(variable, levels = VARS)]
  clhs[, variable := factor(variable, levels = VARS)]

  methods <- c("LOVB_arch","LVA_arch","Shapley_arch",
                "LOVB_cLHS","LVA_cLHS","Shapley_cLHS")
  labels  <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley")
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
  x_labels <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley","Averaged\nScore")
  sep_xs <- c(3.5, 6.5)

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(method == "Averaged Score",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                                 "white", "grey15"),
               fontface = "bold", size = 5.5) +
    geom_vline(xintercept = sep_xs, colour = "white", linewidth = 3) +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(face = "bold", size = 11),
           axis.text.y = element_text(face = "bold", size = 13),
           panel.grid = element_blank(),
           plot.margin = margin(10, 10, 10, 10))
  fig_path <- file.path(OUT, "fig_heatmap_MAE_final.png")
  ggsave(fig_path, p, width = 11, height = 4, dpi = 150)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  # Remove old individual figures to avoid confusion
  for (m in c("LOVB","LVA","Shapley")) {
    for (kind in c("typeB","HOBO")) {
      f <- file.path(OUT, sprintf("fig_%s_%s_final.png", kind, m))
      if (file.exists(f)) file.remove(f)
    }
  }
  build_typeB_final()
  build_hobo_final()
  build_heatmap_final()
  cli_alert_success("3 final figures in {.path {OUT}}")
}
