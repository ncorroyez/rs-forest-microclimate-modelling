# ==============================================================================
# MEB 2026 v4 :
#   (1) NEW : fig_archetypes_profiles.png : 4 archetypes × LAD profile + LAI/Hmax/fCover stats
#   (2) Heatmap : square format, stronger separators
#   (3) Type B : larger fonts + alternative version (axes inverted)
#   (4) HOBO : larger fonts + alternative version (axes inverted)
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
# (1) Archetype profiles : LAD profile (mean ± SD) + LAI/Hmax/fCover stats per cluster
# ============================================================================
build_archetype_profiles <- function() {
  cli_h1("Archetype profiles (4 clusters)")
  df_sample <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds")))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample))
    df_sample[, Cluster := Archetype]
  lad_cols <- grep("^LAD_Layer_", names(df_sample), value = TRUE)
  z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

  # Compute mean LAD profile + SD per cluster
  profiles <- list()
  for (cl in c("1","2","3","4")) {
    sub <- df_sample[as.character(Cluster) == cl, ..lad_cols]
    sub_mat <- as.matrix(sub); sub_mat[is.na(sub_mat)] <- 0
    profiles[[cl]] <- data.table(
      Cluster = paste0("C", cl),
      height = z_breaks,
      lad_mean = colMeans(sub_mat),
      lad_sd   = apply(sub_mat, 2, sd))
  }
  prof_long <- rbindlist(profiles)
  prof_long[, Cluster := factor(Cluster, levels = paste0("C", 1:4))]

  # Trait stats per cluster
  stats <- df_sample[, .(LAI_mean = mean(LAI, na.rm=TRUE), LAI_sd = sd(LAI, na.rm=TRUE),
                            Hmax_mean = mean(Hmax, na.rm=TRUE), Hmax_sd = sd(Hmax, na.rm=TRUE),
                            fCover_mean = mean(fCover, na.rm=TRUE), fCover_sd = sd(fCover, na.rm=TRUE)),
                        by = Cluster]
  stats[, Cluster := paste0("C", Cluster)]
  stats[, Cluster := factor(Cluster, levels = paste0("C", 1:4))]
  setorder(stats, Cluster)
  cli_alert("Cluster stats :"); print(stats)

  # Build annotation text per cluster
  stats[, ann := sprintf("LAI = %.2f ± %.2f\nHmax = %.1f ± %.1f m\nfCover = %.2f ± %.2f",
                            LAI_mean, LAI_sd, Hmax_mean, Hmax_sd, fCover_mean, fCover_sd)]

  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)
  max_lad <- max(prof_long$lad_mean + prof_long$lad_sd, na.rm = TRUE) * 1.05

  p <- ggplot(prof_long, aes(x = lad_mean, y = height, colour = Cluster, fill = Cluster)) +
    geom_ribbon(aes(xmin = pmax(0, lad_mean - lad_sd), xmax = lad_mean + lad_sd),
                  alpha = 0.25, colour = NA) +
    geom_path(linewidth = 1.2, orientation = "y") +
    geom_text(data = stats, aes(x = max_lad*0.55, y = 38, label = ann, colour = Cluster),
               inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4.2,
               fontface = "bold", lineheight = 1.15, show.legend = FALSE) +
    facet_wrap(~ Cluster, nrow = 1L) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(0, max_lad),
                        expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, 40),
                        breaks = c(0, 10, 20, 30, 40)) +
    labs(x = bquote("LAD" ~ "(m"^"2"~"m"^"-3"*")"),
          y = "Height (m)") +
    theme_bw(base_size = 15) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 14),
           legend.position = "none",
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_archetypes_profiles.png")
  ggsave(fig_path, p, width = 14, height = 6, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (2) Heatmap square, stronger separators
# ============================================================================
build_heatmap_square <- function() {
  cli_h1("Heatmap MAE square")
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
    # Stronger separators
    geom_vline(xintercept = sep_xs, colour = "white", linewidth = 5) +
    geom_vline(xintercept = sep_xs, colour = "grey30", linewidth = 1.2,
                 linetype = "solid") +
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
  # Square : 8 × 8 in
  ggsave(fig_path, p, width = 8, height = 8, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}} (square)")
}

# ============================================================================
# (3) Type B with larger fonts + inverted version
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

make_typeB_row <- function(method, show_x_axis = TRUE, show_col_titles = FALSE,
                              inverted = FALSE) {
  long <- build_long_typeB(method)
  yr <- range(long$delta, na.rm = TRUE)
  delta_lim <- c(min(yr[1], 0)*1.05, max(yr[2], 0)*1.05)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)

  scatter <- function(v) {
    sub <- long[variable == v]
    if (!inverted) {
      p <- ggplot(sub, aes(x = trait, y = delta)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.8) +
        geom_smooth(method = "loess", span = 0.75, se = TRUE,
                      colour = "grey20", fill = "grey70", alpha = 0.3) +
        scale_y_continuous(limits = delta_lim) +
        labs(x = if (show_x_axis) v else NULL,
              y = if (v == "LAI") method_yexpr(method) else NULL,
              title = if (show_col_titles) v else NULL)
    } else {
      p <- ggplot(sub, aes(x = delta, y = trait)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
        geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.8) +
        geom_smooth(method = "loess", span = 0.75, se = TRUE,
                      colour = "grey20", fill = "grey70", alpha = 0.3,
                      orientation = "y") +
        scale_x_continuous(limits = delta_lim) +
        labs(y = if (show_x_axis) v else NULL,
              x = if (v == "LAI") method_yexpr(method) else NULL,
              title = if (show_col_titles) v else NULL)
    }
    p +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
             legend.position = "none",
             axis.title = element_text(size = 14),
             axis.text = element_text(size = 12),
             panel.grid.minor = element_blank())
  }
  sub_lad <- long[variable == "LAD"]
  if (!inverted) {
    p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_violin(alpha = 0.45, colour = "grey30", scale = "width", trim = FALSE) +
      geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.5, size = 1.0) +
      geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                    fill = "white", outlier.shape = NA) +
      scale_y_continuous(limits = delta_lim) +
      labs(x = if (show_x_axis) "Cluster" else NULL, y = NULL,
            title = if (show_col_titles) "LAD" else NULL) +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
             legend.position = "none",
             axis.text.x = if (show_x_axis)
                              element_text(size = 12, face = "bold")
                            else element_blank(),
             axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
             axis.title.y = element_text(size = 14),
             axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_blank())
  } else {
    # Inverted : violin oriented horizontally, cluster on Y
    p_lad <- ggplot(sub_lad, aes(y = Cluster_lab, x = delta, fill = Cluster_lab)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_violin(alpha = 0.45, colour = "grey30", scale = "width", trim = FALSE) +
      geom_jitter(aes(colour = Cluster_lab), height = 0.18, alpha = 0.5, size = 1.0) +
      geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                    fill = "white", outlier.shape = NA) +
      scale_x_continuous(limits = delta_lim) +
      labs(y = if (show_x_axis) "Cluster" else NULL, x = NULL,
            title = if (show_col_titles) "LAD" else NULL) +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
             legend.position = "none",
             axis.text.y = if (show_x_axis)
                              element_text(size = 12, face = "bold")
                            else element_blank(),
             axis.title.x = element_blank(),
             axis.text.x = element_blank(), axis.ticks.x = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.y = element_blank())
  }
  p_lad <- p_lad + scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE)
  (scatter("LAI") + scatter("fCover") + scatter("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

build_typeB_v <- function(inverted = FALSE) {
  suffix <- if (inverted) "_inverted" else ""
  cli_h1(sprintf("Type B final%s", suffix))
  rows <- list(
    LOVB    = make_typeB_row("LOVB",    show_x_axis = FALSE, show_col_titles = TRUE,  inverted = inverted),
    LVA     = make_typeB_row("LVA",     show_x_axis = FALSE, show_col_titles = FALSE, inverted = inverted),
    Shapley = make_typeB_row("Shapley", show_x_axis = TRUE,  show_col_titles = FALSE, inverted = inverted)
  )
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_SHORT, levels=.CLUSTER_SHORT)),
            aes(x, x, colour = c)) + geom_point(size = 3.5) +
      scale_colour_manual(values = pal, name = NULL) +
      theme_bw(base_size = 15) +
      theme(legend.position = "bottom", legend.text = element_text(size = 13)))
  stacked <- (rows$LOVB / rows$LVA / rows$Shapley)
  final <- cowplot::plot_grid(stacked, leg, ncol = 1, rel_heights = c(1, 0.04))
  fig_path <- file.path(OUT, sprintf("fig_typeB_final%s.png", suffix))
  ggsave(fig_path, final, width = 14, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# (4) HOBO with larger fonts + inverted version
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

build_hobo_v <- function(inverted = FALSE) {
  suffix <- if (inverted) "_inverted" else ""
  cli_h1(sprintf("HOBO final%s", suffix))
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

  if (!inverted) {
    p <- ggplot(long, aes(x = dT, y = y))
    aes_lab <- list(x_lab = bquote(Delta * T[max] ~ "(°C)"),
                      y_lab = bquote("Simulated  " * Delta[v] ~ "(°C)"))
  } else {
    p <- ggplot(long, aes(x = y, y = dT))
    aes_lab <- list(x_lab = bquote("Simulated  " * Delta[v] ~ "(°C)"),
                      y_lab = bquote(Delta * T[max] ~ "(°C)"))
  }
  p <- p +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.7, alpha = 0.18, span = 0.9,
                  orientation = if (inverted) "y" else "x") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = rho_lab),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = p_lab),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 4.2,
               inherit.aes = FALSE, colour = "grey25") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = ci_lab),
               parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.7,
               inherit.aes = FALSE, colour = "grey35") +
    facet_grid(method ~ Variable, scales = "free", switch = "y") +
    labs(x = aes_lab$x_lab, y = aes_lab$y_lab) +
    theme_bw(base_size = 16) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 15),
           strip.placement = "outside",
           axis.title = element_text(size = 15),
           axis.text = element_text(size = 12),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, sprintf("fig_HOBO_final%s.png", suffix))
  ggsave(fig_path, p, width = 16, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_archetype_profiles()
  build_heatmap_square()
  build_typeB_v(inverted = FALSE)
  build_typeB_v(inverted = TRUE)
  build_hobo_v(inverted = FALSE)
  build_hobo_v(inverted = TRUE)
  cli_alert_success("All MEB v4 updates done")
}
