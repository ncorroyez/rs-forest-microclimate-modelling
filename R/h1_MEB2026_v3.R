# ==============================================================================
# MEB 2026 v3 : add Type A global + restore arch/cLHS n labels on heatmap.
# Outputs:
#   fig_typeA_final.png       : 3 rows × 4 cols, scatter Tmax_LOVB/LVA/Shapley vs Tmax_REF
#   fig_heatmap_MAE_final.png : with "n=4 archetypes" / "n=400 cLHS" tags
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
  "LOVB"    = bquote(T[max]^"LOVB,v" ~ "(°C)"),
  "LVA"     = bquote(T[max]^"LVA,v"  ~ "(°C)"),
  "Shapley" = bquote(T[max]^"REF" - varphi[v] ~ "(°C)"))

# ============================================================================
# Type A merged : 3 rows × 4 cols, scatter sim_method_v vs sim_REF
# ============================================================================
build_typeA_data <- function() {
  DT_c <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
  DT_phi <- readRDS(file.path(DATA, "DT_shapley_per_plot_floor05.rds"))
  phi_w <- dcast(DT_phi, x + y ~ variable, value.var = "phi", fun.aggregate = mean)
  setnames(phi_w, c("LAI","Hmax","fCover","LAD"),
            c("phi_LAI","phi_Hmax","phi_fCover","phi_LAD"))
  DT_c <- merge(DT_c, phi_w, by = c("x","y"), all.x = TRUE)
  for (v in VARS) DT_c[, paste0("Tmax_Shapley_", v) :=
                              get("Tmax_mean_REF") - get(paste0("phi_", v))]
  long <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m) {
    rbindlist(lapply(VARS, function(v) {
      ycol <- if (m == "LOVB")    paste0("Tmax_mean_LOVB_", v)
               else if (m == "LVA") paste0("Tmax_mean_LVA_", v)
               else                  paste0("Tmax_Shapley_", v)
      data.table(method = m, Variable = v,
                  Cluster = DT_c$Cluster,
                  Tmax_REF = DT_c$Tmax_mean_REF,
                  y = DT_c[[ycol]])
    }))
  }))
  long[, method := factor(method, levels = c("LOVB","LVA","Shapley"))]
  long[, Variable := factor(Variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_SHORT[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_SHORT)]
  long
}

build_typeA_final <- function() {
  cli_h1("Type A final (merged 3 rows × 4 cols)")
  long <- build_typeA_data()
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_SHORT)

  # Annotations : RMSE + MAE + R² per facet
  compute_stats <- function(y, x) {
    ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
    rmse <- sqrt(mean((y - x)^2))
    mae  <- mean(abs(y - x))
    ss_res <- sum((y - x)^2); ss_tot <- sum((y - mean(y))^2)
    R2 <- 1 - ss_res / ss_tot
    list(rmse = rmse, mae = mae, R2 = R2)
  }
  ann <- long[, compute_stats(y, Tmax_REF), by = .(method, Variable)]
  ann[, lab_rmse := sprintf("RMSE == %.2f", rmse)]
  ann[, lab_mae  := sprintf("MAE == %.2f", mae)]
  ann[, lab_r2   := sprintf("R^2 == %.2f", R2)]

  p <- ggplot(long, aes(x = Tmax_REF, y = y, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_rmse),
               parse = TRUE, hjust = -0.08, vjust = 1.3, size = 3.4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_mae),
               parse = TRUE, hjust = -0.08, vjust = 2.8, size = 3.4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab_r2),
               parse = TRUE, hjust = -0.08, vjust = 4.3, size = 3.4,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_grid(method ~ Variable, scales = "free", switch = "y") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(x = bquote(T[max]^"REF" ~ "(°C)"),
          y = bquote(T[max] ~ "simulated (method, v) (°C)")) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text = element_text(face = "bold", size = 13),
           strip.placement = "outside",
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, "fig_typeA_final.png")
  ggsave(fig_path, p, width = 16, height = 11, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
}

# ============================================================================
# Heatmap with n labels (smaller header)
# ============================================================================
build_heatmap_v3 <- function() {
  cli_h1("Heatmap MAE with n labels")
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
  # Small n labels at top
  top_labels <- data.table(x = c(2, 5, 7),
                              label = c("n = 4 archetypes",
                                          "n = 400 cLHS",
                                          ""))
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
    geom_text(data = top_labels[label != ""],
               aes(x = x, y = 4.7, label = label),
               inherit.aes = FALSE, fontface = "italic", size = 3.2,
               colour = "grey30") +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    coord_cartesian(ylim = c(0.5, 4.95), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(face = "bold", size = 11),
           axis.text.y = element_text(face = "bold", size = 13),
           panel.grid = element_blank(),
           plot.margin = margin(20, 10, 10, 10))
  fig_path <- file.path(OUT, "fig_heatmap_MAE_final.png")
  ggsave(fig_path, p, width = 11, height = 4.3, dpi = 150)
  cli_alert_success("Saved {.path {fig_path}}")
}

if (sys.nframe() == 0) {
  build_typeA_final()
  build_heatmap_v3()
  cli_alert_success("Type A + heatmap updated in {.path {OUT}}")
}
