# ==============================================================================
# Final figures for MEB 2026 (Silvilaser / Blois deck).
# Output : outputs/figs_MEB2026_final/
#   fig_typeB_{LOVB,LVA,Shapley}_final.png
#   fig_HOBO_{LOVB,LVA,Shapley}_final.png
#   fig_heatmap_MAE_final.png
#   tab_MAE_values.csv
# Conventions:
#   Δ_v positive = v contributes to buffering (positive when removing/adding v warms/cools sim)
#   x HOBO = Tmicro - Tmacro = obs - macro (negative when buffered)
#   MAE = mean(|.|) for all attribution ranks
# Dataset : floor05 (DT_*_floor05.rds)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
OUT  <- here::here("outputs/figs_MEB2026_final")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

VARS <- c("LAI","fCover","Hmax","LAD")
.CLUSTER_LABS <- c("1"="C1 Open","2"="C2 Dense","3"="C3 Low cover","4"="C4 Inter")
.PAL_CLUSTER <- c("C1 Open"="#E69F00","C2 Dense"="#0072B2",
                   "C3 Low cover"="#009E73","C4 Inter"="#CC79A7")

method_yexpr <- function(m) switch(m,
  "LOVB"    = bquote(Delta[v]^"LOVB"    ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA"     ~ "(°C)"),
  "Shapley" = bquote(varphi[v]          ~ "(°C)"))

# ===========================================================================
# Type B figures (4 panels: 3 scatter + LAD violin)
# ===========================================================================
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
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  long
}

build_typeB <- function(method) {
  cli_h2(sprintf("Type B %s", method))
  long <- build_long_typeB(method)
  yr <- range(long$delta, na.rm = TRUE)
  y_lim <- c(min(yr[1], 0)*1.05, max(yr[2], 0)*1.05)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  scatter <- function(v) {
    sub <- long[variable == v]
    ggplot(sub, aes(x=trait, y=delta)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
      geom_point(aes(colour=Cluster_lab), alpha=0.55, size=1.8) +
      geom_smooth(method="loess", span=0.75, se=TRUE,
                    colour="grey20", fill="grey70", alpha=0.3) +
      scale_colour_manual(values=pal, name=NULL, drop=FALSE) +
      scale_y_continuous(limits=y_lim) +
      labs(x = v,
            y = if (v == "LAI") method_yexpr(method) else NULL,
            title = paste("Contribution of", v)) +
      theme_bw(base_size=14) +
      theme(plot.title=element_text(face="bold", hjust=0.5, size=12),
             legend.position="none",
             axis.text.y = if (v=="LAI") element_text() else element_blank(),
             axis.ticks.y = if (v=="LAI") element_line() else element_blank(),
             panel.grid.minor=element_blank())
  }
  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x=Cluster_lab, y=delta, fill=Cluster_lab)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
    geom_violin(alpha=0.45, colour="grey30", scale="width", trim=FALSE) +
    geom_jitter(aes(colour=Cluster_lab), width=0.18, alpha=0.5, size=1.0) +
    geom_boxplot(width=0.12, alpha=0.7, colour="grey20",
                  fill="white", outlier.shape=NA) +
    scale_fill_manual(values=pal, name=NULL, drop=FALSE) +
    scale_colour_manual(values=pal, name=NULL, drop=FALSE) +
    scale_y_continuous(limits=y_lim) +
    labs(x = NULL, y = NULL,
          title = "Contribution of LAD (distribution)") +
    theme_bw(base_size=14) +
    theme(plot.title=element_text(face="bold", hjust=0.5, size=12),
           legend.position="none",
           axis.text.x=element_blank(), axis.ticks.x=element_blank(),
           axis.text.y=element_blank(), axis.ticks.y=element_blank(),
           panel.grid.minor=element_blank(),
           panel.grid.major.x=element_blank())

  row_plot <- (scatter("LAI") + scatter("fCover") + scatter("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))

  # Legend at bottom
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_LABS, levels=.CLUSTER_LABS)),
            aes(x, x, colour=c)) + geom_point(size=3) +
      scale_colour_manual(values=pal, name=NULL) +
      theme_bw(base_size=14) + theme(legend.position="bottom"))

  sub <- "x-axis = real value of v in the plot. LAD shown as distribution (functional object, no scalar projection). | LOESS span = 0.75 with 95%CI"
  title_full <- if (method == "Shapley") {
    bquote("Contribution"~varphi[v]~"vs trait (MEAN) — Shapley")
  } else {
    bquote("Contribution"~Delta[v]~"vs trait (MEAN) —"~.(method))
  }

  combined <- row_plot + plot_annotation(
    title = title_full,
    subtitle = sub,
    theme = theme(plot.title=element_text(face="bold", size=14),
                   plot.subtitle=element_text(colour="grey25", size=10)))
  final <- cowplot::plot_grid(combined, leg, ncol=1, rel_heights=c(1, 0.06))

  fig_path <- file.path(OUT, sprintf("fig_typeB_%s_final.png", method))
  # 4200 x 1650 px at 300 dpi = 14 x 5.5 inches
  ggsave(fig_path, final, width=14, height=5.5, dpi=300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(long)
}

# ===========================================================================
# HOBO figures (4 facets per method, scatter + Spearman)
# ===========================================================================
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

build_hobo <- function(method) {
  cli_h2(sprintf("HOBO %s", method))
  DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  # x = Tmicro - Tmacro = -dTmax_obs
  hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT_h <- merge(DT_h, hobo_t, by = "id_plot")

  # Build y per method
  for (v in VARS) {
    DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
    DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
  }
  if (method == "Shapley") {
    DT_shap <- readRDS(file.path(DATA, "DT_shapley_per_HOBO_floor05.rds"))
    shap_w <- dcast(DT_shap, id_plot ~ variable, value.var = "phi")
    setnames(shap_w, VARS, paste0("Shapley_", VARS))
    DT_h <- merge(DT_h, shap_w, by = "id_plot", all.x = TRUE)
  }

  vars_panel <- c("LAI","Hmax","fCover","LAD")
  long <- rbindlist(lapply(vars_panel, function(v) data.table(
    Variable = v, dT = DT_h$dT, y = DT_h[[paste0(method, "_", v)]])))
  long[, Variable := factor(Variable, levels = vars_panel)]
  long <- long[!is.na(y) & !is.na(dT)]

  ann <- long[, {r <- spearman_boot_ci(dT, y);
                  .(rho = r$rho, p = r$p, ci_lo = r$ci_lo, ci_hi = r$ci_hi)},
              by = Variable]
  ann[, rho_lab := sprintf("rho == %+.2f", rho)]
  ann[, p_lab   := ifelse(p < 1e-3, "italic(p) < 0.001",
                            sprintf("italic(p) == %.3f", p))]
  ann[, ci_lab  := sprintf("CI[95] *' '* '[' * %+.2f * '; ' * %+.2f * ']'",
                              ci_lo, ci_hi)]
  cli_alert("Spearman per variable (HOBO {method}, n = {length(unique(long$dT[long$Variable==vars_panel[1]]))}) :")
  print(ann[, .(Variable, rho = round(rho, 3), p = signif(p,3),
                  ci_lo = round(ci_lo,3), ci_hi = round(ci_hi,3))])

  p <- ggplot(long, aes(x = dT, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
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
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = paste0("HOBO validation (dTmax_obs) — ", method),
      subtitle = sprintf("n = %d sensors | Spearman ρ + bootstrap CI95 (1000 reps) | LOESS span = 0.9",
                          nrow(unique(long[, .(dT)]))),
      x = bquote(T[micro*","*max] - T[macro*","*max] ~ "(°C)"),
      y = method_yexpr(method),
      caption = "Non-parametric rank correlation. LOESS shown for visual reference only."
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=12),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           plot.caption = element_text(size=8, colour="grey35", hjust=0),
           panel.grid.minor = element_blank())
  fig_path <- file.path(OUT, sprintf("fig_HOBO_%s_final.png", method))
  # 3420 x 1290 px at 300 dpi = 11.4 x 4.3 inches
  ggsave(fig_path, p, width = 11.4, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(ann)
}

# ===========================================================================
# Heatmap MAE (rank ordering uniform across methods x scales)
# ===========================================================================
build_heatmap_MAE <- function() {
  cli_h2("Heatmap MAE")
  arch <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
  clhs <- fread(file.path(TAB, "tab_attribution_cLHS_mean_floor05.csv"))
  arch[, variable := factor(variable, levels = VARS)]
  clhs[, variable := factor(variable, levels = VARS)]

  # Build MAE table for export
  mae_table <- merge(arch[, .(variable,
                                  LOVB_arch    = LOVB_arch_value,
                                  LVA_arch     = LVA_arch_value,
                                  Shapley_arch = Shapley_arch_value)],
                       clhs[, .(variable,
                                  LOVB_cLHS    = LOVB_cLHS_value,
                                  LVA_cLHS     = LVA_cLHS_value,
                                  Shapley_cLHS = Shapley_cLHS_value)],
                       by = "variable")
  setorder(mae_table, variable)
  fwrite(mae_table, file.path(OUT, "tab_MAE_values.csv"))
  cli_alert("MAE values :"); print(mae_table)

  # Build long ranking table with order LOVB → LVA → Shapley
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
  cons[, method := "Consensus"]
  long_full <- rbind(long[, .(variable, method, col_idx, Rank)],
                      cons[, .(variable, method, col_idx, Rank)])
  long_full[, fill_rank := as.character(round(Rank))]

  pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
  x_labels <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley","Consensus")
  top_labels <- data.table(
    x = c(2, 5, 7),
    label = c("Level 1 — archetypes (n = 4 K-means centroids, 100 plots/cluster)",
                "Level 2 — cLHS (n = 400 plots)",
                "Consensus"))
  sep_xs <- c(3.5, 6.5)

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(method == "Consensus",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                                 "white", "grey15"),
               fontface = "bold", size = 5) +
    geom_vline(xintercept = sep_xs, colour = "white", linewidth = 3) +
    geom_text(data = top_labels,
               aes(x = x, y = 4.85, label = label),
               inherit.aes = FALSE, fontface = "bold", size = 2.9,
               colour = "grey20") +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
    labs(
      title    = "Attribution ranking (MAE) — 2 levels × 3 methods (LOVB → LVA → Shapley)",
      subtitle = "Rank based on MAE = mean(|Δ_v|) uniformly across all methods and scales",
      x = NULL, y = NULL,
      caption = "Dark = rank 1 (strongest)  |  light = rank 4 (weakest). n=400 plots (cLHS); n=4 centroids (archetypes, 100 plots/cluster)."
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 13),
           plot.subtitle = element_text(colour = "grey25", size = 10),
           plot.caption  = element_text(size = 7, colour = "grey35", hjust = 0),
           axis.text.x   = element_text(face = "bold", size = 9),
           axis.text.y   = element_text(face = "bold", size = 12),
           plot.margin   = margin(15, 10, 10, 10),
           panel.grid    = element_blank())

  fig_path <- file.path(OUT, "fig_heatmap_MAE_final.png")
  # 1642 x 748 px : use 11 x 5 in @ 150 dpi  ~  1650 x 750 px
  ggsave(fig_path, p, width = 11, height = 5, dpi = 150)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(list(long = long_full, mae = mae_table))
}

# ===========================================================================
# Main
# ===========================================================================
if (sys.nframe() == 0) {
  cli_h1("MEB 2026 final figures")
  build_typeB("LOVB")
  build_typeB("LVA")
  build_typeB("Shapley")
  build_hobo("LOVB")
  build_hobo("LVA")
  build_hobo("Shapley")
  build_heatmap_MAE()
  cli_alert_success("All 7 figures + 1 table in {.path {OUT}}")
}
