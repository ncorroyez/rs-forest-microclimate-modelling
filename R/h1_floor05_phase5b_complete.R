# ==============================================================================
# Phase 5b — Generate the complete floor05 figure suite using floor05 caches.
# All outputs in outputs/lovb_floor05/figures/ with PLAIN filenames (no suffix).
#
# Figures produced :
#   fig_typeB_facet.png                 (3 methods × 4 vars, LAD violin)
#   fig_typeB_{LOVB,LVA,Shapley}.png    (individual, LAD violin)
#   fig_typeA_facet.png                 (3 methods × 4 vars scatter Tmax_REF vs sim_v)
#   fig_archetypes_{LOVB,LVA,Shapley}_mean.png  (4 archetypes × 4 bars)
#   fig_HOBO_{LOVB,LVA,Shapley}_spearman_dTmax.png  (individual scatter facet)
#   fig_HOBO_jonckheere_bins_dTmax.png
#   fig_HOBO_spearman_daily_block.png   (visualisation block bootstrap journalier)
#   fig_heatmap_4col_simulation.png     (heatmap simplifiée 4 cols cLHS)
#   fig_HOBO_LOVB_typeA_validation.png  (Δ_v sim vs slope_obs scatter)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork); library(cowplot)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

VARS         <- c("LAI", "fCover", "Hmax", "LAD")
.PAL_VAR     <- c(LAI = "#440154", fCover = "#31688e",
                   Hmax = "#35b779", LAD = "#fde725")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")
.CLUSTER_ARCH <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Low cover", "Arch_C4" = "C4 Inter")

method_yexpr <- function(method) switch(method,
  "LOVB"    = bquote(Delta[v]^"LOVB" ~ "(°C)"),
  "LVA"     = bquote(Delta[v]^"LVA"  ~ "(°C)"),
  "Shapley" = bquote(Delta[v]^"Shapley" ~ "(°C)"))

spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 3L) return(list(rho = NA_real_, p = NA_real_, ci_lo = NA, ci_hi = NA))
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

# JT helpers
jt_test <- function(x, group, alternative = c("two.sided","increasing","decreasing")) {
  alternative <- match.arg(alternative)
  group <- as.factor(group); g_levels <- levels(group); k <- length(g_levels)
  J <- 0
  for (i in seq_len(k-1L)) for (j in seq.int(i+1L, k)) {
    x_i <- x[group == g_levels[i]]; x_j <- x[group == g_levels[j]]
    for (a in x_i) J <- J + sum(x_j > a) + 0.5 * sum(x_j == a)
  }
  n <- length(x); n_g <- tabulate(group, nbins = k)
  mu <- (n^2 - sum(n_g^2))/4
  sigma2 <- (n^2*(2*n+3) - sum(n_g^2*(2*n_g+3)))/72
  z <- (J - mu)/sqrt(sigma2)
  p <- switch(alternative,
                "two.sided"  = 2*(1 - pnorm(abs(z))),
                "increasing" = 1 - pnorm(z),
                "decreasing" = pnorm(z))
  list(z = z, p = p, statistic = J)
}
jt_boot_ci <- function(x, group, alternative, B = 1000L, seed = 42L) {
  set.seed(seed); n <- length(x); zb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    zb[b] <- tryCatch(jt_test(x[idx], group[idx], alternative)$z,
                       error = function(e) NA_real_)
  }
  quantile(zb, c(0.025, 0.975), na.rm = TRUE, type = 7)
}

# ============================================================================
# 1. Type B facet (3 methods × 4 vars, LAD violin)
# ============================================================================
build_typeB_data <- function(method) {
  DT_c <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
  trait_map <- c(LAI="LAI", fCover="fCover", Hmax="Hmax", LAD="FPC1")
  if (method == "LOVB") {
    rows <- lapply(VARS, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait = DT_c[[trait_map[v]]], delta = DT_c[[paste0("Delta_",v,"_mean")]]))
  } else if (method == "LVA") {
    rows <- lapply(VARS, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait = DT_c[[trait_map[v]]], delta = DT_c[[paste0("LVA_",v,"_mean")]]))
  } else if (method == "Shapley") {
    DT_phi <- readRDS(file.path(DATA, "DT_shapley_per_plot_floor05.rds"))
    DT_phi[, trait := NA_real_]
    for (v in VARS) DT_phi[variable == v, trait := get(trait_map[v])]
    rows <- list(DT_phi[, .(Cluster, variable, trait, delta = phi)])
  }
  long <- rbindlist(rows)
  long[, variable := factor(variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  long
}

make_typeB_row <- function(method, show_x = TRUE, show_titles = FALSE) {
  long <- build_typeB_data(method)
  yr <- range(long$delta, na.rm = TRUE)
  y_lim <- c(min(yr[1], 0)*1.05, max(yr[2], 0)*1.05)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  scatter <- function(v) {
    sub <- long[variable == v]
    ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.8) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      scale_y_continuous(limits = y_lim) +
      labs(x = if (show_x) v else NULL,
            y = if (v == "LAI") method_yexpr(method) else NULL,
            title = if (show_titles) v else NULL) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
             legend.position = "none",
             axis.text.y = if (v == "LAI") element_text() else element_blank(),
             axis.ticks.y = if (v == "LAI") element_line() else element_blank())
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
    scale_y_continuous(limits = y_lim) +
    labs(x = if (show_x) "Cluster" else NULL, y = NULL,
          title = if (show_titles) "LAD (distrib.)" else NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
           legend.position = "none",
           axis.text.x = if (show_x)
                            element_text(angle = 25, hjust = 1, size = 10)
                          else element_blank(),
           axis.ticks.x = if (show_x) element_line() else element_blank(),
           axis.text.y = element_blank(), axis.ticks.y = element_blank(),
           panel.grid.major.x = element_blank())
  (scatter("LAI") + scatter("fCover") + scatter("Hmax") + p_lad) +
    plot_layout(nrow = 1, widths = c(1,1,1,1))
}

build_typeB_facet <- function() {
  cli_h2("Type B facet (floor05)")
  rows <- list(
    LOVB    = make_typeB_row("LOVB",    show_x = FALSE, show_titles = TRUE),
    LVA     = make_typeB_row("LVA",     show_x = FALSE, show_titles = FALSE),
    Shapley = make_typeB_row("Shapley", show_x = TRUE,  show_titles = FALSE)
  )
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_LABS, levels=.CLUSTER_LABS)),
            aes(x, x, colour = c)) + geom_point(size=3) +
      scale_colour_manual(values=pal, name=NULL) +
      theme_bw(base_size=14) + theme(legend.position="bottom"))
  stacked <- (rows$LOVB / rows$LVA / rows$Shapley) +
    plot_annotation(
      title    = bquote("Type B — Per-plot attribution" ~ Delta[v]
                          ~ "vs LiDAR traits (n = 400 cLHS)"),
      subtitle = "Rows : LOVB / LVA / Shapley     |     Col 4 : LAD distribution by archetype",
      theme = theme(plot.title = element_text(face="bold", size=16),
                     plot.subtitle = element_text(colour="grey25", size=12)))
  final <- cowplot::plot_grid(stacked, leg, ncol = 1, rel_heights = c(1, 0.045))
  ggsave(file.path(FIG, "fig_typeB_facet.png"), final,
          width = 16, height = 14, dpi = 300)
  cli_alert_success("Saved fig_typeB_facet.png")
}

build_typeB_individual <- function(method) {
  cli_h2("Type B individual — {method}")
  row_plot <- make_typeB_row(method, show_x = TRUE, show_titles = TRUE)
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  leg <- cowplot::get_legend(
    ggplot(data.table(x=1:4, c=factor(.CLUSTER_LABS, levels=.CLUSTER_LABS)),
            aes(x, x, colour = c)) + geom_point(size=3) +
      scale_colour_manual(values=pal, name=NULL) +
      theme_bw(base_size=14) + theme(legend.position="bottom"))
  combined <- row_plot +
    plot_annotation(
      title    = bquote("Contribution" ~ Delta[v]
                          ~ "vs trait (MEAN) — " ~ .(method)),
      subtitle = "LAI / fCover / Hmax = scatter + LOESS span 0.75    |    LAD = distribution by archetype",
      theme = theme(plot.title = element_text(face="bold", size=15),
                     plot.subtitle = element_text(colour="grey25", size=11)))
  final <- cowplot::plot_grid(combined, leg, ncol = 1, rel_heights = c(1, 0.06))
  ggsave(file.path(FIG, sprintf("fig_typeB_%s.png", method)), final,
          width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_typeB_{method}.png")
}

# ============================================================================
# 2. Type A facet (3 methods × 4 vars scatter Tmax_REF vs sim_v)
# ============================================================================
build_typeA_facet <- function() {
  cli_h2("Type A scatter facet (floor05)")
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
                  Cluster = DT_c$Cluster, Tmax_REF = DT_c$Tmax_mean_REF,
                  y = DT_c[[ycol]])
    }))
  }))
  long[, method := factor(method, levels = c("LOVB","LVA","Shapley"))]
  long[, Variable := factor(Variable, levels = VARS)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
  ann <- long[, .(rmse = round(sqrt(mean((y - Tmax_REF)^2, na.rm = TRUE)), 2)),
                by = .(method, Variable)]
  ann[, label := paste0("RMSE == ", sprintf("%.2f", rmse), " * '°C'")]

  p <- ggplot(long, aes(x = Tmax_REF, y = y, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               parse = TRUE, hjust = -0.06, vjust = 1.4, size = 3.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
    facet_grid(method ~ Variable, scales = "free", switch = "y") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(
      title    = bquote("Type A — Tmax simulated (method) vs Tmax simulated (REF)"
                          ~ ", n = 400 cLHS"),
      subtitle = "Rows : LOVB / LVA / Shapley imputed     |     Dashed line = 1:1 identity     |     Deviation = v's contribution",
      x = bquote(T[max] ~ "buffering REF (mean over 122 d)" ~ "(°C)"),
      y = bquote(T[max] ~ "buffering method v" ~ "(°C)")
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 13),
           strip.placement  = "outside",
           plot.title       = element_text(face = "bold", size = 16),
           plot.subtitle    = element_text(colour = "grey25", size = 12),
           legend.position  = "bottom",
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_typeA_facet.png"), p,
          width = 16, height = 13, dpi = 300)
  cli_alert_success("Saved fig_typeA_facet.png")
}

# ============================================================================
# 3. HOBO Spearman individual figures (LOVB, LVA, Shapley)
# ============================================================================
build_HOBO_individual <- function(method) {
  cli_h2("HOBO Spearman individual — {method}")
  DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
  cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
  DT_h <- merge(DT_h, hobo, by = "id_plot")
  for (v in VARS) {
    DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
    DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
    DT_h[, paste0("Shapley_", v) := (get(paste0("LVA_", v)) + get(paste0("LOVB_", v)))/2]
  }
  long <- rbindlist(lapply(VARS, function(v) data.table(
    Variable = v, dT = DT_h$dT, y = DT_h[[paste0(method, "_", v)]])))
  long[, Variable := factor(Variable, levels = VARS)]
  ann <- long[, {r <- spearman_boot_ci(dT, y);
                  .(rho = r$rho, p = r$p, ci_lo = r$ci_lo, ci_hi = r$ci_hi)},
              by = Variable]
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
      title    = bquote("HOBO validation (n = 53) — " ~ .(method)
                          ~ "   " ~ rho[Spearman] ~ "vs" ~ Delta * T[max]),
      subtitle = "LOESS span = 0.9     |     Spearman ρ + bootstrap CI95 (1000 reps)",
      x = bquote(Delta * T[max]^"obs" ~ "(°C)"),
      y = method_yexpr(method)
    ) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 12),
           plot.title       = element_text(face = "bold", size = 15),
           plot.subtitle    = element_text(colour = "grey25", size = 11),
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, sprintf("fig_HOBO_%s_spearman_dTmax.png", method)),
          p, width = 14, height = 4.5, dpi = 300)
  cli_alert_success("Saved fig_HOBO_{method}_spearman_dTmax.png")
}

# ============================================================================
# 4. Jonckheere-Terpstra bins dTmax
# ============================================================================
build_jonckheere_bins <- function() {
  cli_h2("Jonckheere-Terpstra bins (dTmax, floor05)")
  DT_scal <- readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds"))
  cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
  DT <- merge(DT_scal, hobo, by = "id_plot")
  DT[, dT_signed := -dTmax_obs]   # Tmicro - Tmacro

  vars_jt <- c("LAI","Hmax","fCover","LAD_FPC1","LAD_L2")
  trait_cols <- c(LAI="LAI", Hmax="Hmax", fCover="fCover",
                   LAD_FPC1="FPC1", LAD_L2="LAD_L2")
  alt_lookup <- c(LAI="decreasing", Hmax="decreasing", fCover="decreasing",
                   LAD_FPC1="two.sided", LAD_L2="two.sided")
  rows <- list(); bins_long <- list()
  for (v in vars_jt) {
    trait <- DT[[trait_cols[v]]]; y <- DT$dT_signed
    keep <- !is.na(trait) & !is.na(y); trait <- trait[keep]; y <- y[keep]
    n <- length(trait); rnk <- rank(trait, ties.method = "first")
    bin <- cut(rnk, breaks = c(0, n/3, 2*n/3, n),
                include.lowest = TRUE, labels = c("low","mid","high"))
    bin <- factor(bin, levels = c("low","mid","high"))
    res <- jt_test(y, bin, alternative = alt_lookup[v])
    ci  <- jt_boot_ci(y, bin, alt_lookup[v], B = 1000L,
                        seed = 42L + which(vars_jt == v))
    rows[[length(rows)+1L]] <- data.table(Variable = v, JT_z = round(res$z, 3),
                                              JT_pvalue = signif(res$p, 4),
                                              JT_z_ci_lo = round(ci[1], 3),
                                              JT_z_ci_hi = round(ci[2], 3))
    bins_long[[v]] <- data.table(Variable = v, bin = bin, y = y,
                                    p_value = res$p, z = res$z,
                                    ci_lo = ci[1], ci_hi = ci[2])
  }
  tab <- rbindlist(rows)
  fwrite(tab, file.path(TAB, "tab_HOBO_jonckheere_dTmax.csv"))
  long <- rbindlist(bins_long)
  long[, Variable := factor(Variable, levels = vars_jt)]
  ann <- tab[, .(Variable, label = sprintf("JT p = %s\nz = %+.2f [%+.2f ; %+.2f]",
                                              ifelse(JT_pvalue < 1e-3, "<0.001",
                                                       sprintf("%.3f", JT_pvalue)),
                                              JT_z, JT_z_ci_lo, JT_z_ci_hi))]
  ann[, Variable := factor(Variable, levels = vars_jt)]
  p <- ggplot(long, aes(x = bin, y = y)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_boxplot(fill = "#2C5F2D", alpha = 0.55, colour = "grey20",
                  outlier.size = 1.3) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1, colour = "grey30") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.3, size = 3.0, fontface = "bold",
               inherit.aes = FALSE) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = bquote("Jonckheere-Terpstra trend test (" ~ Delta * T[max]^"obs" ~ ")"),
      subtitle = "Tertile binning by trait  |  n = 53 HOBO  |  JT z bootstrap CI95 (1000 reps)",
      x = NULL,
      y = bquote(Delta * T[max]^"obs" ~ "(°C)")
    ) +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold"),
           plot.title = element_text(face="bold", size=14),
           plot.subtitle = element_text(colour="grey25", size=11),
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_HOBO_jonckheere_bins_dTmax.png"),
          p, width = 13, height = 4.5, dpi = 300)
  cli_alert_success("Saved fig_HOBO_jonckheere_bins_dTmax.png")
  cli_alert("JT (floor05) :"); print(tab)
}

# ============================================================================
# 5. Daily block bootstrap (LOVB only, since LOVB is the natural daily method)
# ============================================================================
build_daily_block <- function() {
  cli_h2("Daily Spearman + block bootstrap by sensor (floor05)")
  # Observed daily dTmax (unchanged from original)
  obs <- readRDS(here::here("outputs/lovb/data/DT_HOBO_daily_obs.rds"))
  obs[, dT_signed := -dTmax_obs]
  sim <- as.data.table(readRDS(file.path(DATA, "DT_daily_HOBO_floor05.rds")))
  bit_map <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110", REF="1111")

  sim_ref <- sim[bit_code == bit_map[["REF"]],
                   .(id_plot, date, dTmax_REF = Delta_Tmax)]
  sim_w <- merge(obs[, .(id_plot, date, dT_signed)],
                   sim_ref, by = c("id_plot","date"))
  for (v in VARS) {
    sim_v <- sim[bit_code == bit_map[[v]],
                   .(id_plot, date, dTmax_LOVB = Delta_Tmax)]
    setnames(sim_v, "dTmax_LOVB", paste0("dTmax_LOVB_", v))
    sim_w <- merge(sim_w, sim_v, by = c("id_plot","date"))
  }
  for (v in VARS) sim_w[, paste0("Delta_",v,"_LOVB") :=
                              dTmax_REF - get(paste0("dTmax_LOVB_", v))]

  spear_pooled <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
  sensors <- unique(sim_w$id_plot)
  set.seed(42); res_daily <- list()
  for (v in VARS) {
    dcol <- paste0("Delta_",v,"_LOVB")
    rho_hat <- spear_pooled(sim_w$dT_signed, sim_w[[dcol]])
    rb <- numeric(1000L)
    for (b in seq_len(1000L)) {
      idx <- sample(sensors, length(sensors), replace = TRUE)
      sub <- sim_w[id_plot %in% idx]
      rb[b] <- spear_pooled(sub$dT_signed, sub[[dcol]])
    }
    ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
    res_daily[[v]] <- data.table(Variable = v, rho_daily = rho_hat,
                                    ci_lo_block = unname(ci[1]),
                                    ci_hi_block = unname(ci[2]),
                                    n_pairs = nrow(sim_w),
                                    n_sensors = length(sensors))
  }
  tab_d <- rbindlist(res_daily)
  tab_d[, Variable := factor(Variable, levels = VARS)]
  fwrite(tab_d, file.path(TAB, "tab_HOBO_spearman_daily_block.csv"))

  # Visual
  tab_d_long <- tab_d[, .(Variable, rho = rho_daily,
                            ci_lo = ci_lo_block, ci_hi = ci_hi_block)]
  p <- ggplot(tab_d_long, aes(x = Variable, y = rho)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.6,
                      colour = "#2C5F2D") +
    geom_text(aes(label = sprintf("%+.2f\n[%+.2f, %+.2f]", rho, ci_lo, ci_hi),
                    y = ci_hi + 0.05),
               size = 3.5, fontface = "bold") +
    labs(
      title    = bquote("Daily Spearman with block bootstrap by sensor — LOVB"),
      subtitle = sprintf("n = %d daily pairs  |  %d sensors block-bootstrapped 1000 reps",
                          tab_d$n_pairs[1], tab_d$n_sensors[1]),
      x = NULL,
      y = bquote(rho[Spearman] ~ "daily  (" * Delta * T[max]^"obs"
                  * ", " * Delta[v]^"LOVB"  * ")"),
      caption = "Block bootstrap by sensor handles temporal autocorrelation within sensors."
    ) +
    scale_y_continuous(limits = c(-0.4, 0.0)) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face="bold", size=14),
           plot.subtitle = element_text(colour="grey25", size=11),
           panel.grid.minor = element_blank())
  ggsave(file.path(FIG, "fig_HOBO_spearman_daily_block.png"),
          p, width = 9, height = 5, dpi = 300)
  cli_alert_success("Saved fig_HOBO_spearman_daily_block.png")
  cli_alert("Daily block bootstrap :"); print(tab_d)
}

# ============================================================================
# 6. Archetype bars (LOVB / LVA / Shapley)
# ============================================================================
build_arch_bars <- function() {
  cli_h2("Archetype bars (LOVB, LVA, Shapley) — floor05")
  DT_a <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
  arch_shap <- fread(file.path(TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
  # Need per-archetype Shapley for the bars — recompute from 16 coal daily
  DT_daily_a <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
  # Load missing 6 coalitions for C2/C3/C4 if not present
  bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
              "1000","1001","1010","1011","1100","1101","1110","1111")
  if (length(unique(DT_daily_a$bit_code)) < 16) {
    suppressMessages({library(ncdf4); library(tidyverse); library(musica.tools)})
    source(here::here("R/config.R"))
    source(here::here("R/musica.R"))
    source(here::here("R/io.R"))
    source(here::here("R/lovb_01_load.R"))
    MISSING <- c("0011","0101","0110","1001","1010","1100")
    df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
    extra <- list()
    for (arch in c("Arch_C2","Arch_C3","Arch_C4")) {
      dir_a <- here::here("out_files/H1_archetypes", arch)
      for (bit in MISSING) {
        nc <- file.path(dir_a, sprintf("musica_out_%s_ARCH_%s.nc",
                                          arch, lovb_bit_to_suffix(bit)))
        if (!file.exists(nc) || file.size(nc) < 1e6) next
        res <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq),
                         error = function(e) NULL)
        if (is.null(res) || nrow(res) == 0L) next
        dt <- as.data.table(res); dt[, archetype := arch]; dt[, bit_code := bit]
        extra[[length(extra)+1L]] <- dt
      }
    }
    DT_extra <- rbindlist(extra, fill = TRUE)[, .(archetype, date, bit_code, Delta_Tmax)]
    DT_daily_a <- rbind(DT_daily_a, DT_extra, fill = TRUE)
  }
  DT_agg <- DT_daily_a[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                        by = .(archetype, bit_code)]
  DT_w <- dcast(DT_agg, archetype ~ bit_code, value.var = "Tmax_mean")
  tcols <- paste0("Tmax_mean_", bits)
  setnames(DT_w, bits, tcols)

  shap_one <- function(tvec) {
    v_null <- as.numeric(tvec["0000"])
    v <- function(b) as.numeric(tvec[b]) - v_null
    mb <- function(bo, n=4) { b <- rep("0",n); b[bo] <- "1"; paste(b, collapse="") }
    out <- numeric(4L)
    for (vb in 1:4) {
      others <- setdiff(1:4, vb); total <- 0
      for (s in 0:3) {
        ss <- if (s==0) list(integer(0)) else combn(others, s, simplify = FALSE)
        w <- factorial(s) * factorial(4-s-1) / factorial(4)
        for (S in ss) total <- total + w * (v(mb(c(S, vb))) - v(mb(S)))
      }
      out[vb] <- total
    }
    setNames(out, c("LAI","Hmax","fCover","LAD"))
  }
  # Build per-archetype Shapley
  shap_rows <- list()
  for (i in seq_len(nrow(DT_w))) {
    arch <- DT_w$archetype[i]
    tvec <- as.numeric(unlist(DT_w[i, ..tcols])); names(tvec) <- bits
    if (any(is.na(tvec))) next
    phi <- shap_one(tvec)
    shap_rows[[arch]] <- data.table(archetype = arch,
                                       variable = names(phi),
                                       value = unname(phi))
  }
  shap_df <- rbindlist(shap_rows)

  # Build figure builder
  make_arch_bar <- function(long_df, method_name, ylab_expr) {
    long_df[, variable := factor(variable, levels = VARS)]
    long_df[, archetype_lab := .CLUSTER_ARCH[archetype]]
    long_df[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
    ymax <- max(abs(long_df$value), na.rm = TRUE) * 1.2
    ggplot(long_df, aes(x = variable, y = value, fill = variable)) +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_col(colour = "grey20", width = 0.7) +
      geom_text(aes(label = sprintf("%+.2f", value),
                      vjust = ifelse(value >= 0, -0.3, 1.3)),
                 fontface = "bold", size = 3.6) +
      facet_wrap(~ archetype_lab, nrow = 1L) +
      scale_fill_manual(values = .PAL_VAR, guide = "none") +
      scale_y_continuous(limits = c(-ymax, ymax),
                          expand = expansion(mult = c(0.06, 0.06))) +
      labs(title = method_name, subtitle = bquote("(MEAN)"),
            x = NULL, y = ylab_expr) +
      theme_bw(base_size = 14) +
      theme(strip.background = element_rect(fill="grey92", colour=NA),
             strip.text = element_text(face="bold", size=12),
             plot.title = element_text(face="bold", size=14))
  }

  # LOVB long
  lovb_long <- rbindlist(lapply(VARS, function(v) data.table(
    archetype = DT_a$archetype, variable = v,
    value = DT_a[[paste0("Delta_", v, "_mean")]])))
  ggsave(file.path(FIG, "fig_archetypes_LOVB_mean.png"),
          make_arch_bar(lovb_long, "LOVB contribution by archetype",
                          method_yexpr("LOVB")),
          width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_archetypes_LOVB_mean.png")

  # LVA
  lva_long <- rbindlist(lapply(VARS, function(v) data.table(
    archetype = DT_a$archetype, variable = v,
    value = DT_a[[paste0("LVA_", v, "_mean")]])))
  ggsave(file.path(FIG, "fig_archetypes_LVA_mean.png"),
          make_arch_bar(lva_long, "LVA contribution by archetype",
                          method_yexpr("LVA")),
          width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_archetypes_LVA_mean.png")

  # Shapley
  ggsave(file.path(FIG, "fig_archetypes_Shapley_mean.png"),
          make_arch_bar(shap_df, "Shapley contribution by archetype",
                          method_yexpr("Shapley")),
          width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved fig_archetypes_Shapley_mean.png")
}

# ============================================================================
# 7. Heatmap 4 cols simplified (cLHS only, 3 methods + consensus)
# ============================================================================
build_heatmap_4col <- function() {
  cli_h2("Heatmap 4 cols simplified (cLHS only)")
  DT <- fread(file.path(TAB, "tab_consensus_4methods_floor05.csv"))
  DT[, Variable := factor(Variable, levels = VARS)]
  methods <- c("Shapley_rank","LOVB_rank","LVA_rank")
  labels  <- c("Shapley","LOVB","LVA")
  setorder(DT, Consensus_rank)

  long <- melt(DT, id.vars = "Variable", measure.vars = methods,
                variable.name = "method_raw", value.name = "Rank")
  long[, method := factor(labels[match(method_raw, methods)],
                            levels = c(labels, "Consensus\n(mean rank)"))]
  long[, Rank_int := Rank]
  cons_long <- DT[, .(Variable = factor(Variable, levels = DT$Variable),
                        method = factor("Consensus\n(mean rank)",
                                          levels = c(labels, "Consensus\n(mean rank)")),
                        Rank = round(Consensus_mean, 1),
                        Rank_int = Consensus_rank)]
  all_long <- rbind(long[, .(Variable, method, Rank, Rank_int)], cons_long)
  all_long[, Variable := factor(Variable, levels = DT$Variable)]
  all_long[, fill_rank := as.character(round(Rank_int))]
  pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")

  p <- ggplot(all_long, aes(x = method, y = Variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(grepl("Consensus", method),
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(all_long$fill_rank) <= 2, "white", "grey15"),
               fontface = "bold", size = 6) +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_y_discrete(limits = rev) +
    labs(
      title = "Attribution ranking (4-method triangulation)",
      subtitle = "Dark = rank 1 (strongest)  |  light = rank 4 (weakest)",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25"),
           axis.text.x = element_text(face="bold", size=11),
           axis.text.y = element_text(face="bold", size=14),
           panel.grid = element_blank())
  ggsave(file.path(FIG, "fig_heatmap_4col_simulation.png"),
          p, width = 10, height = 5, dpi = 300)
  cli_alert_success("Saved fig_heatmap_4col_simulation.png")
}

# ---- Main entry point --------------------------------------------------------
if (sys.nframe() == 0) {
  cli_h1("Phase 5b — complete floor05 figure suite (no _floor05 suffix in names)")
  build_typeB_facet()
  for (m in c("LOVB","LVA","Shapley")) build_typeB_individual(m)
  build_typeA_facet()
  for (m in c("LOVB","LVA","Shapley")) build_HOBO_individual(m)
  build_jonckheere_bins()
  build_daily_block()
  build_arch_bars()
  build_heatmap_4col()
  cli_alert_success("ALL floor05 figures saved in {.path {FIG}}")
}
