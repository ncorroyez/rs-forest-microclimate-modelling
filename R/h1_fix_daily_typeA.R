# Fix daily block figure + add R²/MAE alongside RMSE on the 4 Type A figures.
suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS         <- c("LAI", "fCover", "Hmax", "LAD")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")
.CLUSTER_ARCH <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Low cover", "Arch_C4" = "C4 Inter")

# ============================================================================
# Helper : compute metrics + 4 plotmath labels
# ============================================================================
compute_typeA_stats <- function(y, x) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  rmse <- sqrt(mean((y - x)^2))
  mae  <- mean(abs(y - x))
  bias <- mean(y - x)
  ss_res <- sum((y - x)^2); ss_tot <- sum((y - mean(y))^2)
  R2 <- 1 - ss_res / ss_tot
  r2_pearson <- cor(x, y)^2
  list(rmse = rmse, mae = mae, bias = bias, R2 = R2, r2_pearson = r2_pearson)
}
add_stat_labels <- function(ann_dt) {
  ann_dt[, lab1 := sprintf("RMSE == %.2f", rmse)]
  ann_dt[, lab2 := sprintf("MAE == %.2f", mae)]
  ann_dt[, lab3 := sprintf("R^2 == %.2f", R2)]
  ann_dt[, lab4 := sprintf("italic(r)^2 == %.2f", r2_pearson)]
}

# Layered geom_text calls (4 lines stacked)
add_4line_text <- function(ann_dt) {
  list(
    geom_text(data = ann_dt, aes(x = -Inf, y = Inf, label = lab1),
               parse = TRUE, hjust = -0.08, vjust = 1.4, size = 3.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15"),
    geom_text(data = ann_dt, aes(x = -Inf, y = Inf, label = lab2),
               parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15"),
    geom_text(data = ann_dt, aes(x = -Inf, y = Inf, label = lab3),
               parse = TRUE, hjust = -0.08, vjust = 4.6, size = 3.5,
               fontface = "bold", inherit.aes = FALSE, colour = "grey15"),
    geom_text(data = ann_dt, aes(x = -Inf, y = Inf, label = lab4),
               parse = TRUE, hjust = -0.08, vjust = 6.2, size = 3.5,
               inherit.aes = FALSE, colour = "grey35")
  )
}

# ============================================================================
# 1. Fix daily block figure (positive y range)
# ============================================================================
cli_h1("Fix daily block figure (y scale)")
tab_d <- fread(file.path(TAB, "tab_HOBO_spearman_daily_block.csv"))
tab_d[, Variable := factor(Variable, levels = VARS)]
p <- ggplot(tab_d, aes(x = Variable, y = rho_daily)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_pointrange(aes(ymin = ci_lo_block, ymax = ci_hi_block),
                    size = 0.8, colour = "#2C5F2D") +
  geom_text(aes(label = sprintf("%+.2f\n[%+.2f, %+.2f]",
                                  rho_daily, ci_lo_block, ci_hi_block),
                  y = ci_hi_block + 0.04),
             size = 3.8, fontface = "bold") +
  labs(
    title    = "Daily Spearman with block bootstrap by sensor — LOVB",
    subtitle = sprintf("n = %d daily pairs  |  %d sensors block-bootstrapped 1000 reps",
                        tab_d$n_pairs[1], tab_d$n_sensors[1]),
    x = NULL,
    y = bquote(rho[Spearman] ~ "daily  (" * Delta * T[max]^"obs"
                * ", " * Delta[v]^"LOVB"  * ")"),
    caption = "Positive ρ : MuSICA's daily attribution to v correlates with observed daily buffering. LAD borderline."
  ) +
  scale_y_continuous(limits = c(-0.1, 0.6),
                       breaks = seq(-0.1, 0.6, 0.1)) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face="bold", size=14),
         plot.subtitle = element_text(colour="grey25", size=11),
         axis.text.x = element_text(face="bold", size=12),
         panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_HOBO_spearman_daily_block.png"),
        p, width = 9, height = 5, dpi = 300)
cli_alert_success("Fixed fig_HOBO_spearman_daily_block.png")

# ============================================================================
# 2. Type A figures with R²/MAE/RMSE
# ============================================================================
build_typeA <- function(data_source, scope, variant) {
  cli_h2(sprintf("Type A %s %s", scope, variant))
  if (scope == "cLHS") {
    if (variant == "mean") {
      DT <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05.rds"))
      long <- rbindlist(lapply(VARS, function(v) data.table(
        group = DT$Cluster, Variable = v,
        Tmax_REF  = DT$Tmax_mean_REF,
        Tmax_LOVB = DT[[paste0("Tmax_mean_LOVB_", v)]])))
      long[, group_lab := .CLUSTER_LABS[as.character(group)]]
      long[, group_lab := factor(group_lab, levels = .CLUSTER_LABS)]
      pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
      subt <- "n = 400 plots     |     Dashed line = 1:1 identity"
      ttl  <- bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs"
                       ~ T[max]^"REF" ~ "(cLHS, MEAN over 122 d)")
      pt_alpha <- 0.55; pt_size <- 1.6
      x_lab <- bquote(T[max] ~ "buffering REF (°C)")
      y_lab <- bquote(T[max] ~ "buffering LOVB_v (°C)")
    } else {  # daily
      DT_d <- readRDS(file.path(DATA, "DT_daily_cLHS_floor05.rds"))
      ref <- DT_d[bit_code == "1111", .(x, y, date, Tmax_REF = Delta_Tmax)]
      bit_map <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110")
      ll <- list()
      for (v in VARS) {
        lv <- DT_d[bit_code == bit_map[v], .(x, y, date, Tmax_LOVB = Delta_Tmax)]
        m <- merge(ref, lv, by = c("x","y","date"))
        m[, Variable := v]
        ll[[v]] <- m
      }
      long <- rbindlist(ll)
      cl <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))[, .(x, y, group = Cluster)]
      long <- merge(long, cl, by = c("x","y"))
      long[, group_lab := .CLUSTER_LABS[as.character(group)]]
      long[, group_lab := factor(group_lab, levels = .CLUSTER_LABS)]
      pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)
      subt <- "n = 400 plots × 122 days     |     Dashed line = 1:1 identity"
      ttl  <- bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs"
                       ~ T[max]^"REF" ~ "(cLHS, DAILY)")
      pt_alpha <- 0.06; pt_size <- 0.4
      x_lab <- bquote(Delta * T[max] ~ "buffering REF, daily (°C)")
      y_lab <- bquote(Delta * T[max] ~ "buffering LOVB_v, daily (°C)")
    }
  } else {  # archetypes
    if (variant == "mean") {
      DT <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))
      long <- rbindlist(lapply(VARS, function(v) data.table(
        group = DT$archetype, Variable = v,
        Tmax_REF  = DT$Tmax_mean_REF,
        Tmax_LOVB = DT[[paste0("Tmax_mean_LOVB_", v)]])))
      long[, group_lab := .CLUSTER_ARCH[as.character(group)]]
      long[, group_lab := factor(group_lab, levels = .CLUSTER_ARCH)]
      pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)
      subt <- "n = 4 K-means centroids     |     Dashed line = 1:1 identity"
      ttl  <- bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs"
                       ~ T[max]^"REF" ~ "(archetypes, MEAN over 122 d)")
      pt_alpha <- 1.0; pt_size <- 4
      x_lab <- bquote(T[max] ~ "buffering REF (°C)")
      y_lab <- bquote(T[max] ~ "buffering LOVB_v (°C)")
    } else {
      DT_d <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
      ref <- DT_d[bit_code == "1111", .(archetype, date, Tmax_REF = Delta_Tmax)]
      bit_map <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110")
      ll <- list()
      for (v in VARS) {
        lv <- DT_d[bit_code == bit_map[v], .(archetype, date, Tmax_LOVB = Delta_Tmax)]
        m <- merge(ref, lv, by = c("archetype","date"))
        m[, Variable := v]
        ll[[v]] <- m
      }
      long <- rbindlist(ll)
      setnames(long, "archetype", "group")
      long[, group_lab := .CLUSTER_ARCH[as.character(group)]]
      long[, group_lab := factor(group_lab, levels = .CLUSTER_ARCH)]
      pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_ARCH)
      subt <- "n = 4 archetypes × 122 days = 488 daily pairs per panel"
      ttl  <- bquote("Type A — " ~ T[max]^"LOVB,v" ~ "vs"
                       ~ T[max]^"REF" ~ "(archetypes, DAILY)")
      pt_alpha <- 0.45; pt_size <- 1.2
      x_lab <- bquote(Delta * T[max] ~ "buffering REF, daily (°C)")
      y_lab <- bquote(Delta * T[max] ~ "buffering LOVB_v, daily (°C)")
    }
  }
  long[, Variable := factor(Variable, levels = VARS)]
  ann <- long[, c(compute_typeA_stats(Tmax_LOVB, Tmax_REF), .N),
                by = Variable]
  add_stat_labels(ann)

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = group_lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = pt_alpha, size = pt_size) +
    add_4line_text(ann) +
    facet_wrap(~ Variable, nrow = 1L, scales = "free") +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    labs(title = ttl, subtitle = subt, x = x_lab, y = y_lab) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=13),
           plot.title = element_text(face="bold", size=15),
           plot.subtitle = element_text(colour="grey25", size=11),
           legend.position = "bottom",
           panel.grid.minor = element_blank())
  if (variant == "daily" && scope == "cLHS")
    p <- p + guides(colour = guide_legend(override.aes = list(alpha=1, size=3)))

  fig_name <- sprintf("fig_typeA_%s_%s.png",
                         if (scope == "cLHS") "cLHS" else "archetypes", variant)
  ggsave(file.path(FIG, fig_name), p, width = 16, height = 5.5, dpi = 300)
  cli_alert_success("Saved {fig_name}")
}

build_typeA("cLHS", "cLHS", "mean")
build_typeA("cLHS", "cLHS", "daily")
build_typeA("arch", "archetypes", "mean")
build_typeA("arch", "archetypes", "daily")
